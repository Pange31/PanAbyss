#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 22 10:21:02 2025

@author: fgraziani
"""

import os
import shutil
import subprocess
import json
import time

from config import *
from auth_utils import require_authorization
import logging


logger = logging.getLogger("panabyss_logger")

# --- CONSTANTES ---
DOCKER_IMAGE = "neo4j:2025.05-community-bullseye"
NEO4J_BASE_DIR = os.path.abspath("./data")
CONF_FILE = os.path.abspath("./data/conf/neo4j.conf")
CONF_SOURCE_FILE = os.path.abspath("./install/conf/neo4j.conf")
CONF_FILE = os.path.abspath("./conf.json")
DOCKER_COMPOSE_CONF_PATH = os.path.abspath("./docker-compose.yml")
IMPORT_DIR = os.path.abspath("./data/import")
DUMP_FILE = os.path.join(IMPORT_DIR, "neo4j.dump")


#MAX_MEM = "24g"
MAX_SWAP = "25g"
#MAX_CPU = "8"

NEO4J_AUTH = "neo4j/Administrateur"
NEO4J_LOGIN = "neo4j"
NEO4J_PASSWORD = "Administrateur"

MAX_TIME_INDEX = 7200

@require_authorization
def prepare_data_directories_in_container():
    """
    Ensure /data/databases/neo4j exists inside the container (run as root to avoid permission issues).
    """
    logger.info("🛠️ Preparing data directories inside container...")
    subprocess.run([
        "docker", "run", "--rm",
        "--entrypoint", "mkdir",
        "--user=root",
        "-v", f"{NEO4J_BASE_DIR}/data:/data",
        DOCKER_IMAGE,
        "-p", "/data/databases/neo4j"
    ], check=True)
    logger.info("🛠️ ########## Preparing data directories inside container...")

@require_authorization
def remove_directories():
    for folder in ["data", "logs", "plugins"]:
        path = os.path.join(NEO4J_BASE_DIR, folder)
        if os.path.exists(path):
            shutil.rmtree(path)

@require_authorization
def import_dump():
    logger.info("📂 Importing dump...")
    prepare_data_directories_in_container()
    subprocess.run([
        "docker", "run", "--rm",
        #f"--cpus={MAX_CPU}",
        #f"--memory={MAX_MEM}",
        f"--memory-swap={MAX_SWAP}",
        #"-e", f"JAVA_OPTS=-Xmx{MAX_MEM} -Xms1g",
        #"-e", f"NEO4J_dbms.memory.heap.max_size={MAX_MEM}",
        "-v", f"{IMPORT_DIR}:/import",
        "-v", f"{NEO4J_BASE_DIR}/data:/data",
        "-u", f"{os.getuid()}:{os.getgid()}",
        DOCKER_IMAGE,
        "neo4j-admin", "database", "load", "neo4j",
        "--from-path=/import", "--overwrite-destination=true"
    ], check=True)


@require_authorization
def import_csv():
    READ_BUFFER_SIZE = get_conf_read_buffer_size()
    logger.info(f"📂 Importing CSV - data dir : {NEO4J_BASE_DIR}/data with read buffer size :  {READ_BUFFER_SIZE}...")
    prepare_data_directories_in_container()
    subprocess.run([
        "docker", "run", "--rm",
        #f"--cpus={MAX_CPU}",
        "-v", f"{NEO4J_BASE_DIR}/data:/data",
        "-v", f"{IMPORT_DIR}:/import",
        "-e", f"NEO4J_AUTH={NEO4J_AUTH}",
        "-u", f"{os.getuid()}:{os.getgid()}",
        DOCKER_IMAGE,
        "neo4j-admin", "database", "import", "full",
        "--verbose",
        f"--read-buffer-size={READ_BUFFER_SIZE}",
        #f"--max-off-heap-memory={MAX_MEM}",
        "--nodes=Sequence=/import/sequences.csv",
        "--nodes=Node=/import/nodes.csv",
        "--relationships=/import/relations.csv"
    ], check=True)


def start_container():
    if not os.path.exists(CONF_FILE):
        check_conf_file()
    if os.path.exists(CONF_FILE):
        with open(CONF_FILE) as f:
            conf = json.load(f)
        if "container_name" in conf :
            container_name = str(conf["container_name"])
            result = subprocess.run(
                ["docker", "ps", "-a", "--format", "{{.Names}} {{.Status}}"],
                capture_output=True, text=True, check=True
            )
            lines = result.stdout.strip().splitlines()

            for line in lines:
                name, *status_parts = line.split()
                status = " ".join(status_parts)

                if name == container_name:
                    if status.startswith("Up"):
                        return True

            remove_container(container_name)

            HTTP_PORT = int(conf["http_port"])
            BOLT_PORT = int(conf["bolt_port"])
            NEO4J_AUTH = conf["login"]+"/"+conf["password"]
            logger.info("🚀 Starting Neo4j container...")
            subprocess.run([
                "docker", "run", "-d",
                "--name", container_name,
                "-e", f"NEO4J_AUTH={NEO4J_AUTH}",
                "-e", "NEO4J_ACCEPT_LICENSE_AGREEMENT=yes",
                "-e", "NEO4J_apoc_export_file_enabled=true",
                "-e", "NEO4J_apoc_import_file_enabled=true",
                "-e", "NEO4J_apoc_import_file_use__neo4j__config=true",
                "-e", "NEO4J_PLUGINS=[\"apoc\"]",
                "-u", f"{os.getuid()}:{os.getgid()}",
                "-p", f"{HTTP_PORT}:7474",
                "-p", f"{BOLT_PORT}:7687",
                "-v", f"{NEO4J_BASE_DIR}/data:/data",
                "-v", f"{NEO4J_BASE_DIR}/logs:/logs",
                "-v", f"{NEO4J_BASE_DIR}/conf:/conf",
                "-v", f"{NEO4J_BASE_DIR}/import:/import",
                "-v", f"{NEO4J_BASE_DIR}/plugins:/plugins",
                DOCKER_IMAGE
            ], check=True)

            time.sleep(10)
            logger.info(f"✅ Neo4j {container_name} is ready!")
            logger.info(f"🌍 HTTP: http://localhost:{HTTP_PORT}")
            logger.info(f"🔗 BOLT: bolt://localhost:{BOLT_PORT}")
            return True

@require_authorization
def write_config(container_name, HTTP_PORT=7474, BOLT_PORT=7687):
    """
    Create or update the configuration file.
    If the file already exists, only update specific fields:
        - container_name
        - http_port
        - bolt_port
        - login
        - password
    Otherwise, create a new configuration file with default fields.
    """
    logger.info(f"Writing configuration for container: {container_name}")

    # Default configuration structure
    default_config = {
        "container_name": container_name,
        "http_port": HTTP_PORT,
        "bolt_port": BOLT_PORT,
        "login": NEO4J_LOGIN,
        "password": NEO4J_PASSWORD,
        "server_mode": False,
        "admin_mode": False,
        "admin_users": {
            "admin": "1234"
        },
        "server_log_mode": "both",
        "log_retention_days": 7,
        "log_level": "DEBUG"
    }

    # If the configuration file exists, load and update it
    if os.path.exists(CONF_FILE):
        logger.info(f"Configuration file {CONF_FILE} already exists. Updating specific fields...")
        try:
            with open(CONF_FILE, "r") as f:
                config = json.load(f)
        except json.JSONDecodeError:
            logger.warning(f"Existing config file {CONF_FILE} is invalid. Recreating it.")
            config = default_config.copy()
        except Exception as e:
            logger.error(f"Error reading config file: {e}")
            config = default_config.copy()

        # Update only specific fields
        if "http_port" not in config :
            config.update({"http_port": HTTP_PORT})
        if "bolt_port" not in config :
            config.update({"bolt_port": BOLT_PORT})
        if "login" not in config :
            config.update({"login": NEO4J_LOGIN})
        if "password" not in config :
            config.update({"password": NEO4J_PASSWORD})
        config.update({"container_name": container_name})

    else:
        logger.info(f"No existing configuration found. Creating a new one at {CONF_FILE}.")
        config = default_config

    # Write configuration back to file
    try:
        with open(CONF_FILE, "w") as f:
            json.dump(config, f, indent=2)
        logger.info("Configuration successfully written.")
    except Exception as e:
        logger.error(f"Failed to write configuration: {e}")
        
            


def stop_container(container_name=None):
    if not os.path.exists(CONF_FILE):
        check_conf_file()
    if container_name == None:
        with open(CONF_FILE) as f:
            conf = json.load(f)
            if "container_name" in conf and conf["container_name"] is not None and conf["container_name"] != "":
                container_name = conf["container_name"]
    if container_name is not None:
        """
        Stops Docker container if it exists.
        Does nothing if the container does not exist.
        """
        try:
            # List all container names (including stopped ones)
            result = subprocess.run(
                ["docker", "ps", "-a", "--format", "{{.Names}}"],
                capture_output=True, text=True, check=True
            )
            existing_containers = result.stdout.strip().splitlines()

            if container_name in existing_containers:
                logger.info(f"🛑 Stopping existing container: {container_name}")
                subprocess.run(["docker", "stop", container_name], check=True)
            else:
                logger.info(f"ℹ️ Container '{container_name}' does not exist. Nothing to do.")
        except subprocess.CalledProcessError as e:
            logger.error(f"❌ Error while stopping/removing container: {e}")



def remove_container(container_name: str):
    """
    Stops and removes a Docker container if it exists.
    Does nothing if the container does not exist.
    """
    try:
        # List all container names (including stopped ones)
        result = subprocess.run(
            ["docker", "ps", "-a", "--format", "{{.Names}}"],
            capture_output=True, text=True, check=True
        )
        existing_containers = result.stdout.strip().splitlines()

        if container_name in existing_containers:
            logger.info(f"Stopping and removing existing container: {container_name}")
            subprocess.run(["docker", "rm", "-f", container_name], check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"❌ Error while stopping/removing container: {e}")

#This function create a docker neo4j database
#If a dump file exists in ./data/import directory it will load the data into database
#If no dump but nodes.csv, sequences.csv and relations.csv exists in ./data/import directory it will load these data into database
#In other cases it will create an empty database
@require_authorization
def create_db(container_name, docker_image=DOCKER_IMAGE):
    creation_mode = "empty" 
    csv_import_mode = False
    
    if docker_image is not None :
        DOCKER_IMAGE = docker_image
    logger.info(f"🔧 Creating database for container '{container_name}' with {DOCKER_IMAGE} neo4j image")
    
    
    
    # Stop container
    remove_container(container_name)
    
    data_db_dir = os.path.join(NEO4J_BASE_DIR, "data", "databases", "neo4j")
    csv_nodes = os.path.join(IMPORT_DIR, "nodes.csv")
    csv_relations = os.path.join(IMPORT_DIR, "relations.csv")
    csv_sequences = os.path.join(IMPORT_DIR, "sequences.csv")
    
    if os.path.exists(data_db_dir) and os.listdir(data_db_dir):
        remove_directories()
    for d in ["data", "logs", "conf", "import", "plugins", "gfa", "annotations", "data/databases", "data/databases/neo4j"]:
        os.makedirs(os.path.join(NEO4J_BASE_DIR, d), exist_ok=True)
   
    # copy conf file
    if not os.path.isfile(CONF_FILE) and os.path.isfile(CONF_SOURCE_FILE):
        shutil.copy(CONF_SOURCE_FILE, os.path.join(NEO4J_BASE_DIR, "conf"))
        logger.info("🔧 Config file copied")
    else:
        if not os.path.isfile(CONF_SOURCE_FILE):
            logger.warning(f"⚠️ Config file {CONF_SOURCE_FILE} not found")
    
    # --- Import via dump ---
    if os.path.isfile(DUMP_FILE):
        logger.info("📥 Detected dump file for import")
        import_dump()
    
    # --- Import via CSV ---
    elif os.path.isfile(csv_nodes) and os.path.isfile(csv_relations) and os.path.isfile(csv_sequences):
        logger.info("📥 Detected CSV files for import")
        import_csv()
        csv_import_mode = True
    # --- New databse creation --- #
    else :
        remove_directories()
        # Directories creation
        for d in ["data", "logs", "conf", "import", "plugins", "gfa", "annotations"]:
            os.makedirs(os.path.join(NEO4J_BASE_DIR, d), exist_ok=True)
    
    
    # Save container conf in db_conj.json
    logger.info(f"writing conf for container name : {container_name}")
    write_config(container_name)
    
    # Launch the container
    start_container()

    if csv_import_mode:
        creation_mode = "csv"
    
    
    return creation_mode


@require_authorization    
def dump_db(container_name, docker_image=DOCKER_IMAGE):
    logger.info("📦 Creating Neo4j dump...")
    if docker_image is not None :
        DOCKER_IMAGE = docker_image
    # Stop container
    stop_container(container_name)

    try:
        subprocess.run([
            "docker", "run", "--rm",
            #f"--cpus={MAX_CPU}",
            #f"--memory={MAX_MEM}",
            #f"--memory-swap={MAX_SWAP}",
            #"-e", f"JAVA_OPTS=-Xmx{MAX_MEM} -Xms1g",
            #"-e", f"NEO4J_dbms.memory.heap.max_size={MAX_MEM}",
            "-v", f"{NEO4J_BASE_DIR}/data:/data",
            "-v", f"{IMPORT_DIR}:/import",
            "-u", f"{os.getuid()}:{os.getgid()}",
            DOCKER_IMAGE,
            "neo4j-admin", "database", "dump", "neo4j",
            f"--to-path=/import"
        ], check=True)

        logger.info(f"✅ Dump successfully created at: {os.path.join(IMPORT_DIR, 'neo4j.dump')}")
        start_container()

    except subprocess.CalledProcessError as e:
        logger.error(f"❌ Failed to create dump: {e}")