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
NEO4J_LOGS_DIR = os.path.abspath("./data/logs")
NEO4J_RUN_DIR = os.path.abspath("./data/run")


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

    docker_cmd = [
        "docker", "run", "--rm",
        #f"--cpus={MAX_CPU}",
        #f"--memory={MAX_MEM}",
        f"--memory-swap={MAX_SWAP}",
        #"-e", f"JAVA_OPTS=-Xmx{MAX_MEM} -Xms1g",
        #"-e", f"NEO4J_dbms.memory.heap.max_size={MAX_MEM}",
        "-v", f"{IMPORT_DIR}:/import",
        "-v", f"{NEO4J_BASE_DIR}/data:/data"
    ]
    # Linux/macOS :
    if hasattr(os, "getuid") and hasattr(os, "getgid"):
        docker_cmd.extend(["-u", f"{os.getuid()}:{os.getgid()}"])

    docker_cmd.extend([
        DOCKER_IMAGE,
        "neo4j-admin", "database", "load", "neo4j",
        "--from-path=/import", "--overwrite-destination=true"
    ])
    subprocess.run(docker_cmd, check=True)

@require_authorization
def import_csv(docker=True):
    READ_BUFFER_SIZE = get_conf_read_buffer_size()
    logger.info(f"📂 Importing CSV - data dir : {NEO4J_BASE_DIR}/data with read buffer size :  {READ_BUFFER_SIZE}...")

    if docker:
        prepare_data_directories_in_container()

        docker_cmd = [
            "docker", "run", "--rm",
            #f"--cpus={MAX_CPU}",
            "-v", f"{NEO4J_BASE_DIR}/data:/data",
            "-v", f"{IMPORT_DIR}:/import",
            "-e", f"NEO4J_AUTH={NEO4J_AUTH}"
        ]
        # Linux/macOS :
        if hasattr(os, "getuid") and hasattr(os, "getgid"):
            docker_cmd.extend(["-u", f"{os.getuid()}:{os.getgid()}"])

        docker_cmd.extend([
            DOCKER_IMAGE,
            "neo4j-admin", "database", "import", "full",
            "--verbose",
            f"--read-buffer-size={READ_BUFFER_SIZE}",
            #f"--max-off-heap-memory={MAX_MEM}",
            "--nodes=Sequence=/import/sequences.csv",
            "--nodes=Node=/import/nodes.csv",
            "--relationships=/import/relations.csv"
        ])
        subprocess.run(docker_cmd, check=True)
    else:
        data_dir = os.path.join(NEO4J_BASE_DIR, "data")
        neo4j_db_dir = os.path.join(data_dir, "databases", "neo4j")

        logger.info("🛠️ Preparing host directories for Apptainer...")

        os.makedirs(neo4j_db_dir, exist_ok=True)
        os.makedirs(IMPORT_DIR, exist_ok=True)
        os.makedirs(NEO4J_LOGS_DIR, exist_ok=True)
        os.makedirs(NEO4J_RUN_DIR, exist_ok=True)
        logger.info("🛠️ Directories ready on host filesystem.")

        logger.info(
            f"📂 Importing CSV - data dir : {NEO4J_BASE_DIR}/data with read buffer size :  {READ_BUFFER_SIZE}...")

        apptainer_cmd = [
            "apptainer", "exec",

            "--bind", f"{NEO4J_BASE_DIR}/data:/data",
            "--bind", f"{IMPORT_DIR}:/import",
            "--bind", f"{NEO4J_BASE_DIR}/logs:/var/lib/neo4j/logs",
            "--env", f"NEO4J_AUTH={NEO4J_AUTH}",
            f"docker://{DOCKER_IMAGE}",

            "neo4j-admin", "database", "import", "full",
            "--verbose",
            f"--read-buffer-size={READ_BUFFER_SIZE}",
            "--nodes=Sequence=/import/sequences.csv",
            "--nodes=Node=/import/nodes.csv",
            "--relationships=/import/relations.csv"
        ]
        subprocess.run(apptainer_cmd, check=True)

def start_container():
    if not os.path.exists(CONF_FILE):
        check_conf_file()

    if not os.path.exists(CONF_FILE):
        return False

    with open(CONF_FILE) as f:
        conf = json.load(f)

    if "container_name" not in conf:
        return False

    docker = conf.get("docker", True)

    container_name = str(conf["container_name"])

    if docker:
        result = subprocess.run(
            ["docker", "ps", "-a", "--format", "{{.Names}} {{.Status}}"],
            capture_output=True,
            text=True,
            check=True
        )

        lines = result.stdout.strip().splitlines()

        for line in lines:
            name, *status_parts = line.split()
            status = " ".join(status_parts)

            if name == container_name and status.startswith("Up"):
                return True

    else:
        result = subprocess.run(
            ["apptainer", "instance", "list"],
            capture_output=True,
            text=True,
            check=True
        )

        for line in result.stdout.splitlines():
            if container_name in line:
                return True

    remove_container(container_name, docker=docker)

    HTTP_PORT = int(conf["http_port"])
    BOLT_PORT = int(conf["bolt_port"])
    NEO4J_AUTH = conf["login"] + "/" + conf["password"]

    logger.info(
        f"🚀 Starting Neo4j {'docker container' if docker else 'apptainer instance'}..."
    )

    if docker:

        cmd = [
            "docker", "run", "-d",
            "--name", container_name,
            "-e", f"NEO4J_AUTH={NEO4J_AUTH}",
            "-e", "NEO4J_ACCEPT_LICENSE_AGREEMENT=yes",
            "-e", "NEO4J_apoc_export_file_enabled=true",
            "-e", "NEO4J_apoc_import_file_enabled=true",
            "-e", "NEO4J_apoc_import_file_use__neo4j__config=true",
            "-e", 'NEO4J_PLUGINS=[\"apoc\"]'
        ]
        #Only for linux mode
        if hasattr(os, "getuid") and hasattr(os, "getgid"):
            cmd.extend([
                "-u",
                f"{os.getuid()}:{os.getgid()}"
            ])

        cmd.extend([
            "-p", f"{HTTP_PORT}:7474",
            "-p", f"{BOLT_PORT}:7687",
            "-v", f"{NEO4J_BASE_DIR}/data:/data",
            "-v", f"{NEO4J_BASE_DIR}/logs:/logs",
            "-v", f"{NEO4J_BASE_DIR}/conf:/conf",
            "-v", f"{NEO4J_BASE_DIR}/import:/import",
            "-v", f"{NEO4J_BASE_DIR}/plugins:/plugins",
            DOCKER_IMAGE
        ])
        subprocess.run(cmd, check=True)
    else:
        cmd = [
            "apptainer",
            "exec",

            "--env", f"NEO4J_AUTH={NEO4J_AUTH}",
            "--env", "NEO4J_ACCEPT_LICENSE_AGREEMENT=yes",
            "--env", "NEO4J_apoc_export_file_enabled=true",
            "--env", "NEO4J_apoc_import_file_enabled=true",
            "--env", "NEO4J_apoc_import_file_use__neo4j__config=true",
            "--env", 'NEO4J_PLUGINS=["apoc"]',
            "--env", f"NEO4J_server_http_listen__address=0.0.0.0:{HTTP_PORT}",
            "--env", f"NEO4J_server_bolt_listen__address=0.0.0.0:{BOLT_PORT}",
            "--env", f"NEO4J_server_bolt_advertised__address=localhost:{BOLT_PORT}",

            #"--bind", f"{NEO4J_BASE_DIR}/data:/var/lib/neo4j/data",
            "--bind", f"{NEO4J_BASE_DIR}/data:/data",
            "--bind", f"{NEO4J_BASE_DIR}/logs:/var/lib/neo4j/logs",
            "--bind", f"{NEO4J_BASE_DIR}/conf:/conf",
            "--bind", f"{NEO4J_BASE_DIR}/import:/import",
            "--bind", f"{NEO4J_BASE_DIR}/run:/var/lib/neo4j/run",
            "--bind", f"{NEO4J_BASE_DIR}/plugins:/var/lib/neo4j/plugins",

            f"docker://{DOCKER_IMAGE}",
            "neo4j",
            "console"
        ]
        subprocess.Popen(cmd)


    time.sleep(10)

    logger.info(f"✅ Neo4j {container_name} is ready!")
    logger.info(f"🌍 HTTP: http://localhost:{HTTP_PORT}")
    logger.info(f"🔗 BOLT: bolt://localhost:{BOLT_PORT}")

    return True


@require_authorization
def write_config(container_name, HTTP_PORT=7474, BOLT_PORT=7687, docker=True):
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
        "log_level": "DEBUG",
        "docker":docker
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
        if docker == False:
            config.update({"docker": False})
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

    with open(CONF_FILE) as f:
        conf = json.load(f)

    docker = conf.get("docker", True)

    if container_name is None:
        if (
            "container_name" in conf
            and conf["container_name"]
        ):
            container_name = conf["container_name"]

    if container_name is None:
        logger.info("ℹ️ No container name provided. Nothing to stop.")
        return

    try:
        if docker:
            # ----------------------
            # DOCKER MODE
            # ----------------------
            result = subprocess.run(
                ["docker", "ps", "-a", "--format", "{{.Names}}"],
                capture_output=True,
                text=True,
                check=True
            )

            existing = result.stdout.strip().splitlines()

            if container_name in existing:
                logger.info(f"🛑 Stopping Docker container: {container_name}")
                subprocess.run(["docker", "stop", container_name], check=True)
                subprocess.run(["docker", "rm", container_name], check=True)
            else:
                logger.info(f"ℹ️ Docker container '{container_name}' does not exist.")

        else:
            # ----------------------
            # APPTAINER MODE
            # ----------------------
            result = subprocess.run(
                ["apptainer", "instance", "list"],
                capture_output=True,
                text=True,
                check=True
            )

            if container_name in result.stdout:
                logger.info(f"🛑 Stopping Apptainer instance: {container_name}")
                subprocess.run(
                    ["apptainer", "instance", "stop", container_name],
                    check=True
                )
            else:
                logger.info(f"ℹ️ Apptainer instance '{container_name}' does not exist.")

    except subprocess.CalledProcessError as e:
        logger.error(f"❌ Error while stopping container/instance: {e}")



"""
This function Stops and removes a Docker container or Apptainer instance if it exists.
Does nothing if it does not exist.
"""
def remove_container(container_name: str, docker=True):

    # docker = True
    # # --- Load config ---
    # if os.path.exists(CONF_FILE):
    #     with open(CONF_FILE) as f:
    #         conf = json.load(f)
    #     docker = conf.get("docker", True)

    try:
        # -----------------------
        # DOCKER MODE
        # -----------------------
        if docker:
            result = subprocess.run(
                ["docker", "ps", "-a", "--format", "{{.Names}}"],
                capture_output=True,
                text=True,
                check=True
            )
            existing_containers = result.stdout.strip().splitlines()

            if container_name in existing_containers:
                logger.info(f"🧹 Stopping and removing Docker container: {container_name}")
                subprocess.run(
                    ["docker", "rm", "-f", container_name],
                    check=True
                )
        # -----------------------
        # APPTAINER MODE
        # -----------------------
        else:
            result = subprocess.run(
                ["apptainer", "instance", "list"],
                capture_output=True,
                text=True,
                check=True
            )

            existing_instances = result.stdout.strip().splitlines()

            if any(container_name in line for line in existing_instances):
                logger.info(f"🧹 Stopping and removing Apptainer instance: {container_name}")
                subprocess.run(
                    ["apptainer", "instance", "stop", container_name],
                    check=True
                )

    except subprocess.CalledProcessError as e:
        logger.error(f"❌ Error while removing container/instance: {e}")



#This function create a docker neo4j database
#If a dump file exists in ./data/import directory it will load the data into database
#If no dump but nodes.csv, sequences.csv and relations.csv exists in ./data/import directory it will load these data into database
#In other cases it will create an empty database
@require_authorization
def create_db(container_name, docker_image=DOCKER_IMAGE, docker=True):
    creation_mode = "empty" 
    csv_import_mode = False
    
    if docker_image is not None :
        DOCKER_IMAGE = docker_image
    logger.info(f"🔧 Creating database for container '{container_name}' with {DOCKER_IMAGE} neo4j image")

    # Stop container
    remove_container(container_name, docker=docker)
    
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
        import_csv(docker=docker)
        csv_import_mode = True
    # --- New databse creation --- #
    else :
        remove_directories()
        # Directories creation
        for d in ["data", "logs", "conf", "import", "plugins", "gfa", "annotations"]:
            os.makedirs(os.path.join(NEO4J_BASE_DIR, d), exist_ok=True)
    
    
    # Save container conf in db_conj.json
    logger.info(f"writing conf for container name : {container_name}")
    write_config(container_name, docker=docker)
    
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
        docker_cmd = [
            "docker", "run", "--rm",
            #f"--cpus={MAX_CPU}",
            #f"--memory={MAX_MEM}",
            #f"--memory-swap={MAX_SWAP}",
            #"-e", f"JAVA_OPTS=-Xmx{MAX_MEM} -Xms1g",
            #"-e", f"NEO4J_dbms.memory.heap.max_size={MAX_MEM}",
            "-v", f"{NEO4J_BASE_DIR}/data:/data",
            "-v", f"{IMPORT_DIR}:/import"
        ]
        # Linux/macOS :
        if hasattr(os, "getuid") and hasattr(os, "getgid"):
            docker_cmd.extend(["-u", f"{os.getuid()}:{os.getgid()}"])

        docker_cmd.extend([
            DOCKER_IMAGE,
            "neo4j-admin", "database", "dump", "neo4j",
            f"--to-path=/import"
        ])
        subprocess.run(docker_cmd, check=True)




        logger.info(f"✅ Dump successfully created at: {os.path.join(IMPORT_DIR, 'neo4j.dump')}")
        start_container()

    except subprocess.CalledProcessError as e:
        logger.error(f"❌ Failed to create dump: {e}")