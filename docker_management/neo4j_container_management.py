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
from utils.auth_utils import require_authorization
import logging
from pathlib import Path
import re


logger = logging.getLogger("panabyss_logger")

# --- CONSTANTES ---
DOCKER_IMAGE = "neo4j:2025.05-community-bullseye"
APP_DIR = Path(__file__).resolve().parents[1]
DATA_BASE_DIR = APP_DIR / "data"
NEO4J_BASE_DIR = APP_DIR / "data" / "database"
NEO4J_CONF_SOURCE_FILE = APP_DIR / "install" / "conf" / "neo4j.conf"
NEO4J_CONF_FILE = APP_DIR / "data" / "conf" / "neo4j.conf"
CONF_FILE = APP_DIR / "conf.json"
DOCKER_COMPOSE_FILE = APP_DIR / "docker_management" / "docker-compose.yml"
IMPORT_DIR = DATA_BASE_DIR / "import"
DUMP_FILE = IMPORT_DIR / "neo4j.dump"
NEO4J_LOGS_DIR = NEO4J_BASE_DIR / "logs"
NEO4J_PLUGINS_DIR = NEO4J_BASE_DIR / "plugins"
NEO4J_RUN_DIR = NEO4J_BASE_DIR / "run"



#MAX_MEM = "24g"
MAX_SWAP = "25g"
#MAX_CPU = "8"

NEO4J_AUTH = "neo4j/Administrateur"
NEO4J_LOGIN = "neo4j"
NEO4J_PASSWORD = "Administrateur"

MAX_TIME_INDEX = 7200


"""
This function will delete neo4j directories and configuration
"""
def reset_neo4j_data(delete_import_dir=False):

    if os.path.exists(NEO4J_BASE_DIR):
        shutil.rmtree(NEO4J_BASE_DIR)
    if delete_import_dir:
        if os.path.exists(IMPORT_DIR):
            shutil.rmtree(IMPORT_DIR)
    if os.path.exists(DOCKER_COMPOSE_FILE):
        os.remove(DOCKER_COMPOSE_FILE)



"""
This function creates the base neo4j directories
"""
def create_neo4j_base_dir():
    for d in ["database", "conf", "import", "gfa", "annotations"]:
        path = os.path.join(DATA_BASE_DIR, d)
        os.makedirs(path, exist_ok=True)
        os.chmod(path, 0o777)
    for d in ["data", "logs", "plugins"]:
        path = os.path.join(NEO4J_BASE_DIR, d)
        os.makedirs(path, exist_ok=True)
        os.chmod(path, 0o777)

@require_authorization
def import_dump():
    logger.info("📂 Importing dump...")

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


"""
This function uses the neo4j import procedure
"""
@require_authorization
def import_csv(docker=True):
    READ_BUFFER_SIZE = get_conf_read_buffer_size()
    logger.info(f"📂 Importing CSV - data dir : {DATA_BASE_DIR}/data with read buffer size :  {READ_BUFFER_SIZE}...")

    if docker:
        data_path = os.path.abspath(os.path.join(NEO4J_BASE_DIR, "data"))
        import_path = os.path.abspath(IMPORT_DIR)

        docker_cmd = [
            "docker", "run", "--rm",
            #f"--cpus={MAX_CPU}",
            "-v", f"{data_path}:/data",
            "-v", f"{import_path}:/import",
            #"-e", f"NEO4J_AUTH={NEO4J_AUTH}"
        ]
        # Linux/macOS :
        # if hasattr(os, "getuid") and hasattr(os, "getgid"):
        #     docker_cmd.extend(["-u", f"{os.getuid()}:{os.getgid()}"])

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
        data_dir = os.path.join(DATA_BASE_DIR, "../data")
        neo4j_db_dir = os.path.join(data_dir, "databases", "neo4j")

        logger.info("🛠️ Preparing host directories for Apptainer...")

        os.makedirs(NEO4J_BASE_DIR, mode=0o777, exist_ok=True)
        os.makedirs(IMPORT_DIR, mode=0o777, exist_ok=True)
        os.makedirs(NEO4J_LOGS_DIR, mode=0o777, exist_ok=True)
        os.makedirs(NEO4J_PLUGINS_DIR, mode=0o777, exist_ok=True)
        os.makedirs(NEO4J_RUN_DIR, mode=0o777, exist_ok=True)
        logger.info("🛠️ Directories ready on host filesystem.")

        logger.info(
            f"📂 Importing CSV - data dir : {DATA_BASE_DIR}/data with read buffer size :  {READ_BUFFER_SIZE}...")

        apptainer_cmd = [
            "apptainer", "exec",

            "--bind", f"{NEO4J_BASE_DIR}/data:/data",
            "--bind", f"{IMPORT_DIR}:/import",
            "--bind", f"{NEO4J_BASE_DIR}/logs:/var/lib/neo4j/logs",
            #"--env", f"NEO4J_AUTH={NEO4J_AUTH}",
            f"docker://{DOCKER_IMAGE}",

            "neo4j-admin", "database", "import", "full",
            "--verbose",
            f"--read-buffer-size={READ_BUFFER_SIZE}",
            "--nodes=Sequence=/import/sequences.csv",
            "--nodes=Node=/import/nodes.csv",
            "--relationships=/import/relations.csv"
        ]
        subprocess.run(apptainer_cmd, check=True)


"""
This function checks the existence of docker-compose.yml file
If not found it create it from conf data
"""
def create_docker_compose_file(
    container_name,
    http_port,
    bolt_port,
    neo4j_auth=None,
):
    compose_file = Path(DOCKER_COMPOSE_FILE)
    if compose_file.exists():
        return compose_file

    if neo4j_auth:
        auth = neo4j_auth
    else:
        auth = "none"
    remove_container(container_name, docker=True)
    compose_file.parent.mkdir(parents=True, exist_ok=True)

    compose_lines = [
        "services:",
        "  neo4j:",
        f"    container_name: {container_name}",
        f"    image: {DOCKER_IMAGE}",
    ]

    # if hasattr(os, "getuid") and hasattr(os, "getgid"):
    #     compose_lines.append(
    #         f'    user: "{os.getuid()}:{os.getgid()}"'
    #     )
    compose_lines.extend([
        "    environment:",
        f'      NEO4J_AUTH: "{auth}"',
        '      NEO4J_ACCEPT_LICENSE_AGREEMENT: "yes"',
        '      NEO4J_apoc_export_file_enabled: "true"',
        '      NEO4J_apoc_import_file_enabled: "true"',
        '      NEO4J_apoc_import_file_use__neo4j__config: "true"',
        """      NEO4J_PLUGINS: '["apoc"]'""",
        "    ports:",
        f'      - "{http_port}:7474"',
        f'      - "{bolt_port}:7687"',
        "    volumes:",
        '      - "../data/database/data:/data"',
        '      - "../data/database/logs:/logs"',
        '      - "../data/conf:/conf"',
        '      - "../data/import:/import"',
        '      - "../data/database/plugins:/plugins"',
    ])

    compose_content = "\n".join(compose_lines) + "\n"

    compose_file.write_text(
        compose_content,
        encoding="utf-8",
    )

    return compose_file


"""
This function checks if the container is already running
"""
def is_container_running(container_name: str, docker: bool = True) -> bool:
    if docker:
        result = subprocess.run(
            [
                "docker",
                "ps",
                "--filter",
                f"name=^{container_name}$",
                "--format",
                "{{.Names}}",
            ],
            capture_output=True,
            text=True,
            check=True,
        )

        return container_name in result.stdout.splitlines()

    else:
        result = subprocess.run(
            ["apptainer", "instance", "list"],
            capture_output=True,
            text=True,
            check=True,
        )

        return any(
            container_name in line
            for line in result.stdout.splitlines()
        )




def get_compose_project_name(container_name: str) -> str:
    project_name = container_name.lower()
    project_name = re.sub(r"[^a-z0-9_-]", "_", project_name)
    return project_name

"""
This function checks if the container is already running
If not, it launches the container
According to the configuration file, the container is launched on docker or on apptainer
"""
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
    HTTP_PORT = int(conf["http_port"])
    BOLT_PORT = int(conf["bolt_port"])
    NEO4J_AUTH = conf["login"] + "/" + conf["password"]

    if docker:
        # Make sure the docker-compose.yml exists
        compose_file = create_docker_compose_file(
            container_name,
            HTTP_PORT,
            BOLT_PORT
        )


    if is_container_running(container_name, docker):
        logger.info(f"✅ Neo4j container '{container_name}' is already running.")
        return True


    logger.info(
        f"🚀 Starting Neo4j "
        f"{'docker compose' if docker else 'apptainer instance'}..."
    )

    if docker:
        result = subprocess.run(
            [
                "docker",
                "compose",
                "-p",
                str(get_compose_project_name(container_name)),
                "-f",
                str(compose_file),
                "up",
                "-d",
            ],
            capture_output=True,
            text=True,
        )

        if result.returncode != 0:
            logger.error("❌ Docker Compose failed")
            logger.error(f"STDOUT:\n{result.stdout}")
            logger.error(f"STDERR:\n{result.stderr}")
            return False



    else:
        cmd = [
            "apptainer",
            "exec",

            #"--env", f"NEO4J_AUTH={NEO4J_AUTH}",
            "--env", "NEO4J_ACCEPT_LICENSE_AGREEMENT=yes",
            "--env", "NEO4J_apoc_export_file_enabled=true",
            "--env", "NEO4J_apoc_import_file_enabled=true",
            "--env", "NEO4J_apoc_import_file_use__neo4j__config=true",
            "--env", 'NEO4J_PLUGINS=["apoc"]',
            "--env", f"NEO4J_server_http_listen__address=0.0.0.0:{HTTP_PORT}",
            "--env", f"NEO4J_server_bolt_listen__address=0.0.0.0:{BOLT_PORT}",
            "--env", f"NEO4J_server_bolt_advertised__address=localhost:{BOLT_PORT}",

            "--bind", f"{NEO4J_BASE_DIR}/data:/var/lib/neo4j/data",
            "--bind", f"{NEO4J_BASE_DIR}/logs:/var/lib/neo4j/logs",
            "--bind", f"{DATA_BASE_DIR}/conf:/var/lib/neo4j/conf",
            "--bind", f"{DATA_BASE_DIR}/import:/import",
            "--bind", f"{NEO4J_RUN_DIR}:/var/lib/neo4j/run",
            "--bind", f"{NEO4J_BASE_DIR}/plugins:/var/lib/neo4j/plugins",

            f"docker://{DOCKER_IMAGE}",
            "neo4j",
            "console"
        ]

        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

        logger.info(f"Neo4j Apptainer PID: {process.pid}")

        time.sleep(2)

        if process.poll() is not None:
            stdout, stderr = process.communicate()

            logger.error(
                f"Neo4j exited immediately with code {process.returncode}"
            )
            logger.error(f"STDOUT:\n{stdout}")
            logger.error(f"STDERR:\n{stderr}")

            return False

    time.sleep(10)

    logger.info(f"✅ Neo4j {container_name} is ready!")
    logger.info(f"🌍 HTTP: http://localhost:{HTTP_PORT}")
    logger.info(f"🔗 BOLT: bolt://localhost:{BOLT_PORT}")

    return True

"""
Create or update the configuration file.
If the file already exists, only update specific fields:
    - container_name
    - http_port
    - bolt_port
    - login
    - password
Otherwise, create a new configuration file with default fields.
Create the docker compose file.
"""
@require_authorization
def write_config(container_name, HTTP_PORT=7474, BOLT_PORT=7687, docker=True):
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
        NEO4J_AUTH = NEO4J_LOGIN + "/" + NEO4J_PASSWORD
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
        NEO4J_AUTH = f"{config['login']}/{config['password']}"
    #Create docker compose file
    if docker :
        compose_file = create_docker_compose_file(
            container_name,
            HTTP_PORT,
            BOLT_PORT
        )

    # Write configuration back to file
    try:
        with open(CONF_FILE, "w") as f:
            json.dump(config, f, indent=2)
        logger.info("Configuration successfully written.")
    except Exception as e:
        logger.error(f"Failed to write configuration: {e}")
        

"""
This function stops the container (docker or apptainer according to the conf file)
"""
def stop_container(container_name=None):
    if not os.path.exists(CONF_FILE):
        check_conf_file()

    with open(CONF_FILE) as f:
        conf = json.load(f)

    docker = conf.get("docker", True)

    if container_name is None:
        if "container_name" in conf and conf["container_name"]:
            container_name = conf["container_name"]

    if container_name is None:
        logger.info("ℹ️ No container name provided. Nothing to stop.")
        return

    try:
        if docker:
            # ----------------------
            # DOCKER MODE
            # ----------------------
            compose_file = Path(DOCKER_COMPOSE_FILE)

            if not compose_file.exists():
                logger.info(
                    f"ℹ️ Docker Compose file '{compose_file}' does not exist."
                )
                return

            if is_container_running(container_name, docker):
                logger.info(f"🛑 Stopping Docker Compose service: {container_name}")
                result = subprocess.run(
                    [
                        "docker",
                        "compose",
                        "-p",
                        str(get_compose_project_name(container_name)),
                        "-f",
                        str(compose_file),
                        "stop",
                    ],
                    capture_output=True,
                    text=True,
                )

                if result.returncode != 0:
                    logger.error("❌ Docker Compose failed")
                    logger.error(f"STDOUT:\n{result.stdout}")
                    logger.error(f"STDERR:\n{result.stderr}")
                    return False
            else:
                logger.info(
                    f"ℹ️ Docker Compose service '{container_name}' "
                    f"is not running."
                )

        else:
            # ----------------------
            # APPTAINER MODE
            # ----------------------
            if is_container_running(container_name, docker):
                logger.info(f"🛑 Stopping Apptainer instance: {container_name}")

                subprocess.run(
                    [
                        "apptainer",
                        "instance",
                        "stop",
                        container_name,
                    ],
                    check=True
                )
            else:
                logger.info(
                    f"ℹ️ Apptainer instance "
                    f"'{container_name}' does not exist."
                )

    except subprocess.CalledProcessError as e:
        logger.error(
            f"❌ Error while stopping container/instance: {e}"
        )



"""
This function Stops and removes a Docker container or Apptainer instance if it exists.
Does nothing if it does not exist.
"""
def remove_container(container_name: str, docker=True):
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


"""
This function creates a docker neo4j database
If a dump file exists in ./data/import directory it will load the data into database
If no dump but nodes.csv, sequences.csv and relations.csv exists in ./data/import directory it will load these data into database
In other cases it will create an empty database
"""
@require_authorization
def create_db(container_name, docker_image=DOCKER_IMAGE, docker=True):
    creation_mode = "empty" 
    csv_import_mode = False
    
    if docker_image is not None :
        DOCKER_IMAGE = docker_image
    logger.info(f"🔧 Creating database for container '{container_name}' with {DOCKER_IMAGE} neo4j image")
    # Stop container
    remove_container(container_name, docker=docker)
    
    #data_db_dir = os.path.join(DATA_BASE_DIR, "../data", "databases", "neo4j")
    csv_nodes = os.path.join(IMPORT_DIR, "nodes.csv")
    csv_relations = os.path.join(IMPORT_DIR, "relations.csv")
    csv_sequences = os.path.join(IMPORT_DIR, "sequences.csv")
    
    #if os.path.exists(data_db_dir) and os.listdir(data_db_dir):
    reset_neo4j_data()
    create_neo4j_base_dir()

    # copy conf file
    if not os.path.isfile(NEO4J_CONF_FILE) and os.path.isfile(NEO4J_CONF_SOURCE_FILE):
        shutil.copy(NEO4J_CONF_SOURCE_FILE, os.path.join(DATA_BASE_DIR, "conf"))
        logger.info("🔧 Config file copied")
    else:
        if not os.path.isfile(NEO4J_CONF_SOURCE_FILE):
            logger.warning(f"⚠️ Config file {NEO4J_CONF_SOURCE_FILE} not found")
    
    # --- Import via dump ---
    if os.path.isfile(DUMP_FILE):
        logger.info("📥 Detected dump file for import")
        import_dump()
    
    # --- Import via CSV ---
    elif os.path.isfile(csv_nodes) and os.path.isfile(csv_relations) and os.path.isfile(csv_sequences):
        logger.info("📥 Detected CSV files for import")
        import_csv(docker=docker)
        csv_import_mode = True
    
    
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