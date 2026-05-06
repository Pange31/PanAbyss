#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 30 16:45:39 2025

@author: fgraziani
"""



import os
import json
import shutil
import subprocess
import time
from neo4j import GraphDatabase
from neo4j.exceptions import ServiceUnavailable
import logging
import threading

logger = logging.getLogger("panabyss_logger")

DEFAULT_DB_LOAD_GFA_BATCH_SIZE = 2000000

#default maximum nodes number that can be visualized in the IHM
DEFAULT_MAX_NODES_TO_VISUALIZE = 10000

#default maximum nodes number that can be get from DB
DEFAULT_MAX_NODES_FROM_DB = 50000

#Default buffer size to load csv data in neo4j
# it limits the line's size (by default 4 * 1024 * 1024) :
# for long nodes sequences it could be necessary to improve this value
DEFAULT_MAX_READ_BUFFER_SIZE = 64000000

#Maximum number of gwas to store in sqlite database
DEFAULT_MAX_GWAS_STORE = 100
DEFAULT_MAX_GWAS_RUNNING_INACTIVITY_HOURS = 5
#Default max region is not limitated
DEFAULT_MAX_GWAS_REGIONS = 0
#Windows size to look for annotations before or after in gwas
DEFAULT_GWAS_ANNOTATIONS_WINDOWS_SIZE = 100000
#Max number of attempts to get annotations before or after in gwas
# => the total region over which annotations are searched = DEFAULT_GWAS_ANNOTATIONS_WINDOWS_SIZE x DEFAULT_GWAS_ANNOTATIONS_MAX_ATTEMPTS
DEFAULT_GWAS_ANNOTATIONS_MAX_ATTEMPTS = 5

CONF_FILE = os.path.abspath("./conf.json")
OLD_CONF_FILE = os.path.abspath("./db_conf.json")
INSTALL_CONF_FILE = os.path.abspath("./install/conf/conf.json")


#Delay in s to wait between tests of the database connexion
neo4j_start_delay = 10
#max_iter_test_connexion is the maximum number of pooling database to test if it is up
max_iter_test_connexion = 60

# #main driver
# _DRIVER = None

#Default value to add to the old config file
DEFAULT_ADDITIONS = {
    "server_mode": False,
    "admin_mode": False,
    "admin_users": {"admin": "1234"},
    "server_log_mode": "both",
    "log_retention_days": 7,
    "log_level": "DEBUG"
}

#Function to retrocompatibility with old conf file
def check_conf_file():
    """Check if the main config file exists; otherwise copy the old one,
    update it with default fields, and delete the old config file."""
    logger.debug("check conf file")
    # 1️⃣ If the main config file already exists → update it if necessary
    if os.path.exists(CONF_FILE):
        #Check that all the fiels in the install/conf file are present in the conf file (retrocompatibility)
        if os.path.exists(INSTALL_CONF_FILE):
            update_conf = False
            with open(INSTALL_CONF_FILE, "r", encoding="utf-8") as f1, open(CONF_FILE, "r", encoding="utf-8") as f2:
                install_conf = json.load(f1)
                user_conf = json.load(f2)
                for key, value in install_conf.items():
                    if key not in user_conf:
                        update_conf = True
                        user_conf[key] = value
            if update_conf:
                with open(CONF_FILE, "w", encoding="utf-8") as f:
                    json.dump(user_conf, f, indent=4, ensure_ascii=False)
        return
    
    # 2️⃣ If the old config file exists → copy and upgrade it
    if os.path.exists(OLD_CONF_FILE):
        logger.info(f"⚙️ Copying {OLD_CONF_FILE} to {CONF_FILE} ...")
        shutil.copyfile(OLD_CONF_FILE, CONF_FILE)
    
        # Try to load the JSON content
        try:
            with open(CONF_FILE, "r", encoding="utf-8") as f:
                conf_data = json.load(f)
        except json.JSONDecodeError:
            logger.error("❌ Error: old config file is invalid JSON. Creating a new one.")
            conf_data = {}
    
        # Add or update the required fields
        conf_data.update(DEFAULT_ADDITIONS)
    
        # Save the updated config
        with open(CONF_FILE, "w", encoding="utf-8") as f:
            json.dump(conf_data, f, indent=4, ensure_ascii=False)
        logger.info(f"✅ New config file created: {CONF_FILE}")
    
        # Delete the old config file
        try:
            os.remove(OLD_CONF_FILE)
            logger.info(f"🗑️  Old config file removed: {OLD_CONF_FILE}")
        except Exception as e:
            logger.error(f"⚠️ Could not remove old config file: {e}")
        return

#Function for config retrocompatibility
check_conf_file()

def get_conf(log_levels=["INFO", "DEBUG", "WARNING","ERROR", "CRITICAL", "NOTSET"]):
    with open(CONF_FILE) as f:
        conf = json.load(f)
    if "container_name" in conf:
        container_name = str(conf.get("container_name",""))
    else:
        container_name = None
    HTTP_PORT = int(conf.get("http_port", 7474))
    BOLT_PORT = int(conf.get("bolt_port",7687))
    AUTH = (conf.get("login", None),conf.get("password", None))
    DB_URL="bolt://localhost:"+str(BOLT_PORT)
    SERVER_MODE = bool(conf.get("server_mode",False))
    ADMIN_MODE = bool(conf.get("admin_mode",False))
    users_list=conf.get("admin_users", {})
    LOG_RETENTION_DAYS = int(conf.get("log_retention_days",7))
    SERVER_LOG_MODE = str(conf.get("server_log_mode","both"))
    LOG_LEVEL_PARAM = str(conf.get("log_level","INFO")).upper()
    GUNICORN_LOG_LEVEL = str(conf.get("gunicorn_log_level","")).upper()
    DB_LOAD_GFA_BATCH_SIZE = int(conf.get("db_gfa_loading_batch_size",DEFAULT_DB_LOAD_GFA_BATCH_SIZE))
    MAX_NODES_TO_VISUALIZE = int(conf.get("max_nodes_to_visualize", DEFAULT_MAX_NODES_TO_VISUALIZE))
    MAX_NODES_FROM_DB = int(conf.get("max_nodes_from_db", DEFAULT_MAX_NODES_FROM_DB))
    READ_BUFFER_SIZE = int(conf.get("read_buffer_size", DEFAULT_MAX_READ_BUFFER_SIZE))
    MAX_GWAS_STORE = int(conf.get("max_gwas_store", DEFAULT_MAX_GWAS_STORE))
    if MAX_GWAS_STORE <= 0:
        MAX_GWAS_STORE = None
    MAX_GWAS_RUNNING_INACTIVITY_HOURS = int(conf.get("max_gwas_running_inactivity_hours",DEFAULT_MAX_GWAS_RUNNING_INACTIVITY_HOURS))
    MAX_GWAS_REGIONS = int(conf.get("max_gwas_regions",DEFAULT_MAX_GWAS_REGIONS))
    if MAX_GWAS_REGIONS <= 0:
        MAX_GWAS_REGIONS = None
    GWAS_ANNOTATIONS_WINDOWS_SIZE = int(conf.get("gwas_annotations_windows_size",DEFAULT_GWAS_ANNOTATIONS_WINDOWS_SIZE))
    GWAS_ANNOTATIONS_MAX_ATTEMPTS = int(conf.get("gwas_annotations_max_attempts",DEFAULT_GWAS_ANNOTATIONS_MAX_ATTEMPTS))
    if LOG_LEVEL_PARAM in log_levels:
        LOG_LEVEL= "logging."+LOG_LEVEL_PARAM
    else:
        LOG_LEVEL= "logging.INFO"
    USERS = {}
    if len(users_list) > 0:
        for k, v in users_list.items():
            if not isinstance(k, str) or not isinstance(v, str):
                raise ValueError("Password must be string")
            if k.strip() == "" or v.strip() == "":
                raise ValueError("Usernames and password can't be empty")
            USERS[k.strip()] = v.strip()
    return {"container_name":container_name, "HTTP_PORT":HTTP_PORT, "BOLT_PORT":BOLT_PORT, 
            "AUTH":AUTH, "DB_URL":DB_URL, "SERVER_MODE":SERVER_MODE, "USERS":USERS, "ADMIN_MODE":ADMIN_MODE,
            "SERVER_LOG_MODE":SERVER_LOG_MODE,"LOG_RETENTION_DAYS":LOG_RETENTION_DAYS, "LOG_LEVEL":LOG_LEVEL,
            "GUNICORN_LOG_LEVEL":GUNICORN_LOG_LEVEL, "DB_LOAD_GFA_BATCH_SIZE":DB_LOAD_GFA_BATCH_SIZE,
            "MAX_NODES_TO_VISUALIZE":MAX_NODES_TO_VISUALIZE, "MAX_NODES_FROM_DB":MAX_NODES_FROM_DB, "READ_BUFFER_SIZE":READ_BUFFER_SIZE,
            "MAX_GWAS_STORE":MAX_GWAS_STORE, "MAX_GWAS_RUNNING_INACTIVITY_HOURS":MAX_GWAS_RUNNING_INACTIVITY_HOURS,
            "MAX_GWAS_REGIONS":MAX_GWAS_REGIONS, "GWAS_ANNOTATIONS_WINDOWS_SIZE":GWAS_ANNOTATIONS_WINDOWS_SIZE,
            "GWAS_ANNOTATIONS_MAX_ATTEMPTS":GWAS_ANNOTATIONS_MAX_ATTEMPTS}


CONF = get_conf()


def test_connection(DB_URL, AUTH):
    try:
        with GraphDatabase.driver(DB_URL, auth=AUTH) as driver:
            with driver.session() as session:
                session.run("RETURN 1")
        return True
    except ServiceUnavailable:
        return False

def load_config_from_json():
    if not os.path.exists(CONF_FILE):
        return None

    with open(CONF_FILE, "r") as f:
        return json.load(f)

def is_server_mode():
    if not os.path.exists(CONF_FILE):
        check_conf_file()
    server_mode = False
    if os.path.exists(CONF_FILE):
        conf=get_conf()
        server_mode = conf["SERVER_MODE"]
    return server_mode


def is_admin_mode():
    if not os.path.exists(CONF_FILE):
        check_conf_file()
    admin_mode = False
    if os.path.exists(CONF_FILE):
        conf=get_conf()
        admin_mode = conf["ADMIN_MODE"]
    return admin_mode

def get_conf_read_buffer_size():
    if CONF:
        return CONF["READ_BUFFER_SIZE"]
    else:
        return DEFAULT_MAX_READ_BUFFER_SIZE

def set_conf_value(key, value):
    with open(CONF_FILE) as f:
        conf = json.load(f)
        conf[key] = value
    with open(CONF_FILE, "w") as f:
        json.dump(conf, f, indent=4)

def get_db_load_gfa_batch_size():
    if CONF:
        return CONF["DB_LOAD_GFA_BATCH_SIZE"]
    else :
        return DEFAULT_DB_LOAD_GFA_BATCH_SIZE

def get_max_nodes_to_visualize():
    conf=get_conf()
    return conf["MAX_NODES_TO_VISUALIZE"]

def get_max_nodes_from_db():
    conf=get_conf()
    return conf["MAX_NODES_FROM_DB"]


def get_gwas_conf():
    if CONF:
        return CONF["MAX_GWAS_STORE"], CONF["MAX_GWAS_RUNNING_INACTIVITY_HOURS"], CONF["MAX_GWAS_REGIONS"], CONF["GWAS_ANNOTATIONS_WINDOWS_SIZE"], CONF["GWAS_ANNOTATIONS_MAX_ATTEMPTS"]
    else :
        return DEFAULT_MAX_GWAS_STORE, DEFAULT_MAX_GWAS_RUNNING_INACTIVITY_HOURS, DEFAULT_MAX_GWAS_REGIONS, DEFAULT_GWAS_ANNOTATIONS_WINDOWS_SIZE, DEFAULT_GWAS_ANNOTATIONS_MAX_ATTEMPTS

#Authorization is set to True for local installation of PanAbyss
#or if the mode admin is set to True
def check_authorization():
    authorization = True
    if is_server_mode() and not is_admin_mode():
        authorization = False
    return authorization
    
def get_users(): 
     if not os.path.exists(CONF_FILE):
        check_conf_file()
     users = {}
     if os.path.exists(CONF_FILE):
         conf=get_conf()
         users = conf["USERS"]
     return users

    
# def get_driver(max_retries=5, retry_delay=5):
#     if CONF:
#         if "container_name" in CONF and CONF["container_name"] != None and CONF["container_name"] != "":
#             DB_URL = CONF["DB_URL"]
#             AUTH= CONF["AUTH"]
#             for attempt in range(1, max_retries + 1):
#                 if test_connection(DB_URL, AUTH):
#                     driver = GraphDatabase.driver(DB_URL, auth=AUTH)
#                     return driver
#                 else:
#                     logger.warn(f"Retry {attempt}/{max_retries} to connect to driver")
#                     time.sleep(retry_delay)
#
#
#             if test_connection(DB_URL, AUTH):
#                     driver = GraphDatabase.driver(DB_URL, auth=AUTH)
#                     return driver
#             logger.error(f"❌ Fail to access Database of container {conf['container_name']}.")
#             return None
#         else:
#             return None
#     else:
#         return None