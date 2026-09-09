#!/bin/bash
# This file implement specific rules to keep compatibility between some updates
set -e

echo "Checking project updates..."

########### Rules from v1.4.1 to upper version #######################
# neo4j_db_construction.py
if [ -f "./database/construction/neo4j_DB_construction.py" ]; then
    if [ -f "./neo4j_DB_construction.py" ]; then
        echo "Removing obsolete file: ./neo4j_DB_construction.py"
        rm "./neo4j_DB_construction.py"
    fi
fi

# neo4j_driver.py
if [ -f "./database/driver/neo4j_driver.py" ]; then
    if [ -f "./neo4j_driver.py" ]; then
        echo "Removing obsolete file: ./neo4j_driver.py"
        rm "./neo4j_driver.py"
    fi
fi

# neo4j_requetes.py
if [ -f "./database/services/neo4j_requests.py" ]; then
    if [ -f "./neo4j_requests.py" ]; then
        echo "Removing obsolete file: ./neo4j_requests.py"
        rm "./neo4j_requests.py"
    fi
fi

# auth_utils.py
if [ -f "./utils/auth_utils.py" ]; then
    if [ -f "./auth_utils.py" ]; then
        echo "Removing obsolete file: ./auth_utils.py"
        rm "./auth_utils.py"
    fi
fi

#neo4j_container_management.py
if [ -f "./neo4j_management/neo4j_container_management.py" ]; then
    if [ -f "./neo4j_container_management.py" ]; then
        echo "Removing obsolete file: ./neo4j_container_management.py"
        rm "./neo4j_container_management.py"
    fi
fi

#neo4j_available_docker_images_conf.py
if [ -f "./neo4j_management/neo4j_available_docker_images_conf.py" ]; then
    if [ -f "./neo4j_available_docker_images_conf.py" ]; then
        echo "Removing obsolete file: ./neo4j_available_docker_images_conf.py"
        rm "./neo4j_available_docker_images_conf.py"
    fi
fi

echo "Project update check completed."