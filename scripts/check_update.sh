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

echo "Project update check completed."