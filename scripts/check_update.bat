@echo off

REM This file implement specific rules to keep compatibility between some updates

echo Checking project updates...

REM ========================== Rules from v1.4.1 to upper version ==================================
REM Files migration

IF EXIST ".\database\construction\neo4j_DB_construction.py" (
    IF EXIST ".\neo4j_DB_construction.py" (
        echo Removing obsolete file: .\neo4j_DB_construction.py
        del /F /Q ".\neo4j_DB_construction.py"
    )
)

IF EXIST ".\database\driver\neo4j_driver.py" (
    IF EXIST ".\neo4j_driver.py" (
        echo Removing obsolete file: .\neo4j_driver.py
        del /F /Q ".\neo4j_driver.py"
    )
)

IF EXIST ".\database\services\neo4j_requests.py" (
    IF EXIST ".\neo4j_requests.py" (
        echo Removing obsolete file: .\neo4j_requests.py
        del /F /Q ".\neo4j_requests.py"
    )
)

IF EXIST ".\utils\auth_utils.py" (
    IF EXIST ".\auth_utils.py" (
        echo Removing obsolete file: .\auth_utils.py
        del /F /Q ".\auth_utils.py"
    )
)

IF EXIST ".\docker_management\neo4j_container_management.py" (
    IF EXIST ".\neo4j_container_management.py" (
        echo Removing obsolete file: .\neo4j_container_management.py
        del /F /Q ".\neo4j_container_management.py"
    )
)

IF EXIST ".\docker_management\neo4j_available_docker_images_conf.py" (
    IF EXIST ".\neo4j_available_docker_images_conf.py" (
        echo Removing obsolete file: .\neo4j_available_docker_images_conf.py
        del /F /Q ".\neo4j_available_docker_images_conf.py"
    )
)


echo Project update check completed.