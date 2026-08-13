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

echo Project update check completed.