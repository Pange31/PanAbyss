@echo off
SETLOCAL ENABLEDELAYEDEXPANSION

:: --- Parameters ---
SET "ENV_NAME=panabyss"
SET "ENV_FILE=panabyss_win.yaml"
SET "HASH_FILE=.%ENV_NAME%_env_hash"
SET "DASH_PORT=8050"
IF NOT "%~1"=="" SET "DASH_PORT=%~1"


REM --- Check and apply project updates ---
call .\scripts\check_update.bat

if errorlevel 1 (
    echo Error while checking project updates.
    exit /b 1
)

:: --- Check if conda is available ---
where conda >nul 2>&1
IF ERRORLEVEL 1 (
    echo [ERROR] Conda not found in PATH.
    echo Please run this script from Anaconda Prompt or add conda to PATH.
    EXIT /B 1
)

:: --- Compute clean SHA1 hash ---
SET "CURRENT_HASH="
FOR /F "usebackq skip=1 tokens=*" %%H IN (certutil -hashfile "%ENV_FILE%" SHA1 ^| findstr /R /V "hash") DO (
    SET "LINE=%%H"
    SET "LINE=!LINE: =!"
    IF NOT "!LINE!"=="" (
        SET "CURRENT_HASH=!LINE!"
        GOTO done_hash
    )
)
:done_hash
SET "CURRENT_HASH=!CURRENT_HASH:~0,40!"
FOR /F "delims=" %%A IN ("!CURRENT_HASH!") DO SET "CURRENT_HASH=%%~A"

echo [DEBUG] Current hash: '!CURRENT_HASH!'

:: --- Check if environment exists ---
SET "ENV_EXISTS=0"
FOR /F "tokens=1" %%E IN ('conda env list ^| findstr /B "%ENV_NAME%"') DO (
    IF "%%E"=="%ENV_NAME%" SET "ENV_EXISTS=1"
)

:: --- Determine if update or creation needed ---
SET "UPDATE_ENV=0"
IF "!ENV_EXISTS!"=="0" (
    echo [INFO] Environment %ENV_NAME% not found. Creating it...
    call conda env create --name %ENV_NAME% --file "%ENV_FILE%"
    IF ERRORLEVEL 1 (
        echo [ERROR] Failed to create environment.
        EXIT /B 1
    )
    SET "UPDATE_ENV=1"
) ELSE (
    IF NOT EXIST "%HASH_FILE%" (
        echo [INFO] No previous hash found. Will update environment.
        SET "UPDATE_ENV=1"
    ) ELSE (
        SET /P PREV_HASH=<"%HASH_FILE%"
        SET "PREV_HASH=!PREV_HASH: =!"
        SET "PREV_HASH=!PREV_HASH:~0,40!"
        FOR /F "delims=" %%A IN ("!PREV_HASH!") DO SET "PREV_HASH=%%~A"
        echo [DEBUG] Previous hash: '!PREV_HASH!'
        IF /I NOT "!CURRENT_HASH!"=="!PREV_HASH!" (
            echo [INFO] YAML has changed. Updating environment...
            SET "UPDATE_ENV=1"
        ) ELSE (
            echo [INFO] Environment is up to date.
        )
    )
)

:: --- Update environment if needed ---
IF "!UPDATE_ENV!"=="1" (
    echo [INFO] Running conda env update...
    call conda env update --name %ENV_NAME% --file "%ENV_FILE%" --prune
    IF ERRORLEVEL 1 (
        echo [ERROR] Failed to update environment.
        EXIT /B 1
    )

>"%HASH_FILE%" echo !CURRENT_HASH!
)

:: --- Activate environment ---
echo [INFO] Activating environment %ENV_NAME%...
call conda activate %ENV_NAME%
IF ERRORLEVEL 1 (
    echo [ERROR] Failed to activate environment %ENV_NAME%.
    EXIT /B 1
)

REM Check if db_conf.json exists
if exist "db_conf.json" (
	REM file exists nothing to do
) else (
    REM Check if conf.json exists
    if exist "conf.json" (
		REM file exists nothing to do
    ) else (
        REM No existing conf file => copy default conf file
        echo Copy file conf.json to conf.json...
        copy "install\conf\conf.json" "conf.json"
        if %errorlevel% neq 0 (
            echo Error when copying conf file.
            exit /b 1
        )
    )
)


REM Check if neo4j.conf exists
if exist "data\conf\neo4j.conf" (
	REM file exists nothing to do
) else (

	REM No existing neo4j conf file => copy default conf file
	echo Copy file install\conf\neo4j.conf to data\conf\conf.json...
	copy "install\conf\neo4j.conf" "data\conf\neo4j.conf"
	if %errorlevel% neq 0 (
		echo Error when copying neo4j conf file.
		exit /b 1
	)

)

:: --- Launch Dash app ---
echo [INFO] Launching PanAbyss on port %DASH_PORT%...
python index.py --port %DASH_PORT%
IF ERRORLEVEL 1 (
    echo [ERROR] Dash app failed to start.
    EXIT /B 1
)

ENDLOCAL
