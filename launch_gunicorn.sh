#!/bin/bash
# --- Parameters ---
DASH_PORT=${1:-8050}
ENV_NAME="panabyss"
ENV_FILE="panabyss.yaml"
HASH_FILE="./.${ENV_NAME}_env_hash"
CONF_FILE="./conf.json"

set -e

# --- Initialize conda ---
CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"

# --- Ensure mamba is installed in base ---
if ! command -v mamba &> /dev/null; then
    echo "Mamba not found, installing in base environment..."
    conda install -y -n base mamba --channel conda-forge
fi

# --- Compute current YAML hash ---
CURRENT_HASH=$(sha1sum "$ENV_FILE" | awk '{print $1}')
PREV_HASH=""
if [ -f "$HASH_FILE" ]; then
    PREV_HASH=$(cat "$HASH_FILE")
fi

# --- Create or update environment with mamba ---
if conda info --envs | grep -q "^$ENV_NAME\s"; then
    if [ "$CURRENT_HASH" != "$PREV_HASH" ]; then
        echo "Environment '$ENV_NAME' exists but YAML changed. Updating..."
        mamba env update --name "$ENV_NAME" --file "$ENV_FILE" --prune
        echo "$CURRENT_HASH" > "$HASH_FILE"
    else
        echo "Environment '$ENV_NAME' exists and is up to date."
    fi
else
    echo "Environment '$ENV_NAME' not found. Creating..."
    mamba env create --file "$ENV_FILE"
    echo "$CURRENT_HASH" > "$HASH_FILE"
fi

# Activate environment
echo "Activating environment '$ENV_NAME'"
conda activate "$ENV_NAME"

# --- Config checks ---
# Check conf file
#First check old config file
if [ -f "db_conf.json" ]; then
	:
else
    # Check if conf.json file exists
    if [ -f "conf.json" ]; then
		:
    else
        # Old conf file and new config file don't exist => create default conf file
        echo "Copy conf file install/conf/conf.json to conf.json..."
        cp "install/conf/conf.json" "conf.json"
        if [ $? -ne 0 ]; then
            echo "Error when copying conf file."
            exit 1
        fi
    fi
fi

# Check neo4j conf file
if [ ! -f "data/conf/neo4j.conf" ]; then
    # neo4j conf file doesn't exists => create default neo4j conf file
	echo "Copy neo4j conf file install/conf/neo4j.conf to data/conf/neo4j.conf..."
	cp "install/conf/neo4j.conf" "data/conf/neo4j.conf"
	if [ $? -ne 0 ]; then
		echo "Error when copying neo4j conf file."
		exit 1
	fi
fi

echo $DASH_PORT > gunicorn.port

# --- Launch application ---
ENV_GUNICORN=$(conda info --base)/envs/$ENV_NAME/bin/gunicorn
echo "Launching PanAbyss on port $DASH_PORT with command $ENV_GUNICORN..."
nohup $ENV_GUNICORN wsgi:application -k gevent -w 1 -b 0.0.0.0:$DASH_PORT \
    --timeout 180000 --config gunicorn_config.py --pid gunicorn.pid --access-logfile None&

