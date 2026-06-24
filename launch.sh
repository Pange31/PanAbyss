#!/bin/bash
# --- Parameters ---
# --- Parameters ---
DASH_PORT=8050
GENERATE_CSV_IMPORT=0
CREATE_DATABASE=0
DATABASE_NAME=""
BATCH_SIZE=""

while [ $# -gt 0 ]; do
    case "$1" in
        --generate_csv_import)
            GENERATE_CSV_IMPORT=1
            shift
            ;;
        --create_database)
            CREATE_DATABASE=1
            DATABASE_NAME="DB_1.0.0_$2"
            shift 2
            ;;
        --batch_size)
            BATCH_SIZE="$2"
            shift 2
            ;;
        --help|-h)
            echo ""
            echo "PanAbyss launcher"
            echo ""
            echo "Usage:"
            echo "  ./launch.sh [PORT]"
            echo "  ./launch.sh --generate_csv_import"
            echo "  ./launch.sh --create_database <database_name>"
            echo "  ./launch.sh --batch_size <size> --generate_csv_import"
            echo "  ./launch.sh --batch_size <size> --create_database <database_name>"
            echo "  ./launch.sh --help"
            echo ""
            echo "Options:"
            echo "  PORT"
            echo "      Launch Dash application on specified port"
            echo "      Default: 8050"
            echo ""
            echo "  --create_database <database_name>"
            echo "      Create a Neo4j database named <database_name>"
            echo ""
            echo "  --generate_csv_import"
            echo "      Convert all .gfa files from:"
            echo "          ./data/gfa/"
            echo "      into CSV import files in:"
            echo "          ./data/import/"
            echo ""
            echo "      If multiple .gfa files are present in ./data/gfa/,"
            echo "      a file named:"
            echo "          ./data/gfa/chromosomes_file.csv"
            echo "      must also be present."
            echo ""
            echo "      Supported separators:"
            echo "          comma (,)"
            echo "          semicolon (;)"
            echo "          tab"
            echo ""
            echo "      Expected structure:"
            echo ""
            echo "          filename,chromosome"
            echo "          chr1.gfa,chr1"
            echo "          chr2.gfa,chr2"
            echo ""
            echo "  --batch_size <size>"
            echo "      Set batch size used for CSV generation from GFA files"
            echo "      Default: 2000000"
            echo ""
            echo "  --help, -h"
            echo "      Show this help message"
            echo ""
            exit 0
            ;;
        *)
            DASH_PORT="$1"
            shift
            ;;
    esac
done

echo "GENERATE_CSV_IMPORT=$GENERATE_CSV_IMPORT"
echo "CREATE_DATABASE=$CREATE_DATABASE"
echo "DATABASE_NAME=$DATABASE_NAME"
echo "BATCH_SIZE=$BATCH_SIZE"

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
if [ -f "data/conf/neo4j.conf" ]; then
	:
else
	# neo4j conf file doesn't exists => create default neo4j conf file
	echo "Copy neo4j conf file install/conf/neo4j.conf to data/conf/neo4j.conf..."
	cp "install/conf/neo4j.conf" "data/conf/neo4j.conf"
	if [ $? -ne 0 ]; then
		echo "Error when copying neo4j conf file."
		exit 1
	fi
fi


# --- Special mode: generate CSV import from GFA files ---
if [ "$GENERATE_CSV_IMPORT" -eq 1 ]; then
    if [ -n "$BATCH_SIZE" ]; then
      bash ./scripts/generate_csv_import_file.sh "$BATCH_SIZE"
    else
        bash ./scripts/generate_csv_import_file.sh
    fi
    exit $?
fi

# --- Special mode: create database ---
if [ "${CREATE_DATABASE:-0}" -eq 1 ]; then
    echo "WARNING: This step requires Apptainer to be loaded."

    if [ -n "$BATCH_SIZE" ]; then
        bash ./scripts/generate_csv_import_file.sh "$BATCH_SIZE" || exit $?
    else
        bash ./scripts/generate_csv_import_file.sh || exit $?
    fi

    python3 - << EOF
from neo4j_container_management import create_db
create_db("${DATABASE_NAME}", docker=False)
EOF

    exit $?
fi

# --- Launch application ---
echo "Launching PanAbyss on port $DASH_PORT..."
python index.py --port $DASH_PORT
