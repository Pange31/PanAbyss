#!/bin/bash

set -e

BATCH_SIZE=${1:-2000000}

echo "Using batch size: $BATCH_SIZE"

if ! [[ "$BATCH_SIZE" =~ ^[0-9]+$ ]]; then
    echo "ERROR: batch size must be an integer"
    exit 1
fi

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

GFA_DIR="$SCRIPT_DIR/data/gfa"
IMPORT_DIR="$SCRIPT_DIR/data/import"
CHROMOSOME_CSV="$SCRIPT_DIR/data/gfa/chromosomes_file.csv"

echo "Generating CSV import files from GFA..."

# ------------------------------------------------------------
# Detect GFA files
# ------------------------------------------------------------
shopt -s nullglob
gfa_files=("$GFA_DIR"/*.gfa)

gfa_count=${#gfa_files[@]}

if [ "$gfa_count" -eq 0 ]; then
    echo "ERROR: No .gfa files found in $GFA_DIR"
    exit 1
fi

# ------------------------------------------------------------
# Single GFA file
# ------------------------------------------------------------
if [ "$gfa_count" -eq 1 ]; then

    file_path="${gfa_files[0]}"

    echo "Processing single GFA file:"
    echo "  $file_path"

    python <<EOF
import sys
from pathlib import Path

project_dir = Path(r'''$SCRIPT_DIR''')
sys.path.insert(0, str(project_dir))

from neo4j_DB_construction import load_gfa_data_to_csv

load_gfa_data_to_csv(
    r'''$file_path''',
    import_dir=r'''$IMPORT_DIR''',
    chromosome_file="",
    chromosome_prefix=False,
    batch_size=$BATCH_SIZE,
    start_chromosome=None,
    haplotype=True
)
EOF

    echo "CSV generation completed."
    exit 0
fi

# ------------------------------------------------------------
# Multiple GFA files
# ------------------------------------------------------------
echo "Multiple GFA files detected."

if [ ! -f "$CHROMOSOME_CSV" ]; then
    echo "ERROR: $CHROMOSOME_CSV is required when multiple GFA files are present."
    echo ""
    echo "Expected structure:"
    echo "filename,chromosome"
    echo "chr1.gfa,chr1"
    echo "chr2.gfa,chr2"
    echo ""
    echo "Supported separators: comma, semicolon, tab"
    exit 1
fi

echo "Using chromosome mapping file:"
echo "  $CHROMOSOME_CSV"

# ------------------------------------------------------------
# Detect separator
# ------------------------------------------------------------
header=$(head -n 1 "$CHROMOSOME_CSV")

if [[ "$header" == *$'\t'* ]]; then
    separator=$'\t'
elif [[ "$header" == *";"* ]]; then
    separator=";"
elif [[ "$header" == *","* ]]; then
    separator=","
else
    echo "ERROR: Unable to detect CSV separator."
    echo "Supported separators: comma, semicolon, tab"
    exit 1
fi

echo "Detected separator: '$separator'"

# ------------------------------------------------------------
# Detect header
# ------------------------------------------------------------
first_line=$(head -n 1 "$CHROMOSOME_CSV" | tr -d '\r')

if [[ "$first_line" =~ filename ]] && [[ "$first_line" =~ chromosome ]]; then
    echo "Header detected in chromosome mapping file."
    input_lines=$(tail -n +2 "$CHROMOSOME_CSV")
else
    echo "No header detected in chromosome mapping file."
    input_lines=$(cat "$CHROMOSOME_CSV")
fi

# ------------------------------------------------------------
# Validation pass
# ------------------------------------------------------------
declare -a valid_files
declare -a valid_chromosomes

validation_error=0

while IFS="$separator" read -r filename chromosome
do

    # Remove CR characters
    filename=$(echo "$filename" | tr -d '\r')
    chromosome=$(echo "$chromosome" | tr -d '\r')

    # Skip empty lines
    if [ -z "$filename" ] && [ -z "$chromosome" ]; then
        continue
    fi

    if [ -z "$filename" ]; then
        echo "ERROR: Missing filename in CSV."
        validation_error=1
        continue
    fi

    if [ -z "$chromosome" ]; then
        echo "ERROR: Missing chromosome for file '$filename'."
        validation_error=1
        continue
    fi

    file_path="$GFA_DIR/$filename"

    if [ ! -f "$file_path" ]; then
        echo "ERROR: File not found:"
        echo "  $file_path"
        validation_error=1
        continue
    fi

    valid_files+=("$file_path")
    valid_chromosomes+=("$chromosome")

done <<< "$input_lines"

if [ "$validation_error" -ne 0 ]; then
    echo ""
    echo "CSV validation failed."
    exit 1
fi

# ------------------------------------------------------------
# Processing pass
# ------------------------------------------------------------
for i in "${!valid_files[@]}"; do

    file_path="${valid_files[$i]}"
    chromosome="${valid_chromosomes[$i]}"

    echo ""
    echo "Processing:"
    echo "  File       : $file_path"
    echo "  Chromosome : $chromosome"

    python <<EOF
import sys
from pathlib import Path

project_dir = Path(r'''$SCRIPT_DIR''')
sys.path.insert(0, str(project_dir))

from neo4j_DB_construction import load_gfa_data_to_csv

load_gfa_data_to_csv(
    r'''$file_path''',
    import_dir=r'''$IMPORT_DIR''',
    chromosome_file=r'''$chromosome''',
    chromosome_prefix=False,
    batch_size=$BATCH_SIZE,
    start_chromosome=None,
    haplotype=True
)
EOF

done

echo ""
echo "CSV generation completed."
