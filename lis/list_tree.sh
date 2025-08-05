#!/bin/bash
# File: list_tree.sh
# Run: ./list_tree.sh > list_tree.txt

# Set base directory relative to where script is run
BASE_DIR="$(pwd)"

echo "Scanning LIS directory structure from: $BASE_DIR"

# Find and list relevant directories/files
find "$BASE_DIR" -type d | while read dir; do
    echo -e "\nDIR: $dir"
    ls -1 "$dir" | grep -E 'GPS|GRACE'  # Show GPS/GRACE-related files
done
