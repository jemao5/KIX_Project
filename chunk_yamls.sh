#!/bin/bash

YAML_DIR=/scratch/jem9759/ZhangWork/KIX_Project/yaml_outputs/out_yaml_full_library
CHUNK_BASE=/scratch/jem9759/ZhangWork/KIX_Project/yaml_chunks
CHUNK_SIZE=500   # tune based on your time limit

mkdir -p "$CHUNK_BASE"
i=0
chunk=0
find "$YAML_DIR" -maxdepth 1 -name "*.yaml" | sort | while read f; do
    if (( i % CHUNK_SIZE == 0 )); then
        chunk_dir=$(printf "%s/chunk_%04d" "$CHUNK_BASE" $((i / CHUNK_SIZE)))
        mkdir -p "$chunk_dir"
    fi
    ln -s "$f" "$chunk_dir/"
    ((i++))
done
echo "Split into $(( (31392 + CHUNK_SIZE - 1) / CHUNK_SIZE )) chunks"
