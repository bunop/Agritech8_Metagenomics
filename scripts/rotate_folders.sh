#!/bin/bash
set -euxo pipefail

# check for received argument
if [ -z "${1:-}" ]; then
    echo "Usage: $0 <base-folder-name>"
    exit 1
fi

# Ottieni il nome base della cartella da argomento
BASE_NAME="$1"

# Trova tutte le cartelle che iniziano con il nome base passato come argomento e ordinale in ordine decrescente
folders=$(ls -d ${BASE_NAME}* 2>/dev/null | sort -r)

# Verifica se ci sono cartelle da rinominare
if [ -z "$folders" ]; then
    echo "No folders found starting with $BASE_NAME"
    exit 1
fi

# Rinomina le cartelle in ordine decrescente per evitare sovrascrizioni
for folder in $folders; do
    # Estrai il numero alla fine del nome della cartella (se esiste)
    if [[ $folder =~ ${BASE_NAME}\.([0-9]+)$ ]]; then
        num=${BASH_REMATCH[1]}
        new_num=$((num + 1))
        mv "$folder" "${BASE_NAME}.$new_num"
    elif [[ $folder == "$BASE_NAME" ]]; then
        mv "$folder" "${BASE_NAME}.1"
    fi
done
