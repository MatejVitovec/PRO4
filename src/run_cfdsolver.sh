#!/bin/bash

# Soubor s argumenty
ARGUMENTS_FILE="arguments.txt"

# Kontrola, zda soubor s argumenty existuje
if [[ ! -f "$ARGUMENTS_FILE" ]]; then
  echo "Soubor $ARGUMENTS_FILE neexistuje."
  exit 1
fi

# Čtení souboru řádek po řádku
while IFS= read -r line
do
  # Spuštění programu s argumenty z aktuálního řádku
  ./main $line

  # Zkontrolovat výsledek spuštění programu (volitelné)
  if [[ $? -ne 0 ]]; then
    echo "Chyba při spuštění programu s argumenty: $line"
  fi
done < "$ARGUMENTS_FILE"