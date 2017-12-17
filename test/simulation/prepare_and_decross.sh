#!/bin/bash

# requires:
# (1) base simulation folder (ex. contamination_simulation/)
# (2) first dataset id
# (3) second dataset id
# (4) percentage

set -e

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ]
  then
    echo "Wrong arguments supplied"
    exit 1
fi

BASE_FOLDER=$(readlink -f $1)
PERCENTAGE=$4

ORGS=($2 $3)
readarray -t SORTED_ORGS < <(for e in "${ORGS[@]}"; do echo "$e"; done | sort)
CONT_FOLDER_NAME=${SORTED_ORGS[0]}_vs_${SORTED_ORGS[1]}_${PERCENTAGE}_pct
CONT_FOLDER_NAME=$(echo $CONT_FOLDER_NAME | sed 's/\./_/')
CONT_FOLDER=$BASE_FOLDER/$CONT_FOLDER_NAME
mkdir $CONT_FOLDER

test/simulation/make_all_dirty.sh $BASE_FOLDER $CONT_FOLDER $2 $3 $4

DEFAULT_SETTINGS_PATH=test/simulation/settings.yml.default
SETTINGS_PATH=test/simulation/conf_$(date | md5sum | awk '{print $1}').yml

INPUT_FOLDER=$CONT_FOLDER/for_decross
OUTPUT_FOLDER=$INPUT_FOLDER/output
mkdir -p $OUTPUT_FOLDER

sed "s#%data_folder%#${INPUT_FOLDER}#" $DEFAULT_SETTINGS_PATH | \
  sed "s#%output_folder%#${OUTPUT_FOLDER}#" > $SETTINGS_PATH

bin/prepare_data.py --config_path $SETTINGS_PATH

sed -i "s/CLOSE/REGULAR/" $OUTPUT_FOLDER/types.csv
bin/find_contaminations.py --config_path $SETTINGS_PATH

rm $SETTINGS_PATH
echo "run_decross.sh finished succesfully"
