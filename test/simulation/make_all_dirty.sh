#!/bin/bash

# requires:
# (1) base simulation folder (ex. contamination_simulation/)
# (2) contamination folder (ex. brucei_vs_giardia_10_pct/)
# (3) first dataset id
# (4) second dataset id
# (5) percentage

set -e

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ -z "$5" ]
  then
    echo "Wrong arguments supplied"
    exit 1
fi

BASE_FOLDER=$(readlink -f $1)
CONT_FOLDER=$(readlink -f $2)
PERCENTAGE=$5

test/simulation/make_dirty.sh $BASE_FOLDER $CONT_FOLDER $3 $4 $PERCENTAGE
test/simulation/make_dirty.sh $BASE_FOLDER $CONT_FOLDER $4 $3 $PERCENTAGE

echo "make_all_dirty.sh finished succesfully"
