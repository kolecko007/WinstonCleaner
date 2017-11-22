#!/bin/bash

if [ "$1" == "--local" ]; then
  pip2 install --user -r requirements.txt
else
  pip2 install -r requirements.txt
fi

if [ ! -f config/settings.yml ]; then
    echo "Generating config.yml"
    cp config/settings.yml.default config/settings.yml
fi
