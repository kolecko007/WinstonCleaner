#!/bin/bash
set -e

rm -rf test/integration/output/*
bin/prepare_data.py --config_path test/integration/settings.yml

bin/find_contaminations.py --config_path test/integration/settings.yml
