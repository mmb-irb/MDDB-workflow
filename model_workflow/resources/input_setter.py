#!/usr/bin/env python
# coding: utf-8
# Temporary until mwf reads YML

import sys
import json
from yaml import load, CLoader


try:
    with open(sys.argv[1])  as config_file:
        inputs = load(config_file, Loader=CLoader)
except Exception as e:
    sys.exit('Usage: input_setter.py inputs.yml inputs.json')


# Export it to json
inputs_filename = sys.argv[2]
with open(inputs_filename, 'w') as file:
    json.dump(inputs, file, indent=4)

