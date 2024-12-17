#!/bin/bash

set -e
set -u

script=${1}

zip -r ${script}_freeze.zip _freeze

