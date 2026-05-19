#!/usr/bin/env bash

# Copy this file to private_env.sh, fill in your own keys, and source it before
# running the workflow. private_env.sh is git-ignored. Do not commit real API
# keys to GitHub.
export OPENTOPOGRAPHY_API_KEY="<your-opentopography-api-key>"

# cdsapi accepts these environment variables instead of ~/.cdsapirc.
export CDSAPI_URL="https://cds.climate.copernicus.eu/api"
export CDSAPI_KEY="<your-cds-api-key>"
