#!/bin/bash
## Basic script to install the ENSEMBL perl api-client
# This API is used in the workflow to download GERP scores
# While the GERP scores themselves are also found on the ftp-site
# the additonal scores, such as the constained element scores are not.

# Author: Job van Schipstal
# based on https://www.ensembl.org/info/docs/api/api_installation.html

wget https://github.com/bioperl/bioperl-live/archive/release-1-6-924.zip
unzip -qo release-1-6-924.zip
rm release-1-6-924.zip
mv bioperl-live-release-1-6-924 bioperl-1.6.924
wget https://github.com/Ensembl/ensembl-compara/archive/refs/heads/release/110.zip
unzip -qo 110.zip
rm 110.zip
wget https://github.com/Ensembl/ensembl-io/archive/refs/heads/release/110.zip
unzip -qo 110.zip
rm 110.zip
wget https://github.com/Ensembl/ensembl/archive/refs/heads/release/110.zip
unzip -qo 110.zip
rm 110.zip