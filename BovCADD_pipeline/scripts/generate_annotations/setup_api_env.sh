#!/bin/bash
# Adds manually downloaded perl api modules to the Perl path.
# This way they will be found without specifying the path in the script.
# The script takes one parameter: The base path were the modules are installed.
# Author: Job van Schipstal
# based on https://www.ensembl.org/info/docs/api/api_installation.html

# Allow execution even when PERL5LIB is undefined
set +u

# Add manually installed libraries to path
PERL5LIB=${PERL5LIB}:"$1"/bioperl-1.6.924
PERL5LIB=${PERL5LIB}:"$1"/ensembl-release-110/modules
PERL5LIB=${PERL5LIB}:"$1"/ensembl-compara-release-110/modules
PERL5LIB=${PERL5LIB}:"$1"/ensembl-io-release-110/modules
export PERL5LIB