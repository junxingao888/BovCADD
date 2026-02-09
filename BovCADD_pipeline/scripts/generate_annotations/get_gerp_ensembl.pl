#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# :Edited: 20-10-2023
# Changes:
# - re-use code for both gerp scores and constrained elements, only 1 script
# - Perform polling in batches, helps for the longest chromosomes.
# - Make alignment configurable (can use others instead of mammals)

use strict;
use warnings;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Exception qw(throw);

#
# Simple example to show how to get conservation scores for a slice.
#
# Imports Gepopt package which is doing the command line handling
use Getopt::Long;

# Save arguments in the variable names utilized by this script
GetOptions("chr=s" => \my $chr
    , "start=i"    => \my $region_start
    , "end=i"      => \my $region_end
    , "sp=s"       => \my $species
    , "set=s"      => \my $set
    , "elem"       => \my $is_elem
);

my $batch_size = $is_elem ? 50000000 : 1000000; # keeps memory usage below 1GB

my $reg = "Bio::EnsEMBL::Registry";

$reg->load_registry_from_db(
    -host => "ensembldb.ensembl.org",
    -user => "anonymous",
    #    -VERBOSE => 1
);

my $mlss_adaptor = $reg->get_adaptor("Multi", "compara", "MethodLinkSpeciesSet");
my $mlss = $mlss_adaptor->fetch_by_method_link_type_species_set_name($is_elem ? "GERP_CONSTRAINED_ELEMENT" : "GERP_CONSERVATION_SCORE", $set);
throw("Unable to find method_link_species_set") if (!defined($mlss));

my $slice_adaptor = $reg->get_adaptor($species, 'core', 'Slice');
throw("Registry configuration file has no data for connecting to <$species>") if (!$slice_adaptor);

my $cs_adaptor = $reg->get_adaptor("Multi", 'compara', $is_elem ? 'ConstrainedElement' : 'ConservationScore');

# Calculate the total length of the slice
my $slice_length = $region_end - $region_start + 1;

# Determine the number of batches needed
my $num_batches = int(($slice_length + $batch_size - 1) / $batch_size);

print STDERR "Fetching GERP scores from $region_start to $region_end, $chr in $num_batches batches.\n";

# Fetching scores for 1Mb at a time, write to file in 0-based bed-graph style
for my $batch_num (1 .. $num_batches) {
    my $batch_start = ($batch_num - 1) * $batch_size;
    my $batch_end = $batch_num * $batch_size;
    $batch_end = $slice_length if $batch_end > $slice_length;

    my $batch_slice = $slice_adaptor->fetch_by_region('toplevel', $chr, $region_start + $batch_start, $region_start + $batch_end - 1);
    throw("No Slice can be created with coordinates $chr:$region_start-$region_end") if (!$batch_slice);

    my $display_size = $batch_slice->end - $batch_slice->start + 1;
    my $scores = $cs_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $batch_slice, $display_size);

    print STDERR "Batch $batch_num: #nr of scores " . @$scores . "\n";

    # Print out the region/chr, position and found scores
    if (defined $is_elem && $is_elem) {
        foreach my $score (@$scores) {
            printf("%s\t%d\t%d\t%.4f\t%.4e\n",
                $chr,
                $score->start + $batch_start - 2,
                $score->end + $batch_start - 1,
                $score->score,
                $score->p_value);
        }
    }
    else {
        foreach my $score (@$scores) {
            if (defined $score->diff_score) {
                printf("%s\t%d\t%d\t%.4f\t%.4f\t%.4f\n",
                    $chr,
                    $score->position + $batch_start - 2,
                    $score->position + $batch_start - 1,
                    $score->observed_score,
                    $score->expected_score,
                    $score->diff_score);
            }
        }
    }
}