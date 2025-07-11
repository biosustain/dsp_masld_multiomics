#!/usr/bin/env perl
# Juliana A
# count Start and Stop Codons in nucleotide sequence
#
use strict;
use warnings;
use Bio::SeqIO;
use File::Basename;

# Check for correct number of arguments
if (@ARGV != 1) {
    die "Usage: $0 <input_fasta_file>\n";
}

# Input file passed as argument
my ($f1) = @ARGV;

# Generate the output file name: input filename + "_stats.tsv"
my $output_base = basename($f1, ".fa", ".fasta");  # Remove .fa or .fasta if present
my $tsv_out = "${output_base}_stats.tsv";
my $internal_stop_out = "${output_base}_INTERNAL_STOPS.fa";

# Create SeqIO objects for input and output
my $in = Bio::SeqIO->new(-file => "<$f1");
my $out = Bio::SeqIO->new(-file => ">$internal_stop_out", -format => "fasta");

# Open the TSV output file for writing
open(my $tsv_fh, '>', $tsv_out) or die "Could not open file '$tsv_out' $!";

# Write the header to the TSV file
print $tsv_fh join("\t", qw/ID FirstAA Length StartCheck StopCheck InternalStops/) . "\n";

# Process each sequence
while (my $s = $in->next_seq) {
    my $id = $s->primary_id;
    my $alpha = $s->alphabet;
    my $aa = $s->seq;

    # Get the full sequence length before removing start/stop codons
    my $length = length($aa);

    # Check for start and stop codons
    my $start = ($aa =~ /^ATG/) ? '1' : '0';
    my $end = ($aa =~ /TAA$|TAG$|TGA$/) ? '1' : '0';

    # Remove the first and last characters from the sequence (start and stop codons)
    $aa =~ s/^.{1}//;  # Remove the first character (start codon)
    $aa =~ s/.{1}$//;  # Remove the last character (stop codon)

    # Count internal stop codons
    my @stops = ($aa =~ /(\.|\*)/g);
    my $stop = scalar @stops;

    # Write results to the TSV file, using the full length before codon removal
    print $tsv_fh join("\t", $id, $alpha || '.', $length, $start, $end, $stop || '0') . "\n";

    # Write sequences with internal stop codons to the output FASTA file
    $out->write_seq($s) if $stop;
}

# Close the TSV file
close($tsv_fh);

# Print completion message with output filenames
print "Results written to $tsv_out and $internal_stop_out\n";

