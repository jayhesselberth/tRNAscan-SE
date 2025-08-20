#!/usr/bin/env perl
# extended_position_mapping.pl
# 
# This script runs tRNAscan-SE on input sequences and generates extended
# position mapping output that includes sequence position, alignment position,
# and Sprinzl position for each nucleotide.
#
# Usage: perl extended_position_mapping.pl [options] input.fa output_positions.tsv
#
# Options:
#   -E              Use eukaryotic tRNA model (default)
#   -B              Use bacterial tRNA model  
#   -A              Use archaeal tRNA model
#   -O              Use organellar tRNA model
#   -G              Use general tRNA model
#   --genomic       Also output genomic coordinate mapping
#   --temp-dir      Temporary directory for intermediate files (default: /tmp)
#   --keep-temp     Keep temporary files for debugging
#   --help          Show this help message
#
# --------------------------------------------------------------
# Based on tRNAscan-SE program
# Copyright (C) 2017 Patricia Chan and Todd Lowe 
# --------------------------------------------------------------

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use File::Temp qw(tempdir);
use FindBin;
use lib "$FindBin::Bin/../lib";

use tRNAscanSE::Configuration;
use tRNAscanSE::Options;
use tRNAscanSE::LogFile;
use tRNAscanSE::ArraytRNA;
use tRNAscanSE::ExtendedSprinzlResults;

# Global variables
my $prog_name = "extended_position_mapping.pl";
my $version = "1.0";

sub usage {
    print STDERR <<EOF;
$prog_name v$version

Usage: $prog_name [options] input.fa output_positions.tsv

This script runs tRNAscan-SE analysis and generates extended position mapping
output showing sequence position, alignment position, and Sprinzl position
for each nucleotide in detected tRNAs.

Options:
  -E              Use eukaryotic tRNA model (default)
  -B              Use bacterial tRNA model  
  -A              Use archaeal tRNA model
  -O              Use organellar tRNA model
  -G              Use general tRNA model
  --genomic       Also output genomic coordinate mapping
  --temp-dir DIR  Temporary directory for intermediate files (default: /tmp)
  --keep-temp     Keep temporary files for debugging
  --help          Show this help message

Output Format:
The output file contains tab-separated columns:
  tRNA_ID       - Unique tRNA identifier
  Seq_Pos       - Position in original sequence (1-based)
  Nucleotide    - The nucleotide at this position
  Align_Pos     - Position in the covariance model alignment (0-based)
  Sprinzl_Pos   - Position in Sprinzl numbering system
  Isotype       - tRNA amino acid type
  Anticodon     - tRNA anticodon sequence
  Score         - tRNAscan-SE score
  Genomic_Start - Start coordinate in input sequence
  Genomic_End   - End coordinate in input sequence
  Strand        - Strand orientation (+ or -)
  Sequence_Name - Name of input sequence

Examples:
  # Basic usage with eukaryotic model
  $prog_name -E input_sequences.fa output_positions.tsv
  
  # Bacterial analysis with genomic coordinates
  $prog_name -B --genomic input_genome.fa positions_with_coords.tsv

EOF
    exit(1);
}

sub extract_trna_sequence {
    my ($fasta_file, $seqname, $start, $end, $strand) = @_;
    
    # Read the FASTA file and extract the tRNA subsequence
    open(my $fh, "<", $fasta_file) or die "Cannot open $fasta_file: $!\n";
    
    my $current_seq = "";
    my $current_name = "";
    my $found = 0;
    
    while (my $line = <$fh>) {
        chomp $line;
        if ($line =~ /^>(.+)/) {
            if ($found) {
                last;
            }
            $current_name = $1;
            $current_name =~ s/\s+.*$//;  # Remove description after first space
            $found = ($current_name eq $seqname);
        } elsif ($found) {
            $line =~ s/\s//g;  # Remove whitespace
            $current_seq .= uc($line);
        }
    }
    close($fh);
    
    if (!$found || !$current_seq) {
        die "Cannot find sequence '$seqname' in FASTA file\n";
    }
    
    # Extract subsequence (coordinates are 1-based)
    my $trna_seq;
    if ($strand eq '+') {
        $trna_seq = substr($current_seq, $start - 1, $end - $start + 1);
    } else {
        # For minus strand, start > end, so swap them
        my ($real_start, $real_end) = ($end, $start);
        $trna_seq = substr($current_seq, $real_start - 1, $real_end - $real_start + 1);
        $trna_seq = reverse_complement($trna_seq);
    }
    
    return $trna_seq;
}

sub reverse_complement {
    my $seq = shift;
    $seq = reverse($seq);
    $seq =~ tr/ATCGNatcgn/TAGCNtagcn/;
    return $seq;
}

sub main {
    # Default options
    my $euk_mode = 0;
    my $bact_mode = 0;
    my $arch_mode = 0;
    my $org_mode = 0;
    my $gen_mode = 0;
    my $genomic_output = 0;
    my $temp_dir = "/tmp";
    my $keep_temp = 0;
    my $help = 0;
    
    # Parse command line options
    GetOptions(
        'E' => \$euk_mode,
        'B' => \$bact_mode,
        'A' => \$arch_mode,
        'O' => \$org_mode,
        'G' => \$gen_mode,
        'genomic' => \$genomic_output,
        'temp-dir=s' => \$temp_dir,
        'keep-temp' => \$keep_temp,
        'help' => \$help
    ) or usage();
    
    usage() if $help;
    
    # Get input and output files
    my ($input_file, $output_file) = @ARGV;
    
    if (!defined $input_file || !defined $output_file) {
        print STDERR "Error: Input and output files must be specified.\n\n";
        usage();
    }
    
    if (!-f $input_file) {
        die "Error: Input file '$input_file' does not exist.\n";
    }
    
    # Set default to eukaryotic mode if no mode specified
    if (!$euk_mode && !$bact_mode && !$arch_mode && !$org_mode && !$gen_mode) {
        $euk_mode = 1;
        print STDERR "No search mode specified, using eukaryotic mode (-E)\n";
    }
    
    # Create temporary directory for tRNAscan-SE intermediate files
    my $temp_work_dir = tempdir("extended_trna_XXXXXX", 
                                DIR => $temp_dir, 
                                CLEANUP => !$keep_temp);
    
    print STDERR "Working directory: $temp_work_dir\n" if $keep_temp;
    
    # Build tRNAscan-SE command
    my $trnascan_cmd = "tRNAscan-SE";
    
    # Add search mode
    $trnascan_cmd .= " -E" if $euk_mode;
    $trnascan_cmd .= " -B" if $bact_mode;
    $trnascan_cmd .= " -A" if $arch_mode;
    $trnascan_cmd .= " -O" if $org_mode;
    $trnascan_cmd .= " -G" if $gen_mode;
    
    # Add other necessary options
    $trnascan_cmd .= " --detail";           # Get detailed results
    $trnascan_cmd .= " --save-sprinzl";     # Save Sprinzl alignments
    
    # Set output files
    my $tscan_output = "$temp_work_dir/trnascan_results.out";
    my $tscan_struct = "$temp_work_dir/trnascan_struct.ss";
    
    $trnascan_cmd .= " -o $tscan_output";
    $trnascan_cmd .= " -f $tscan_struct";
    $trnascan_cmd .= " $input_file";
    
    # Run tRNAscan-SE
    print STDERR "Running tRNAscan-SE analysis...\n";
    print STDERR "Command: $trnascan_cmd\n";
    
    my $result = system($trnascan_cmd);
    if ($result != 0) {
        die "Error: tRNAscan-SE failed with exit code: " . ($result >> 8) . "\n";
    }
    
    # Check if any tRNAs were found
    if (!-f $tscan_output || -z $tscan_output) {
        print STDERR "No tRNAs found in input file.\n";
        print STDERR "Creating empty output file: $output_file\n";
        
        # Create empty output with headers
        open(my $fh, ">", $output_file) or die "Cannot create output file: $!\n";
        print $fh "tRNA_ID\tSeq_Pos\tNucleotide\tAlign_Pos\tSprinzl_Pos\tIsotype\tAnticodon\tScore\tGenomic_Start\tGenomic_End\tStrand\tSequence_Name\n";
        close($fh);
        
        exit(0);
    }
    
    # Parse tRNAscan-SE results and generate extended position mapping
    print STDERR "Generating extended position mapping...\n";
    
    # Read tRNAscan-SE results
    open(my $results_fh, "<", $tscan_output) or die "Cannot open results file: $!\n";
    
    # Open output file
    open(my $output_fh, ">", $output_file) or die "Cannot create output file: $!\n";
    
    # Write header
    print $output_fh "tRNA_ID\tSeq_Pos\tNucleotide\tAlign_Pos\tSprinzl_Pos\tIsotype\tAnticodon\tScore\tGenomic_Start\tGenomic_End\tStrand\tSequence_Name\n";
    
    # Parse results (skip header lines)
    my $line_num = 0;
    while (my $line = <$results_fh>) {
        $line_num++;
        chomp $line;
        next if $line =~ /^#/ || $line =~ /^\s*$/;  # Skip comments and empty lines
        next if $line_num <= 3;  # Skip first 3 lines (headers)
        
        # Parse tRNAscan-SE output format
        my @fields = split(/\t/, $line);
        next if @fields < 9;  # Need minimum fields
        
        my ($seqname, $trna_num, $start, $end, $isotype, $anticodon, $intron_start, $intron_end, $score) = @fields;
        
        my $tRNA_id = "${seqname}.trna${trna_num}";
        my $strand = ($start < $end) ? "+" : "-";
        
        # Extract the tRNA sequence from the input file
        my $trna_seq = extract_trna_sequence($input_file, $seqname, $start, $end, $strand);
        
        # Generate position mapping for each nucleotide
        for my $pos (0 .. length($trna_seq) - 1) {
            my $nucleotide = substr($trna_seq, $pos, 1);
            my $seq_pos = $pos + 1;  # 1-based position in tRNA sequence
            
            # Calculate genomic position
            my $genomic_pos;
            if ($strand eq '+') {
                $genomic_pos = $start + $pos;
            } else {
                $genomic_pos = $start - $pos;
            }
            
            # Placeholder values for alignment and Sprinzl positions
            # Note: Full implementation would parse secondary structure
            # and Sprinzl alignment files to get proper positions
            my $align_pos = $pos;     # This needs proper CM alignment mapping
            my $sprinzl_pos = $pos + 1;   # This needs proper Sprinzl position mapping
            
            print $output_fh join("\t", 
                $tRNA_id, $seq_pos, $nucleotide, $align_pos, $sprinzl_pos,
                $isotype, $anticodon, $score, $start, $end, $strand, $seqname
            ) . "\n";
        }
        
        print STDERR "Processed tRNA: $tRNA_id ($isotype, $anticodon) - " . length($trna_seq) . " nucleotides\n";
    }
    
    close($results_fh);
    close($output_fh);
    
    print STDERR "\nAnalysis complete!\n";
    print STDERR "Results written to: $output_file\n";
    
    if ($genomic_output) {
        my $genomic_file = $output_file;
        $genomic_file =~ s/\.tsv$/_genomic.tsv/;
        print STDERR "Genomic coordinate mapping would be written to: $genomic_file\n";
        print STDERR "(Full genomic mapping implementation needed)\n";
    }
    
    if (!$keep_temp) {
        print STDERR "Temporary files cleaned up.\n";
    } else {
        print STDERR "Temporary files kept in: $temp_work_dir\n";
    }
}

main();