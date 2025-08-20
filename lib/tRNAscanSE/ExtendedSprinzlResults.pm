# tRNAscanSE/ExtendedSprinzlResults.pm
# This class contains functions to write extended Sprinzl position alignment results
# with comprehensive position mapping including sequence, alignment, and Sprinzl coordinates.
#
# --------------------------------------------------------------
# This module extends the tRNAscan-SE program.
# Based on original SprinzlAlignResults.pm
# Copyright (C) 2017 Patricia Chan and Todd Lowe 
# --------------------------------------------------------------
#

package tRNAscanSE::ExtendedSprinzlResults;

use strict;
use tRNAscanSE::tRNA;
use tRNAscanSE::ArraytRNA;
use tRNAscanSE::LogFile;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(write_extended_position_map write_extended_position_db write_genomic_coordinate_map);

sub write_extended_position_map
{
    my ($global_vars, $file, $tRNA) = @_;
    my $log = $global_vars->{log_file};
    
    $log->status("Writing extended position map $file");
    open(FILE_OUT, ">$file") or die "Error: Fail to open $file\n";
    
    # Write header
    print FILE_OUT "tRNA_ID\tSeq_Pos\tNucleotide\tAlign_Pos\tSprinzl_Pos\tIsotype\tAnticodon\tScore\n";
    
    # Get tRNA data
    my $tRNA_id = $tRNA->tRNAscan_id();
    my $isotype = $tRNA->isotype();
    my $anticodon = $tRNA->anticodon();
    my $score = sprintf("%.2f", $tRNA->score());
    my $seq = $tRNA->seq();
    my $sprinzl_align = $tRNA->sprinzl_align();
    my @ar_sprinzl_map = $tRNA->ar_pos_sprinzl_map();
    
    # Map each position
    my $seq_pos = 1;
    for (my $align_pos = 0; $align_pos < length($sprinzl_align); $align_pos++)
    {
        my $nucleotide = substr($sprinzl_align, $align_pos, 1);
        my $sprinzl_pos = $ar_sprinzl_map[$align_pos];
        
        # Only output non-gap positions from original sequence
        if ($nucleotide ne '-')
        {
            print FILE_OUT "$tRNA_id\t$seq_pos\t$nucleotide\t$align_pos\t$sprinzl_pos\t$isotype\t$anticodon\t$score\n";
            $seq_pos++;
        }
    }
    
    close(FILE_OUT);
}

sub write_extended_position_db
{
    my ($global_vars, $file) = @_;
    my $log = $global_vars->{log_file};
    my $tRNAs = $global_vars->{tRNAs};
    
    $tRNAs->sort_array("tRNAscan_id");
    
    $log->status("Writing extended position database $file");
    open(FILE_OUT, ">$file") or die "Error: Fail to open $file\n";
    
    # Write header
    print FILE_OUT "tRNA_ID\tSeq_Pos\tNucleotide\tAlign_Pos\tSprinzl_Pos\tIsotype\tAnticodon\tScore\tGenomic_Start\tGenomic_End\tStrand\tSequence_Name\n";
    
    for (my $i = 0; $i < $tRNAs->get_count(); $i++)
    {
        my $tRNA = $tRNAs->get($i);
        
        # Get tRNA metadata
        my $tRNA_id = $tRNA->tRNAscan_id();
        my $isotype = $tRNA->isotype();
        my $anticodon = $tRNA->anticodon();
        my $score = sprintf("%.2f", $tRNA->score());
        my $genomic_start = $tRNA->start();
        my $genomic_end = $tRNA->end();
        my $strand = $tRNA->strand();
        my $seqname = $tRNA->seqname();
        
        my $seq = $tRNA->seq();
        my $sprinzl_align = $tRNA->sprinzl_align();
        my @ar_sprinzl_map = $tRNA->ar_pos_sprinzl_map();
        
        # Map each position
        my $seq_pos = 1;
        for (my $align_pos = 0; $align_pos < length($sprinzl_align); $align_pos++)
        {
            my $nucleotide = substr($sprinzl_align, $align_pos, 1);
            my $sprinzl_pos = $ar_sprinzl_map[$align_pos];
            
            # Only output non-gap positions from original sequence
            if ($nucleotide ne '-')
            {
                print FILE_OUT "$tRNA_id\t$seq_pos\t$nucleotide\t$align_pos\t$sprinzl_pos\t$isotype\t$anticodon\t$score\t$genomic_start\t$genomic_end\t$strand\t$seqname\n";
                $seq_pos++;
            }
        }
    }
    
    close(FILE_OUT);
}

# Additional function for genomic coordinate mapping
sub write_genomic_coordinate_map
{
    my ($global_vars, $file) = @_;
    my $log = $global_vars->{log_file};
    my $tRNAs = $global_vars->{tRNAs};
    
    $tRNAs->sort_array("tRNAscan_id");
    
    $log->status("Writing genomic coordinate map $file");
    open(FILE_OUT, ">$file") or die "Error: Fail to open $file\n";
    
    # Write header
    print FILE_OUT "tRNA_ID\tGenomic_Coord\tSeq_Pos\tNucleotide\tAlign_Pos\tSprinzl_Pos\tIsotype\tAnticodon\tStrand\tSequence_Name\n";
    
    for (my $i = 0; $i < $tRNAs->get_count(); $i++)
    {
        my $tRNA = $tRNAs->get($i);
        
        # Get tRNA metadata
        my $tRNA_id = $tRNA->tRNAscan_id();
        my $isotype = $tRNA->isotype();
        my $anticodon = $tRNA->anticodon();
        my $strand = $tRNA->strand();
        my $seqname = $tRNA->seqname();
        my $genomic_start = $tRNA->start();
        my $genomic_end = $tRNA->end();
        
        my $seq = $tRNA->seq();
        my $sprinzl_align = $tRNA->sprinzl_align();
        my @ar_sprinzl_map = $tRNA->ar_pos_sprinzl_map();
        
        # Map each position with genomic coordinates
        my $seq_pos = 1;
        for (my $align_pos = 0; $align_pos < length($sprinzl_align); $align_pos++)
        {
            my $nucleotide = substr($sprinzl_align, $align_pos, 1);
            my $sprinzl_pos = $ar_sprinzl_map[$align_pos];
            
            # Only output non-gap positions from original sequence
            if ($nucleotide ne '-')
            {
                # Calculate genomic coordinate
                my $genomic_coord;
                if ($strand eq '+')
                {
                    $genomic_coord = $genomic_start + $seq_pos - 1;
                }
                else
                {
                    $genomic_coord = $genomic_end - $seq_pos + 1;
                }
                
                print FILE_OUT "$tRNA_id\t$genomic_coord\t$seq_pos\t$nucleotide\t$align_pos\t$sprinzl_pos\t$isotype\t$anticodon\t$strand\t$seqname\n";
                $seq_pos++;
            }
        }
    }
    
    close(FILE_OUT);
}

1;