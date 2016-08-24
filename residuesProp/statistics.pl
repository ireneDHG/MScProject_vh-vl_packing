#!usr/bin/perl -w
#********************************************************************
# MSc Bioinformatics with Systems Biology 2015/16
# 
# Analysis of the VH/VL packing
#
# Author:       Irene del Hierro García
# Script name:  statistics.pl
# Version:      V1.1
# Date:         19/08/16
#
#********************************************************************
#
# Purpose
# =======
#        Calculate the proportions of the amino acids in every position of a given number of positions.
#
#********************************************************************
#
# Assumptions
# ===========
#       1. The program requires an input file in tabular format
#       2. The first column should be the pdb, and the nexts the positions with the residues below as:
#               PDB	L87	L46	L41	H105
#               12E8_1	PHE	LEU	GLY	GLN
#       3. The results will be printed in screen if an output file is not specified.
#
#       Usage: perl statistics.pl < free_Str_H2.txt
#
#********************************************************************   
# Strategy
# ========
#       1. Open and take out the header of the file
#       2. Store the positions and extract the amino acids from the file
#       3. Add the amino acids to their corresponding positions with a hash
#       4. Calculate the number of each amino acid in each position
#       5. Calculate the proportion of the global presence of the residues
#
#********************************************************************
#
# Revision history
# ================ 
#       V1.0:   Original
#       V1.1:   Comments
#
#********************************************************************

use strict;


# Take out the header of the file, and store the names of the positions.
my $header = <>;
my @elements = ($header=~ /\t?(\w+)\t?/g);

my %resCounts = ();
my @counts = ();

# Add the positions into a hash (avoid the PDB header)
for(my $i=1;$i<scalar @elements;$i++)
{      
        my $res = $elements[$i];
        $resCounts{$res}=();
}

# For the rest of the lines of the input file, take the residues
while (my $line = <>)
{
        chomp $line;
        my @residues = ($line=~ /\t?(\w+)\t?/g);
        
        for(my $i=1;$i<scalar @elements;$i++)
        {
                for(my $j=1;$j<scalar @residues;$j++)
                {
                        if($i == $j)
                        {
                                # If the position corresponds with the residue, add it to the hash 
                                foreach my $key (keys %resCounts)
                                {
                                        if($elements[$i] eq $key)
                                        {
                                                push @{$resCounts{$key}},$residues[$i];
                                        }
                                }
                        }
                }
        }
}

# Print the results 
print "RESULTS\n##############\n\n";
print "Frequency of residues in the positions determined\n\n";

my %total = ();

# Once the counts are obtained, make the proportions
foreach my $key (keys %resCounts)
{
        print "$key = ";
        
        my %counts = ();
        my $length = scalar @{$resCounts{$key}};

        foreach my $res (@{$resCounts{$key}})
        {
                $counts{$res}++;
                $total{$res}++;
        }
        # Calculate the fraction of each amino acid present
        foreach my $c (keys %counts)
        {
                my $fraction = $counts{$c} / $length;
                print "\t$c $fraction";
        }
        print "\n";
}         

print "\nTotal amino acids:\n";

my $countT=0;

foreach my $key (sort keys %total)
{
        $countT += $total{$key};
}

foreach my $key (sort keys %total)
{
        my $norm = $total{$key} / $countT;
        # display 3 digits
        my $round = sprintf("%.3f",$norm);
        print "$key = $round\n";
}
exit;
