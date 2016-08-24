#!usr/bin/perl -w
#********************************************************************
# MSc Bioinformatics with Systems Biology 2015/16
# 
# Analysis of the VH/VL packing
#
# Author:       Irene del Hierro García
# Script name:  sdmatrix.pl
# Version:      V1.1
# Date:         19/08/16
#
#********************************************************************
#
# Purpose
# =======
#        Divide the mixed, free and complex redundant groups in a high spread set and a low spread set.
#        The output is given in matrix format.
#
#********************************************************************
#
# Assumptions
# ===========
#       1. The program requires an input txt file with a list of mixed, free or redundant groups with their sd 
#          (ideally this file is the output of the script complexFree.pl).
#       3. If another sd division is required, change the split and name of the output file
#
#       Usage: perl sdmatrix.pl < outputCF.txt > matrix.txt
#
#********************************************************************   
# Strategy
# ========
#       1. Open and read each line of the input file.
#       2. Get the names of the antibodies and open their respective PDB files
#       3. Determine if the antibodies are in free-state or bounded with a ligand
#       4. Determine the class of the redundant group when all its antibodies have been analysed.
#       5. Print the result in screen or in file if specified.
#
#********************************************************************
#
# Revision history
# ================ 
#       V1.0:   Original
#       V1.1:   Comments
#
#********************************************************************

use lib "/home/irene/Documents/MScProject/Modules";     #location of the modules (if not want to use just comment it).

# Divide the result of complexfree.pl into high and low sd sets and create a matrix of relative abundances
# Call as: perl sdmatrix.pl < outputCF.txt

use strict;


my $complexL=0; my $freeL=0; my $mixL=0;
my $complexH=0; my $freeH=0; my $mixH=0;
my $countH=0; my $countL = 0;

# For each line of the file
while(my $line = <>)
{
        chomp $line;
        # Store the class of the group
        my $variable = $line;
        # SKip the line and stored the sd
        my $line = <>; my $sd;
        
        if($line=~ /^(\d+.?\d+)$/)
        {
                $sd = $1;
        }
	
        # Make the sd division
        if($sd <= 2)
        {
                # Update the counters depending on the class of the group
                if($variable =~ /Free set/)
                {
                        $freeL++;
                }
                if($variable =~ /Complex set/)
                {
                        $complexL++;
                }
                if($variable =~ /Mixed set/)
                {
                        $mixL++;
                }

        }
        if($sd > 2)
        {

                if($variable =~ /Free set/)
                {
                        $freeH++;
                }
                if($variable =~ /Complex set/)
                {
                        $complexH++;
                }
                if($variable =~ /Mixed set/)
                {
                        $mixH++;
                }
        }
}

# Print the counters in matrix format
print "SD\tF\tC\tF/C\n";
print "Low\t$freeL\t$complexL\t$mixL\n";
print "High\t$freeH\t$complexH\t$mixH\n";


exit;
