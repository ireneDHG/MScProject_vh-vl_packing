#!usr/bin/perl -w
#********************************************************************
# MSc Bioinformatics with Systems Biology 2015/16
# 
# Analysis of the VH/VL packing
#
# Author:       Irene del Hierro García
# Script name:  checkpdbLowSD.pl
# Version:      V1.0
# Date:         19/08/16
#
#********************************************************************
#
# Purpose
# =======
#        Identify if the proteins belonging to the same redundant group with 
#        low standard deviation come from the same pdb code or not, i.e. are
#        several files were obtained from the original one.
#        The value of the sd studied can be changed.
#        The output file is tabular TXT format that can be analysed by a t-test.
#
#********************************************************************
#
# Assumptions
# ===========
#       1. The program requires an input file with a list of antibodies grouped by sequence.
#       2. The antibodies of the group should be separated by commas in the first line and in the next line 
#          the standard deviation and mean should appear as: "SD: 0.2344455 MEAN: 2.48".
#       3. The results will be printed in screen if an output file is not specified.
#
#       Usage: perl checkpdbLowSD.pl < RedundantChothia_MeanSD.txt > samePDB.txt
#
#       To change the standard deviation threshold applied, just edit this file (line 53)
#
#********************************************************************   
# Strategy
# ========
#       1. Open and read each line of the input file.
#       2. Take the names of the antibodies of the redundant group and extract its sd.
#       3. If sd lower than a given value, extract the pdb codes of the antibodies and compare them
#       4. If all the pdb codes are the same, classify group as "Same", else as "Different"
#       5. Print class and sd in tabular format 
#
#********************************************************************
#
# Revision history
# ================ 
#       V1.0:   Original
#
#********************************************************************

use strict;

# Set SD threshold and call a subroutine to perform the analysis
my $threshold = 1;

&checkSameDiffGroups($threshold);


#********************************************************************
# Purpose: Given a sd threshold, identify the groups below that threshold and
#          determine whether they come from the same antibody or not.
#          Classify as "Same" if they come from the same antibody or as 
#          "Different" if not.
#
# Arguments:
#       string $_[0]: threshold value 
#
# Requirements:
#       1. Precise 1 argument
#       2. $_[0] should be a number
#       3. $_[0] is not empty string
#
# Return:
#       print the output in tabular format where the first column is the class of the 
#       group and the second the sd values
#
# Give error message if string introduced is empty or non digit
#********************************************************************
sub checkSameDiffGroups 
{
        my ($sdThreshold) = @_;

        if($sdThreshold eq "" or $sdThreshold =~ /\D+/)
        {
                print STDERR "Error: sd threshold introduced is not a digit or it's empty\n";
        }
        
        # Print the header of the output file
        print "PDBs\tSD\n";

        # For each line of the file
        while(my $line = <>)
        {
                chomp $line;
        
                # Take the names of the proteins
                my @names = ($line =~ /(\d.{3}\_?\d*)/g);
        
                # Skip the line to get the standard deviation
                $line = <>;       

                if($line =~ /^SD:\s(.+)\sMEAN:.*$/)
                {
                        # Just get the standard deviations that are not set as undefined and that are
                        # lower than a given value i.e., 2.
                        unless($1=~ /undef/){

                                if($sdThreshold > $1)
                                {
                                        my $sd = $1;

                                        # For each name of the proteins, extract the 4-char pdb code
                                        # and store in a hash
                                        my %no_seen=();
                                        foreach my $elem (@names)
                                        {
                                                my $pdb = GetPdbCode($elem);
                                        
                                                $no_seen{$pdb}++;
                                        }

                                        my @key = keys %no_seen;

                                        # If the hash has more than one key, means that the redundant group is
                                        # composed by different antibodies, if not, they are all from the same 
                                        # antibody
                                        if(scalar @key == 1)
                                        {
                                                print "Same\t$sd\n";
                                        }
                                        else 
                                        {
                                                print "Different\t$sd\n";
                                        }
                                }
                        }
                }
        }
}




#********************************************************************
# Purpose: Given a protein name, return the 4-character pdb code
# Arguments:
#       string $_[0]: name of the protein
#
# Requirements:
#       1. Precise 1 argument
#       2. $_[0] should be a alphanumerical string
#       3. $_[0] is not empty string
#
# Return:
#       string  a 4 character string
#
# Give error message if string introduced is empty
#********************************************************************
sub GetPdbCode 
{
        my ($name) = @_;
        
        my $pdb="";
        
        if($name eq "")
        {
                die "The string introduced in &GetPdbCode is empty\n";
        }
        else
        {
                if($name =~ /(\w{4})_\d/)
                {
                        $pdb = $1;
                }
        }
        return($pdb);
}

exit;
