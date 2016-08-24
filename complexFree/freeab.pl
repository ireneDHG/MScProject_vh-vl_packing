#!usr/bin/perl -w
#********************************************************************
# MSc Bioinformatics with Systems Biology 2015/16
# 
# Analysis of the VH/VL packing
#
# Author:       Irene del Hierro García
# Script name:  freeab.pl
# Version:      V1.1
# Date:         19/08/16
#
#********************************************************************
#
# Purpose
# =======
#        Classify the redundant groups into complex, free or mixed if their antibodies
#        are all complexed, in free-state or both complexed and free respectively, and just give
#        a list of the free redundant groups with their respective standard deviations.
#
#********************************************************************
#
# Assumptions
# ===========
#       1. The program requires an input txt file with a list of antibodies grouped by sequence
#          and an input folder that contains the structures of antibodies (as PDB files).
#       2. The antibodies of the group should be separated by commas in the first line and in the next line 
#          the standard deviation and mean should appear as: "SD: 0.2344455 MEAN: 2.48".
#       3. The results will be printed in screen if an output file is not specified.
#
#       Usage: perl freeab.pl < ../RedundantChothia_AllMeanSD.txt > outputCF.txt
#
#********************************************************************   
# Strategy
# ========
#       1. Open and read each line of the input file.
#       2. Get the names of the antibodies and open their respective PDB files
#       3. Determine if the antibodies are in free-state or bounded with a ligand
#       4. Determine the class of the redundant group when all its antibodies have been analysed.
#       5. Print the results in screen or in file if specified.
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

use strict;
use config;
use CFtest;
use packingAngle;

# Get the directory of the input folder
my $dir;
if(!defined($dir))
{
        $dir = $config::propdir;
}

# Get the pdb files from the directory
opendir (DIR, $dir) or die "Unable to open dir $dir\n";
my @pdbFiles = grep(/.*\.pdb$/,readdir(DIR));
closedir(DIR);

my $id = 0; my %hashCF=();

# For each line of the input file,
# get the names of the antibodies, the sd and the mean
while(my $line = <>)
{
        chomp $line;

        my @names = ($line =~ /(\d.{3}\_?\d*)/g);
        $id++;
	
        # Skip line to get the sd and mean
        $line = <>; 
        my $sd = CFtest::GetSDFromFile($line);
        my $mean = CFtest::GetMeanFromFile($line);
	
        # Avoid those groups with undefined sd
        unless($sd eq "undef")
        {	
                my $complex=0; my $free =0;
                foreach my $str (@names)
                {
                        foreach my $name (@pdbFiles)
                        {
                                # For each antibody in the group, check its pdb file
                                # Get the name of the file. If it is the same, look for complexity
                                my $pdb = packingAngle::GetPdbFileName($name);
                                if($str eq $pdb)
                                {
                                        # Update the counters depending on the outcome (above 2 means complexed)
                                        my $number = CFtest::AssignCF($name);
                                        if ($number > 2)
                                        {
                                                $complex++;
                                        }
                                        else
                                        {
                                                $free++;
                                        }
                                }
                        }
                }
	        # Once all antibodies of the group have been analysed, check the class of the group
                my $value = CFtest::CheckCFGroup($complex,$free);
	        # Assign a class depending on the outcome and print just the free redundant groups
                if($value == 1)
                {
                        foreach my $elem (@names)
                        {
                                print "$elem, ";
                        }
                        print "\nSD: $sd MEAN: $mean\n";
                }

        }
}

exit;
