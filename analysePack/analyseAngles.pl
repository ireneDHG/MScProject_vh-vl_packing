#!usr/bin/perl
#********************************************************************
# MSc Bioinformatics with Systems Biology 2015/16
# 
# Analysis of the VH/VL packing
#
# Author:       Irene del Hierro García
# Script name:  analyseAngles.pl
# Version:      V1.1
# Date:         12/08/16
#
#********************************************************************
#
# Purpose
# =======
#       Calculate the mean and standard deviation of each of the redundant groups in the  
#       RedundantChothia_Angles.txt file
#       In the case of groups formed by one element, the standard deviation is left as undefined.
#
#********************************************************************
#
# Assumptions
# ===========
#       1. The program requires an input txt file with the list of antbodies grouped by sequence and
#          their respective packing angle values below
#       2. The antibodies of the groups as well as the angles in the txt file should be separated by commas
#       3. The location of the modules with the subroutines used in this script should be located in the appropriate path
#          (change the location of the Modules folder to correctly used this program)
#       4. Unless a screen printing is required, add an output file to obtain the results.
#
#       perl analyseAngles.pl < ../RedundantChothia_Angles.txt > ../RedundantChothia_MeanSD.txt 
#
#********************************************************************   
# Strategy
# ========
#       1. Open and read each line of the input file
#       2. Get the angle values for each group
#       3. Calculate the mean and sd and print the results in an output file
#
#********************************************************************  
#
# Revision history
# ================ 
#       V1.0:   Original
#       V1.1:   Comments
#
#********************************************************************

use lib "/home/irene/Documents/MScProject/Modules";  

use strict;
use sd;                 # Package to calculate standard deviation and mean
use packingAngle;       # Package with the subroutines for this script

# Get the x (group number) and y vectors (sd values)
my ($xref,$yref) = &ProcessInputFile;

my @xVal = @$xref;
my @yVal = @$yref;

# Call a subroutine that it is simply going to create a diferent file that can be analyzed by R
my $nameOut = "sd_values";
packingAngle::PrintFileToDistro(\@xVal,\@yVal,$nameOut);




#********************************************************************
# Purpose: Process the input file and retrieve a vector with the number of the groups
#          and another vector with the values.
#
# Arguments:
#       none: the input is the input file 
#
# Requirements:
#       1. The file has to have specific format, i.e.:
#       5DFV_1, 5DFV_2, 
#       -48.838169, -45.589493, 
#
# Return:
#       array references of the group numbers and values
#
# If any line is empty, skip it
# If the file does not have that format, there is no guarantee of success
#********************************************************************
sub ProcessInputFile
{
        my $count=0;    # The counter would be the group name for the R file

        my @xValues=(); # To store the redundant group number
        my @yValues=(); # To store the sd values

        while (my $line = <>)
        {
                # Check if the line is empty
                unless($line =~ /^$/)
                {
                        chomp $line;
        
                        # Print the line with the PDB names in the output file
                        print "$line\n";
        
                        # Skip the line to get the line with the angle values
                        $line = <>;
        
                        my @angles=();

                        # Get the values of the angles of a single redundant group
                        @angles = ($line =~ /(\d+\.?\d*)/g);
                        $count++;
	
                        push @xValues, $count;
      	
                        # Call subroutine to retrieve the standard deviation and print the mean and sd in the output file
                        my $sd = packingAngle::CallSDCalculator (\@angles);
	
                        push @yValues, $sd;  
                }      
        }
        return(\@xValues,\@yValues);
}


exit;      
        
