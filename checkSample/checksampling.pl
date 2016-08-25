#!usr/bin/perl -w

#********************************************************************
# MSc Bioinformatics with Systems Biology 2015/16
# 
# Analysis of the VH/VL packing
#
# Author:       Irene del Hierro García
# Script name:  checksampling.pl
# Version:      V1.1
# Date:         12/08/16
#
#********************************************************************
#
# Purpose
# ========
#       Calculate the standard deviation of each of the redundant groups in the  
#       RedundantChothia_Angles.txt file applying three different group thresholds, i.e. the first
#       threshold is 3 which means that the sd is only going to be obtained for groups with three or more 
#       structures. 
#       The thresholds applied are three, four and five.
#       The output files obtained with this script are in tabular format so they can be accessed by an R
#       script that is going to determine the distribution of the standard deviations obtained and therefore,
#       the sampling of the groups.
#
#********************************************************************
#
# Assumptions
# ===========
#       1. The program requires an input file with a list of antibodies grouped by sequence.
#       2. The antibodies of the group should be separated by commas in the first line and in the next line 
#          their respective angles
#       3. The results will be printed in screen if an output file is not specified.
#       4. Changes in the location of the modules should be applied by updating the path un the "use lib" statement
#
#       Usage: perl checksampling.pl < ../RedundantChothia_AllAngles.txt
#
#********************************************************************   
# Strategy
# ========
#       1. Open and read each line of the input file.
#       2. Compute the group size thresholds
#       3. Calculate the sd
#       4. Create an output tabular file for each threshold
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
use packingAngle;       # Package that contains the subroutines for all the angle analysis
use sd;	                # Package that contains the subroutines for the sd calculation

# Set the name of the output files and the group size threshold
my $name3 = "sd_AllvaluesThres_3"; 
my $three = 3;
my $name4 = "sd_AllvaluesThres_4"; 
my $four = 4;
my $name5 = "sd_AllvaluesThres_5"; 
my $five = 5;

my ($x3,$y3,$x4,$y4,$x5,$y5) = ProcessFile($three, $four, $five);

my @xValues3=@$x3; my @yValues3=@$y3;       
my @xValues4=@$x4; my @yValues4=@$y4;
my @xValues5=@$x5; my @yValues5=@$y5;

# Calling to subroutines to print the three tabular output files
packingAngle::PrintFileToDistro(\@xValues3,\@yValues3,$name3);
packingAngle::PrintFileToDistro(\@xValues4,\@yValues4,$name4);
packingAngle::PrintFileToDistro(\@xValues5,\@yValues5,$name5);
       
#********************************************************************
# Purpose: Process the file and add three group size thresholds to it
#
# Arguments:
#       string $_[0],[1],[2]: threshold values
#
# Requirements:
#       1. Precise 3 arguments
#       2. All should be numeric and non empty values
#
# Return:
#       three pairs of arrays (each pair for each threshold: an x array with the number of the group and an
#       y array with the sd values.
#
# Error retrieved if value introduced is empty or non numeric
#********************************************************************
sub ProcessFile
{
        my ($three, $four, $five)=@_;

        if($three eq "" or $four eq "" or $five eq "")
        {
                die "One or more of the thresholds introduced in &ProcessFile are empty\n";
        }
        unless($three =~ /\d+/ or $four =~ /\d+/ or $five =~ /\d+/)
        {
                die "One or more of the thresholds introduced in &ProcessFile are non numeric\n";
        }

        # The counters are going to be the groups number 
        # The "x" arrays are going to store the group numbers and the "y" arrays are going to store the sd values
        my $count3=0;
        my @xValues3=(); my @yValues3=();       

        my $count4=0;
        my @xValues4=(); my @yValues4=();

        my $count5=0;
        my @xValues5=(); my @yValues5=();

        while (my $line = <>)
        {
                # Check if the line is empty
                unless($line =~ /^$/)
                {

                        chomp $line;
        
                        # Skip the line to get the line with the angle values
                        $line =<>;

                        my @angles=();
                        # Get the values of the angles of a single redundant group
                        @angles = ($line =~ /(\d+\.?\d*)/g);
                        
                        # Call a subroutine to apply the threshold
                        my $sd3 = &ApplyThreshold($three,\@angles);
                        # If the sd retrieved is not 0, then add the values to the arrays
                        if($sd3 ne "")
                        {
                        $count3++;
                        push @xValues3,$count3;
                        push @yValues3,$sd3; 
                        }
                        
                        my $sd4 = &ApplyThreshold($four,\@angles);
                        unless($sd4 eq "")
                        {
                        $count4 ++;
                        push @xValues4,$count4;
                        push @yValues4,$sd4;
                        }
        
                        my $sd5 = &ApplyThreshold($five,\@angles);
                        unless($sd5 eq "")
                        {
                        $count5 ++;
                        push @xValues5,$count5;
                        push @yValues5,$sd5;
                        }
                }
        }
        return(\@xValues3,\@yValues3,\@xValues4,\@yValues4,\@xValues5,\@yValues5);
}

#********************************************************************
# Purpose: Check if a group fits a threshold and retrieve sd if it does
#
# Arguments:
#       string $_[0]: threshold value
#       string $_[1]: array reference of angle values
#
# Requirements:
#       1. Precise 2 arguments
#       2. Threshold should be numeric and non empty
#       3. The reference should be an array reference of angle values, non empty
#
# Return:
#       the sd as empty if the group does not fit, as numeric if it fits
#
# Error retrieved if values introduced are empty
# Error retrieved if threshold is non numeric
#********************************************************************
sub ApplyThreshold 
{
        my ($thres,$angle) =@_; 
        
        if($thres eq "" or $angle eq "")
        {
                die "One of the arguments introduced in &ApplyThreshold is empty\n";
        }
        unless($thres =~ /\d+/)
        {
                die "The threshold introduced in &ApplyThreshold is non numeric\n";
        }

        my @angles = @$angle; my $sd="";

        # Call a subroutine that is going to set the threshold for each group
        # In this case, for the threshold of three groups
        # If the output is 1, then the redundant group fits the threshold, if 0, it doesn't.
        my $check= packingAngle::Threshold($thres,\@angles);
        # Add the sd values if the group fits the thereshold
        if($check == 1)
        {
                # Call subroutine to retrieve the standard deviation and print the mean and sd in the output file
                $sd = packingAngle::CallSDCalculator(\@angles);
         }
        return($sd);
}

	
exit;
