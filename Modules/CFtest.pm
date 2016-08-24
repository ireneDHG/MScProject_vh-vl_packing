package CFtest;
#********************************************************************
# MSc Bioinformatics with Systems Biology 2015/16
# 
# Analysis of the VH/VL packing
#
# Author:       Irene del Hierro García
# Module name:  CFtest.pm
#
# Description:  This module contains subroutines that are involved in determining whether
#               antibodies or redundant groups are complexed or in free-state
#
#********************************************************************

use config;

# Define the path of the folder with the PDB files
my $dir;
if(!defined($dir))
{
        $dir = $config::propdir;
}

#********************************************************************
# Purpose: Check if a redundant group is composed only for free, complexed or both 
#          types of antibodies
#
# Date: 12/08/16
# Version: V1.1
#
# Arguments:
#       string $_[0]: number of complexed antibodies
#       string $_[1]: number of free-state antibodies
#
# Requirements:
#       1. Precise 2 arguments
#       2. Both should be numeric
#
# Return:
#       Value = 1 means the group is Complex
#       Value = 2 means the group is Free
#       Value = 3 means the group is Mixed
#
#********************************************************************
sub CheckCFGroup
{
        my ($valueC,$valueF)= @_;

        my $finalVal = 0;

        # Assign the class by determining the number of free or complexed antibodies
        if($valueC == 0)
        {
                $finalVal = 1;
        }
        elsif($valueF == 0)
        {
                $finalVal = 2;
        }
        elsif($valueC >0 and $valueF > 0)
        {
                $finalVal = 3;
        }

        return ($finalVal);
}

#********************************************************************
# Purpose: Determine if a given antibody is in free or complexed state,
#           i.e., bound or not with a ligand
#
# Date: 12/08/16
# Version: V1.1
#
# Arguments:
#       string $_[0]: pdb file
#
# Requirements:
#       1. Precise 1 argument
#       2. $_[0] should be a PDB file that is in the folder defined by config.pm
#
# Return:
#       A value as a number of chains found in the file (if it is more than two
#       means that there is a ligand. Otherwise, the antibody is free.
#
# Error message if the pdb file cannot be opened.
#********************************************************************
sub AssignCF
{
        my ($file)=@_;
	
        # Open the file from the directory defined above
        open(IN,"$dir/$file") or die "Unable to open the pdb file\n";
	
        my $count = 0;

        # Check the number of chains the protein has and retrieve that number
        while(my $line = <IN>)
        {
	        if($line=~/^REMARK\s+\d+\s+CHAIN\s+\w\s+.*$/)
	        {
                        $count++;	
                }
        }
        close(IN);
        return($count);
}

#********************************************************************
# Purpose: Get the standard deviation of a string with a specific format
#
# Date: 12/08/16
# Version: V1.1
#
# Arguments:
#       string $_[0]: a line of a file
#
# Requirements:
#       1. Precise 1 argument
#       2. $_[0] should be a string with this format: SD: 0.3454354 MEAN: 0.789
#
# Return:
#       A string as the standard deviation
#
#********************************************************************
sub GetSDFromFile
{
        my($line)=@_;

        my $sd=0;
	
        if($line =~ /^SD:\s(.+)\sMEAN:.*$/)
        {
                $sd=$1;
        }
	
        return($sd);
}

#********************************************************************
# Purpose: Get the mean of a string with a specific format
#
# Date: 12/08/16
# Version: V1.1
#
# Arguments:
#       string $_[0]: a line of a file
#
# Requirements:
#       1. Precise 1 argument
#       2. $_[0] should be a string with this format: SD: 0.3454354 MEAN: 0.789
#
# Return:
#       A string as the mean
#
#********************************************************************
sub GetMeanFromFile
{
        my($line)=@_;

        my $mean=0;
	
        if($line =~ /^SD:\s.+\sMEAN:(.*)$/)
        {
                $mean=$1;
        }
	
        return($mean);
}   

1;
