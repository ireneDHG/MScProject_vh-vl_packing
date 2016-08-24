package sd;
#********************************************************************
# MSc Bioinformatics with Systems Biology 2015/16
# 
# Analysis of the VH/VL packing
#
# Author:       Extracted from ACRMGroup in Github and modified by Irene del Hierro García
# Module name:  sd.pm
#
# Description:  This module contains a sd calculation subroutine
#
# Date: 12/08/16
#
#********************************************************************

#********************************************************************
# Purpose: Calculate the mean and sd for a given group of values
#
# Date: 12/08/16
# Version: V1.1
#          
# Arguments:
#       string $_[0]: reference of array values
#
# Requirements:
#       1. Precise 1 argument
#       2. $_[0] should be an array reference of numeric values
#
# Return:
#       the sd and mean
#
#********************************************************************
sub CalcSD
{
        my($refValues) = @_;
        my @values = @$refValues;

        my $nvals = @values;
        my $sum = 0;

        # Calculate the mean
        for(my $i=0; $i<$nvals; $i++)
        {
                $sum += $values[$i];
        }
        my $mean = $sum / $nvals;

        # Calculate the sum of squares
        $sum = 0;
        my $sd;

        # Of a sample bigger than 1
        if ($nvals > 1)
        {
                for(my $i=0; $i<$nvals; $i++)
                {
                        $sum += (($values[$i] - $mean)*($values[$i] - $mean));
                }
        # And divide by n-1 then square root
        $sd = sqrt($sum/($nvals-1));
        }
    
        elsif ($nvals <= 1)
        {
        # As calculate the sd of one element does not make any sense, leave as undefined
                $sd = "undef";
        } 
        return($sd, $mean);
}


1;
