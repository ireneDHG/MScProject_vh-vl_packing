package packingAngle;

#********************************************************************
# MSc Bioinformatics with Systems Biology 2015/16
# 
# Analysis of the VH/VL packing
#
# Author:       Irene del Hierro García
# Module name:  packingAngle.pm
#
# Description:  This module contains subroutines that are involved in the
#               calculation and analysis of the packing angle of the interfaces
#               of the antibodies
#
#********************************************************************

use config;     # Module to configure the paths
use sd;         # Module to calculate the sd and mean

#********************************************************************
# Purpose: Given an array reference of PDB files, calculate the packing angle
#          for each one
#
# Date: 12/08/16
# Version: V1.1
#
# Arguments:
#       string $_[0]: array reference of PDB files
#
# Requirements:
#       1. Precise 1 argument
#       2. $_[0] should be an array reference
#
# Return:
#       a hash reference with the PDB code as the key and the angle as the value
#
# The angles that could not be computed by the C program are avoided
#********************************************************************
sub CreateHashAnglePdb
{
        my ($reference) = @_;
        
        # Dereference the array
        my @pdbFiles = @$reference; 
        
        my %pdbRndHash=();

        # For each file in the array:
        foreach my $file (@pdbFiles)
        {	
                # Take out the ".pdb" extension
                my $pdb = GetPdbFileName($file);

                # Call the C program to get the angle
                my $angle = `cta.sh $config::antibodyDir/$file`;

                # If no angle is retrieved, do not update the hash
                if($angle != "")
                {
                        chomp $angle;
                        # Update the hash if an angle is retrieved and it hasn't been
                        # calculated before for the PDB file
                        $pdbRndHash{$pdb}=$angle unless exists $pdbRndHash{$pdb};
                }
        }
        # Return the hash reference
        return (\%pdbRndHash);
}

#********************************************************************
# Purpose: Given a PDB file retrieved the name without the extension
#
# Date: 12/08/16
# Version: V1.1
#          
# Arguments:
#       string $_[0]: PDB file
#
# Requirements:
#       1. Precise 1 argument
#       2. $_[0] should be a PDB file name
#
# Return:
#       the name of the file without the extension
#
# An error message is retrieved if a PDB file is not introduced.
#********************************************************************
sub GetPdbFileName
{
        my ($file)=@_;
        
        my $name = "";
        
        # Take the name of the file
        if($file =~ /(.*)\.pdb/)
        {
                $name = $1;
        }
        # If not a PDB file, give an error message
        else{
                die "The input $file introduced in packingAngle::GetPdbFileName is not a PDB file\n";
        }
        # Return the name of the file
        return ($name);
}


#********************************************************************
# Purpose: Given a hash reference, check the keys of the hash and print the result in 
#          an output file
#
# Date: 12/08/16
# Version: V1.1
#          
# Arguments:
#       string $_[0]: a hash reference
#
# Requirements:
#       1. Precise 2 arguments.
#       2. $_[0] should be a hash reference.
#       3. The keys of the hash should be the PDB codes and the values an angle value.
#       4. $_[1] is introduced in the terminal as an input file <>.
#       5. The input file is a txt file that contains a list of antibodies grouped by sequence
#          and separated by commas.
#
# Return:
#       Print the output in the screen if not output file is defined when calling the program
#
# An error message is retrieved if a PDB file is not introduced.
#********************************************************************
sub AddAngleToFile 
{
        my ($reference) = @_;
        
        my %hash = %$reference;
        
        # For every line of the input file:
        while (my $line = <>)
        {
        	chomp $line;
                my @names =();
                my @secondLine =();
        	
                # Store the names of the antibodies of that line in an array.
                @names = ($line =~ /(\d.{3}\_?\d*)/g);
        	
                my @abs = ();
                # For each antibody in that group, check if it exists in the input hash.
                # If exists, print the antibody separated by commas and add the value of the angle 
                # (obtained from the hash to an array that is going to be the next line).
                foreach my $element (@names)
                {		
                        foreach my $key (keys %hash)
                        {
                                if($element eq $key)
                                {
                                        push @abs, $element;
                                        push @secondLine, $hash{$key};
                                }
                        }
                }
                
                unless(scalar @abs == 0)
                {
                        foreach my $ab (@abs)
                        {
                                print "$ab, ";
                        }
                        print "\n";
                }
                	
                 unless(scalar @secondLine == 0)
                {
                        foreach my $angle (@secondLine)
                        {
                                print "$angle, ";
                        }
                        print "\n";
                }
        	
        }
}

#********************************************************************
# Purpose: Call the sd and mean calculator 
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
#       print the SD and mean in an output file (or the screen by default)
#       retrieve the sd
#
#********************************************************************
sub CallSDCalculator
{
        my ($refArray) =@_;
        
        my @angles = @$refArray;
        
        my $sd=""; my $mean="";
        
        # Call sd package to calculate either the mean and the sd
        if(@angles)
        {
                ($sd, $mean) = sd::CalcSD(\@angles);
                print "SD: $sd MEAN: $mean\n";		
        }
        else
        {
                @angles=(1);
                ($sd, $mean) = sd::CalcSD(\@angles);
                print "SD: $sd MEAN: $mean\n"; 
        } 
        
        return($sd); 
}

#********************************************************************
# Purpose: Create the txt file to allow the display of a distribution in R
#
# Date: 12/08/16
# Version: V1.1
#          
# Arguments:
#       string $_[0]: reference to the array of x values
#       string $_[1]: reference to the array of y values
#       string $_[2]: string with the name of the output file
#
# Return:
#       txt file
#
#********************************************************************
sub PrintFileToDistro 
{
        my($xval,$yval,$outputName)=@_;
        	
        # Create the output file
        open(OUT,">$outputName.txt") or die "Unable to open file";
        
        my @xArray = @$xval;
        my @yArray = @$yval;
        
        print OUT "GROUP\tSD\n";
        
        # Print in the output file the group number and the sd value in tabular format
        for(my $i=0;$i<scalar(@xArray);$i++)
        {
                for(my $j=0; $j<scalar(@yArray);$j++)
                {
                        if($i == $j and $yArray[$j] ne "undef")
                        {
                                print OUT "$xArray[$i]\t$yArray[$j]\n";
                        }
                }
        }
        
        close(OUT);
}	

#********************************************************************
# Purpose: Check if a group size is above or below a given threshold
#
# Date: 12/08/16
# Version: V1.1
#          
# Arguments:
#       string $_[0]: threshold applied
#       string $_[1]: reference of an array
#
# Requirements:
#       1. Precise 2 arguments
#       2. $_[0] should be numeric
#       3. $_[1] should be an array reference
#
# Return:
#       value = 1 if the size of the array is bigger or equal than the threshold
#       value = 0 if not
#
#********************************************************************	
sub Threshold
{
        my ($threshold,$arrayRef)=@_;
        
        my @array=@$arrayRef;
        my $nElem=scalar(@array);
        my $check=0;
        
        # Check if the threshold is numeric
        if($threshold =~ /\D+/)
        {
                print STDERR "The threshold introduced in packingAngle::Threshold is not numeric\n";
        }
        
        # Return 1 if it is equal or bigger than the threshold, or 0 if not
        if($nElem >= $threshold)
        {
                $check=1;
        }
        return($check);
}


1;

