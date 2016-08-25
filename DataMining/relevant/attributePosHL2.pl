#!usr/bin/perl -w
#********************************************************************
# MSc Bioinformatics with Systems Biology 2015/16
# 
# Analysis of the VH/VL packing
#
# Author:       Irene del Hierro García
# Script name:  attributePosHL2.pl
# Version:      V1.1
# Date:         19/08/16
#
#********************************************************************
#
# Purpose
# =======
#        Create a CSV file with the residues in the relevant key determing positions and a nominal class output (input 3).
#
#********************************************************************
#
# Assumptions
# ===========
#       1. The program requires an input file with a list of antibodies grouped by sequence.
#       2. The antibodies of the group should be separated by commas in the first line and in the next line 
#          the standard deviation and mean should appear as: "SD: 0.2344455 MEAN: 2.48".
#
#       Usage: perl attributePosHL2.pl < ../RedundantChothia_MeanSD.txt
#
#********************************************************************   
# Strategy
# ========
#       1. Open and read each line of the input file.
#       2. Get the sd and the interface residues for the file chosen
#       3. Divide between high and low sd to get the classes H and L
#       4. Add the residues to the positions
#       5. Print the results in an output file
#
#********************************************************************
#
# Revision history
# ================ 
#       V1.0:   Original
#       V1.1:   Comments
#
#********************************************************************

use lib "/home/irene/Documents/MScProject/Modules";     #location of the modules


use strict;
use totalProt;
use packingAngle;
use ResiduesProp;
use config;

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

#check the undef
my $name = "relevantPosHL2";
my @pos = ("L41","L87","H105","L38","H60","H45");
my $position = join(",",@pos);
open(OUT,">$name.csv") or die "Unable to open the output file $!\n";  
print OUT "Protein,SD,$position\n";

# For every line of the file, take the names of the group and the sd
while(my $line = <>)
{
        chomp $line;
        my $pdb="";
        my @names = ($line =~ /(\d.{3}\_?\d*)/g);
        my $name=$names[0];
        $line = <>;

        foreach my $file(@pdbFiles)
        {
                $pdb= packingAngle::GetPdbFileName($file);
                if($name eq $pdb)
                {
                         if($line =~ /^SD:\s(.+)\sMEAN:.*$/)
                         {
                                # If the sd is not undefined, get the residues from the interface and add them to their positions
                                unless($1 =~ /undef/)
                                {
                                        my $sd = $1;
                                        my $value="";
                                        # Create the split in the sd
                                        if($sd > 4)
                                        {
                                                $value = "H";
                                                my $hashref= ResiduesProp::GetAaInterfaceVHVLCode($file);
                                                my %aaHash = %$hashref;
                                                print OUT "$pdb,$value";
                                                foreach my $k (@pos)
                                                {
                                                        foreach my $key (keys %aaHash)
                                                        {
                                                                my $res = @{$aaHash{$key}}[1].@{$aaHash{$key}}[2];
                                                                if($k eq $res)
                                                                {
                                                                        print OUT ",",@{$aaHash{$key}}[0];
                                                                }
                                                        }
                                                }
                                                print OUT "\n";
                                        }
                                        elsif($sd < 2)
                                        {
                                                $value = "L";
                                                my $hashref= ResiduesProp::GetAaInterfaceVHVLCode($file);
                                                my %aaHash = %$hashref;
                                                print OUT "$pdb,$value";
                                                foreach my $k (@pos)
                                                {
                                                        foreach my $key (keys %aaHash)
                                                        {
                                                                my $res = @{$aaHash{$key}}[1].@{$aaHash{$key}}[2];
                                                                if($k eq $res)
                                                                {
                                                                        print OUT ",",@{$aaHash{$key}}[0];
                                                                }
                                                        }
                                                }
                                                print OUT "\n";
                                        }
                                        
                                }
                          }
                          
                }
        }
}
      
close(OUT);

exit;



