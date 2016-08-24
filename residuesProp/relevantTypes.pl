#!usr/bin/perl -w
#********************************************************************
# MSc Bioinformatics with Systems Biology 2015/16
# 
# Analysis of the VH/VL packing
#
# Author:       Irene del Hierro García
# Script name:  relevantTypes.pl
# Version:      V1.1
# Date:         19/08/16
#
#********************************************************************
#
# Purpose
# =======
#        Classify the relevant residues in four types: hydrophobic, hydrophilic, amphipathic and glycines
#
#********************************************************************
#
# Assumptions
# ===========
#       1. The program requires an input file with a list of antibodies grouped by sequence.
#       2. The antibodies of the group should be separated by commas in the first line and in the next line 
#          the standard deviation and mean should appear as: "SD: 0.2344455 MEAN: 2.48".
#       3. The results will be printed in screen if an output file is not specified.
#       4. Change the sd split if necessary and the names of the output files
#
#       Usage: perl relevantTypes.pl < ../RedundantChothia_AllMeanSD.txt
#
#********************************************************************   
# Strategy
# ========
#       1. Open and read each line of the input file.
#       2. Split the file into two sets: a high sd set and a low sd set
#       3. For each set get a representantive of each redundant group and classify into types
#       4. Print the results in two output files
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
use ResiduesProp;
use packingAngle;
use config;

my %highSD=();
my %lowSD=();

# Create the two sets of standard deviations from the input file
while(my $line = <>)
{
        chomp $line;

        my @names = ($line =~ /(\d.{3}\_?\d*)/g);
	
        $line=<>;       
        
        # Extract the sd
        if($line =~ /^SD:\s(.+)\sMEAN:.*$/)
        {
                # Avoid those groups with undefined sd
                unless($1=~ /undef/){
                        if(3 < $1)
                        {
                                $highSD{$names[0]}=$1;
                        }
                        elsif(1 > $1)
                        {
                                $lowSD{$names[0]}=$1;
                        }
                }
        }
}


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

my $var="relTypeOutH";
my $rVar="reltypeRHighSD";

open(my $FI1,">$var.txt") or die "Unable to open the file 1\n";
open(my $R1,">$rVar.txt") or die "Unable to open the r file 1\n";

print $R1 "PDB\tSD\tHydrophobic\tHydrophilic\tAmphipathic\tGlycines\n";

# For the high set of types, classify the relevant residues into types
foreach my $name (keys %highSD)
{
        foreach my $file(@pdbFiles)
        {
                my $pdb = packingAngle::GetPdbFileName($file);
                if($name eq $pdb)
                {
                        my $sd=$highSD{$name};
                        my $hashref = ResiduesProp::GetAaRelevantVHVLCode($file);
                        ResiduesProp::SetResidueType($hashref,$pdb,$FI1,$R1,$sd);
                }
        }
}

my $var2="relTypeOutL";
my $rVar2="relTypeRLowSD";

open(my $FI2,">$var2.txt") or die "Unable to open the file 2\n";
open(my $R2,">$rVar2.txt") or die "Unable to open the r file 2\n";

print $R2 "PDB\tSD\tHydrophobic\tHydrophilic\tAmphipathic\tGlycines\n";

# For the low set of types, classify the relevant residues into types
foreach my $name (keys %lowSD)
{
        foreach my $file(@pdbFiles)
        {
                my $pdb = packingAngle::GetPdbFileName($file);
                if($name eq $pdb)
                {
                        my $sd=$lowSD{$name};
                        my $hashref = ResiduesProp::GetAaRelevantVHVLCode($file);
                        ResiduesProp::SetResidueType($hashref,$pdb,$FI2,$R2,$sd);
                }
        }
}

close($FI1);
close($FI2);
close($R1);
close($R2);

exit;
