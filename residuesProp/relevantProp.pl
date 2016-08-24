#!usr/bin/perl -w
#********************************************************************
# MSc Bioinformatics with Systems Biology 2015/16
# 
# Analysis of the VH/VL packing
#
# Author:       Irene del Hierro García
# Script name:  relevantProp.pl
# Version:      V1.1
# Date:         19/08/16
#
#********************************************************************
#
# Purpose
# =======
#        Calculate the average hydrophobicity, average isoelectric point and molecular weight of one element of
#        each redundant group and just for the relevant key determining residues
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
#       Usage: perl relevantProp.pl < ../RedundantChothia_AllMeanSD.txt
#
#********************************************************************   
# Strategy
# ========
#       1. Open and read each line of the input file.
#       2. Split the file into two sets: a high sd set and a low sd set
#       3. For each set get a representantive of each redundant group and calculate the properties
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
                        # Set the division set
                        if(4 < $1)
                        {
                                push @{$highSD{$names[0]}},$1;
                        }
                        elsif(2 > $1)
                        {
                                push @{$lowSD{$names[0]}},$1;
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

# Calculate the properties for the high sd set
foreach my $name (keys %highSD)
{
        foreach my $file(@pdbFiles)
        {
                # Get the name of the file
                my $pdb = packingAngle::GetPdbFileName($file);
                # If it is equal to the name of the high sd set
                if($name eq $pdb)
                {
                        # Get the relevant key determining residues of the interface
                        my $hashref = ResiduesProp::GetAaRelevantVHVLCode($file);
                        # Calculate properties
                        my $hydro = ResiduesProp::CalculateHydropAVG($hashref);
                        my $mw = ResiduesProp::CalculateMW($hashref);
                        my $pi = ResiduesProp::CalculatePIAVG($hashref);
			# Add the results to the high sd hash
                        push @{$highSD{$name}},$hydro,$mw,$pi;
                }
        }
}

# Calculate the properties for the low sd set
foreach my $name (keys %lowSD)
{
        foreach my $file(@pdbFiles)
        {
                my $pdb = packingAngle::GetPdbFileName($file);
                if($name eq $pdb)
                {
                        my $hashref = ResiduesProp::GetAaRelevantVHVLCode($file);
                        my $hydro = ResiduesProp::CalculateHydropAVG($hashref);
                        my $mw = ResiduesProp::CalculateMW($hashref);
                        my $pi = ResiduesProp::CalculatePIAVG($hashref);
                        
                        push @{$lowSD{$name}},$hydro,$mw,$pi;
                }
        }
}

# Change the names of these files if necessary

open(OUT,">relevantH2.txt") or die "Unable to open the output file\n";

print OUT "HIGHpdb\tSD\tHYDRO\tMW\tPI\n";

foreach my $key (keys %highSD)
{
        print OUT "$key";
        foreach my $elem (@{$highSD{$key}})
        {
                print OUT "\t$elem";
        }
        print OUT "\n";
}

close(OUT);

open(LOW,">relevantL2.txt") or die "Unable to open the output file\n";

print LOW "LOWpdb\tSD\tHYDRO\tMW\tPI\n";

foreach my $key (keys %lowSD)
{
        print LOW "$key";
        foreach my $elem (@{$lowSD{$key}}){
                print LOW "\t$elem";
        }
        print LOW "\n";
}

close(LOW);

exit;
