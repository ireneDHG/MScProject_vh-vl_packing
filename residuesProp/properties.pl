#!usr/bin/perl -w
#********************************************************************
# MSc Bioinformatics with Systems Biology 2015/16
# 
# Analysis of the VH/VL packing
#
# Author:       Irene del Hierro García
# Script name:  properties.pl
# Version:      V1.1
# Date:         19/08/16
#
#********************************************************************
#
# Purpose
# =======
#        Calculate the average hydrophobicity, average isoelectric point and molecular weight of one element of
#        each redundant group
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
#       Usage: perl properties.pl < ../RedundantChothia_AllMeanSD.txt
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

# Define the sd split:
my $H = 4;      # High sd set
my $L = 1;      # Low sd set

my ($ref1,$ref2) = &CreateLowHighHashes($H,$L);

# For the Low SD set:
my $nameL = "freediv6propL";
my $hashrefL = &CalculateProp($ref1,\@pdbFiles);
&PrintProp($nameL,$hashrefL);

#For the High SD set
my $nameH = "freediv6propH";
my $hashrefH = &CalculateProp($ref2,\@pdbFiles);
&PrintProp($nameH,$hashrefH);


#********************************************************************
# Purpose: Classify the groups in a low sd hash and a high sd hash
#
# Arguments:
#       string $_[0]: high threshold value
#       string $_[1]: low threshold value
#       input file introduced in terminal 
#
# Requirements:
#       1. Precise 2 arguments
#       2. $_[0] should be a number and non empty
#       3. $_[0] should be a number and non empty
#
# Return:
#       a reference to a low sd hash and a reference to a high sd hash
#
# Give error message if thresholds introduced are empty or non digit
#********************************************************************
sub CreateLowHighHashes
{
        my ($h,$l) = @_;

        if($h eq "" or $l eq "")
        {
                die "H and L thresholds are empty\n";
        }
        unless($h =~ /\d+/ or $l =~ /\d+/)      
        {
                die "H and L thresholds are not numeric\n";
        }

        my %lowSD =();
        my %highSD = ();

        # For each line of the file, get the names of the proteins and the sd
        while (my $line =<>)
        {
                chomp $line;
        
                my @names = ($line =~ /(\d.{3}\_?\d*)/g);
        
                $line = <>;
        
                if($line =~ /^SD:\s(.+)\sMEAN:.*$/)
                {
                        # unless the sd is undefined, add to the high or low sd set
                        unless($1=~ /undef/)
                        {
                                if($h < $1)
                                {
                                        push @{$highSD{$names[0]}},$1;
                                }
                                elsif($l > $1)
                                {
                                        push @{$lowSD{$names[0]}},$1;
                                }                
                        }
                }
        }
        return(\%lowSD,\%highSD);
}

#********************************************************************
# Purpose: Calculate properties for the elements of a hash of arrays and
#          add them to the hash
#
# Arguments:
#       string $_[0]: reference to a hash of arrays
#       string $_[1]: reference to an array 
#
# Requirements:
#       1. Precise 2 arguments
#       2. $_[0] should be a hash of arrays with the name of the PDB file as key and the sd as first and only value of the array
#       3. $_[0] should be an array of PDB files
#
# Return:
#       a reference to the hash of arrays updated
#
#********************************************************************
sub CalculateProp
{
        my($ref1,$ref2) =@_;
        
        my %setSD = %$ref1;

        my @pdbFiles = @$ref2;

        # Calculate the properties for the high sd set
        foreach my $name (keys %setSD)
        {
                foreach my $file(@pdbFiles)
                {
                        # Get the name of the file
                        my $pdb = packingAngle::GetPdbFileName($file);
                        # If it is equal to the name of the high sd set
                        if($name eq $pdb)
                        {
                                # Get the key determining residues of the interface
                                my $hashref = ResiduesProp::GetAaInterfaceVHVLCode($file);
                                # Calculate properties
                                my $hydro = ResiduesProp::CalculateHydropAVG($hashref);
                                my $mw = ResiduesProp::CalculateMW($hashref);
                                my $pi = ResiduesProp::CalculatePIAVG($hashref);
			
                                # Add the results to the high sd hash
                                push @{$setSD{$name}},$hydro,$mw,$pi;
                        }
                }
        }
        return(\%setSD);
}

#********************************************************************
# Purpose: Print a hash of arrays in a file
#
# Arguments:
#       string $_[0]: name of the file
#       string $_[1]: reference to the hash
#
# Requirements:
#       1. Precise 2 arguments
#       2. $_[0] should be non empty
#       3. $_[0] should be non empty
#
# Return:
#       a reference to a low sd hash and a reference to a high sd hash
#
# Give error message if one of the arguments is empty
#********************************************************************
sub PrintProp
{
        my($name, $ref) =@_;

        if($name eq "" or $ref eq "")
        {
                die "At least one of the arguments introduced in &PrintProp() is empty\n";
        }

        my %set = %$ref;

        open(OUT,">$name.txt") or die "Unable to open the output $name file\n";

        print OUT "HIGHpdb\tSD\tHYDRO\tMW\tPI\n";

        foreach my $key (keys %set)
        {
                print OUT "$key";
                foreach my $elem (@{$set{$key}})
                {
                        print OUT "\t$elem";
                }
                print OUT "\n";
        }

        close(OUT);
}


exit;
