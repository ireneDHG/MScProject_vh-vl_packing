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
my $H = 3;      # High sd set
my $L = 1;      # Low sd set

my ($ref1,$ref2) = &CreateLowHighHashes($H,$L);


# Set output file names for high set
my $var="relTypeOutHp";
my $rVar="relTypeRHighSDp";
&CalculateTypesPrint($var,$rVar,$ref2,\@pdbFiles);

# Low set
my $var2="relTypeOutLp";
my $rVar2="relTypeRLowSDp";
&CalculateTypesPrint($var2,$rVar2,$ref1,\@pdbFiles);



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
                                        $highSD{$names[0]}=$1;
                                }
                                elsif($l > $1)
                                {
                                        $lowSD{$names[0]}=$1;
                                }                
                        }
                }
        }
        return(\%lowSD,\%highSD);
}

#********************************************************************
# Purpose: Classify the residues of a hash of arrays in types and print the results
#
# Arguments:
#       string $_[0]: name of txt file
#       string $_[1]: name of R file
#       string $_[2]: reference to the hash of arrays
#       string $_[3]: reference to an array of files 
#
# Requirements:
#       1. Precise 4 arguments
#       2. $_[0,1,2,3] should not be empty
#
# Return:
#       two output files
#
# Error if one of the arguments introduced is empty
#********************************************************************
sub CalculateTypesPrint
{
        my ($name1, $name2,$ref,$files) = @_;

        if($name1 eq "" or $name2 eq "" or $ref eq "" or $files eq "")
        {
                die "At least one of the arguments introduced in CalculateTypesPrint() is empty\n";
        }

        my %set = %$ref;
        my @pdbFiles = @$files;

        open(my $FI2,">$name1.txt") or die "Unable to open the file $name1\n";
        open(my $R2,">$name2.txt") or die "Unable to open the r file $name2\n";

        print $R2 "PDB\tSD\tHydrophobic\tHydrophilic\tAmphipathic\tGlycines\n";

        # For the low set of types, classify the relevant residues into types
        foreach my $name (keys %set)
        {
                foreach my $file(@pdbFiles)
                {
                        my $pdb = packingAngle::GetPdbFileName($file);
                        if($name eq $pdb)
                        {
                                my $sd=$set{$name};
                                my $hashref = ResiduesProp::GetAaRelevantVHVLCode($file);
                                ResiduesProp::SetResidueType($hashref,$pdb,$FI2,$R2,$sd);
                        }
                }
        }
        close($FI2);
        close($R2);
}
exit;
