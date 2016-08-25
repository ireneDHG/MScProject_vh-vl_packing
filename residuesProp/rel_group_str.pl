#!usr/bin/perl -w
#********************************************************************
# MSc Bioinformatics with Systems Biology 2015/16
# 
# Analysis of the VH/VL packing
#
# Author:       Irene del Hierro García
# Script name:  rel_group_str.pl
# Version:      V1.1
# Date:         19/08/16
#
#********************************************************************
#
# Purpose
# =======
#        Display the residues of determined positions for a group of structures.
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
#       Usage: perl rel_group_str.pl < ../freeabSD.txt
#
#********************************************************************   
# Strategy
# ========
#       1. Open and read each line of the input file.
#       2. Split the file into two sets: a high sd set and a low sd set
#       3. For each set get the residues and print in a file
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
use config;
use packingAngle;

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
my $H = 5;      # High sd set
my $L = 2;      # Low sd set

my ($ref1,$ref2) = &CreateLowHighHashes($H,$L);

# Define the names of the output files
my $lowset="rel_Str_L4";
my $highset = "rel_Str_H4";

# Create the low SD file
&GetResiduesAndPrint($lowset,$ref1);
# Create the high SD file
&GetResiduesAndPrint($highset,$ref2);


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
# Purpose: Get the relevant residues from a hash reference and print in an output file
#
# Arguments:
#       string $_[0]: name of the output file
#       string $_[1]: hash reference
#
# Requirements:
#       1. Precise 2 arguments
#       2. $_[0] should be non empty
#       3. $_[0] should be non empty
#
# Return:
#       an output file
#
# Give error message if thresholds introduced are empty or non digit
#********************************************************************
sub GetResiduesAndPrint
{
        my ($name, $ref) = @_;
        my %setSD = %$ref;

        if($name eq "" or $ref eq "")           
        {
                die "At least one of the arguments introduced in &GetResiduesAndPrint() is empty\n";
        }
        
        open(OUT,">$name.txt") or die "Unable to open the file $name $!\n";

        my %residues = (); my %positions = ();

        # For each of the proteins in the set, open its PDB file and get the residues of the interface
        foreach my $name (keys %setSD)
        {
                foreach my $file (@pdbFiles)
                {
                        my $pdb = packingAngle::GetPdbFileName($file);
                        if($name eq $pdb)
                        {

                                my $sd=$setSD{$name};
                                my $hashref = ResiduesProp::GetAaRelevantVHVLCode($file);
                                my %hash = %$hashref;
                                my %res = ();
                                # From the hash, take the position i.e., L38 and the residue type
                                foreach my $key (keys %hash)
                                {
                                        my $relPos = "";
                                        my $aa = @{$hash{$key}}[0];
                                        my $chain = @{$hash{$key}}[1];
                                        my $pos = @{$hash{$key}}[2];
                                        $relPos = $chain . $pos;
                                        $res{$relPos} = $aa;
                                        $positions{$pos} = $chain;
                                }
                                $residues{$name} = \%res;
                        }
                }
        }

        print OUT "PDB";
        my @values = ();
        foreach my $key (keys %positions)
        {
                my $val = $positions{$key} . $key;
                push @values, $val;
                print OUT "\t$val";
        }
        print OUT "\n";

        my $number = scalar @values;

        foreach my $key (keys %residues)
        {
                print OUT "$key";
                for(my $i=0;$i<$number;$i++)
                {
                        foreach my $key2 (keys %{$residues{$key}})
                        {
                                if($values[$i] eq $key2)
                                {
                                        print OUT "\t$residues{$key}{$key2}";
                                }
                        }
                }
                print OUT "\n";
        }
        close(OUT);
}


exit;
