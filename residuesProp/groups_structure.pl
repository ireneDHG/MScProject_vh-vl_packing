#!usr/bin/perl -w
#********************************************************************
# MSc Bioinformatics with Systems Biology 2015/16
# 
# Analysis of the VH/VL packing
#
# Author:       Irene del Hierro García
# Script name:  statistics.pl
# Version:      V1.1
# Date:         19/08/16
#
#********************************************************************
#
# Purpose
# =======
#        Calculate the proportions of the amino acids in every position of a given number of positions.
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
#       Usage: perl groups_structure.pl < ../freeabSD.txt
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
                        if(4 < $1)
                        {
                                $highSD{$names[0]}=$1;
                        }
                        elsif(2 > $1)
                        {
                                $lowSD{$names[0]}=$1;
                        }                
                }
        }
}

my $lowset="free_Str_L2";
my $highset = "free_Str_H2";

open(my $LOW,">$lowset.txt") or die "Unable to open the file low $!\n";
open(my $HIGH,">$highset.txt") or die "Unable to open the file high $!\n";

my %low_residues = ();
my %positions = ();

# For each of the proteins in the low set, open its PDB file and get the residues of the interface
foreach my $name (keys %lowSD)
{
        foreach my $file (@pdbFiles)
        {
                my $pdb = packingAngle::GetPdbFileName($file);
                if($name eq $pdb)
                {

                        my $sd=$lowSD{$name};
                        my $hashref = ResiduesProp::GetAaInterfaceVHVLCode($file);
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
                        $low_residues{$name} = \%res;
                }
        }
}

print $LOW "PDB";
my @values = ();
foreach my $key (keys %positions)
{
        my $val = $positions{$key} . $key;
        push @values, $val;
        print $LOW "\t$val";
}
print $LOW "\n";

my $number = scalar @values;

foreach my $key (keys %low_residues)
{
        print $LOW "$key";
        for(my $i=0;$i<$number;$i++)
        {
                foreach my $key2 (keys %{$low_residues{$key}})
                {
                        if($values[$i] eq $key2)
                        {
                                print $LOW "\t$low_residues{$key}{$key2}";
                        }
                }
        }
        print $LOW "\n";
}
                       
my %high_residues = ();

# Do the same but with the high sd set

foreach my $name (keys %highSD)
{
        foreach my $file (@pdbFiles)
        {
                my $pdb = packingAngle::GetPdbFileName($file);
                if($name eq $pdb)
                {
                        #print "$pdb\n";

                        my $sd=$highSD{$name};
                        my $hashref = ResiduesProp::GetAaInterfaceVHVLCode($file);
                        my %hash = %$hashref;
                        my %res = ();

                        foreach my $key (keys %hash)
                        {
                                my $relPos = "";
                                my $aa = @{$hash{$key}}[0];
                                my $chain = @{$hash{$key}}[1];
                                my $pos = @{$hash{$key}}[2];
                                $relPos = $chain . $pos;
                                $res{$relPos} = $aa;
                                #print "$relPos $aa\n";
                                #$high_rel_positions{$pos} = $chain;
                        }
                        $high_residues{$name} = \%res;
                }
        }
}

print $HIGH "PDB";
foreach my $a (@values)
{
        print $HIGH "\t$a";
}

print $HIGH "\n";

foreach my $key (keys %high_residues)
{
        print $HIGH "$key";
        for(my $i=0;$i<$number;$i++)
        {
                foreach my $key2 (keys %{$high_residues{$key}})
                {
                        if($values[$i] eq $key2)
                        {
                                print $HIGH "\t$high_residues{$key}{$key2}";
                        }
                }
        }
        print $HIGH "\n";
}
                 

close($LOW);
close($HIGH);

exit;
