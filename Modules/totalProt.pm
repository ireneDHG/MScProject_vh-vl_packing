package totalProt;
#********************************************************************
# MSc Bioinformatics with Systems Biology 2015/16
# 
# Analysis of the VH/VL packing
#
# Author:       Irene del Hierro García
# Module name:  totalProt.pm
#
# Description:  This module contains subroutines that are involved in extracting
#               the sequence and counting amino acids
#
#********************************************************************

use strict;
use config;

# Set the path
my $dir;
if(!defined($dir))
{
        $dir = $config::propdir;
}

#********************************************************************
# Purpose: Get the sequence from a protein
#
# Date: 12/08/16
# Version: V1.1
#
# Arguments:
#       string $_[0]: a PDB file
#
# Requirements:
#       1. Precise 1 argument
#       2. $_[0] should be a PDB file
#
# Return:
#       An array reference with the sequence
#
#********************************************************************
sub GetSequenceFromPDB
{
        my ($file)=@_;
        
        open(IN,"$dir/$file") or die "Unable to open the pdb file $!\n";
        
        my @sequence=();

        while(my $line = <IN>)
        {
                # Get the sequence just using the CA atoms
                if($line =~ /^ATOM\s+\d+\s+(\w+)\s+(\w{3})\s(\w)\s+(\w+)\s+.+\s+.+\s+.+\s+.+\s.+\s+\w\s{2}$/)
                {
                        if($1 eq "CA")
                        {
                                push @sequence, $2;
                        }
                }
        }

        close(IN);
        return (\@sequence);
}

#********************************************************************
# Purpose: Count the number of each amino acid in a given sequence
#
# Date: 12/08/16
# Version: V1.1
#
# Arguments:
#       string $_[0]: an array reference
#
# Requirements:
#       1. Precise 1 argument
#       2. $_[0] should be an array reference with the amino acid sequence
#
# Return:
#       A hash reference with the counts for each amino acid
#
#********************************************************************
sub CountAaSequence
{
        my ($reference)=@_;
        
        my @sequence=@$reference;
                
        my %counts=();
        # All the amino acids
        my @totalAa = ("ILE","VAL","LEU","PHE","CYS","MET","ALA","GLY",
                       "THR","SER","TRP","TYR","PRO","HIS","GLU","GLN",
                       "ASP","ASN","LYS","ARG");
         
        # Add to a hash if a new amino acid is found in the sequence, else,
        # update the counter     
        foreach my $aa (@sequence)
        {
                $counts{$aa}++;
        }
        
        # Add those amino acids not present in the sequence to the hash
        foreach my $elem (@totalAa)
        {
                unless(exists $counts{$elem})
                {
                        $counts{$elem}=0;
                }
        }
        
        return(\%counts);
}


        
        
              
        

1;

