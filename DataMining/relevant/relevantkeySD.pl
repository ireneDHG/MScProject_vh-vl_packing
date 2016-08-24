#!usr/bin/perl -w
#********************************************************************
# MSc Bioinformatics with Systems Biology 2015/16
# 
# Analysis of the VH/VL packing
#
# Author:       Irene del Hierro García
# Script name:  relevantkeySD.pl
# Version:      V1.1
# Date:         19/08/16
#
#********************************************************************
#
# Purpose
# =======
#        Create a CSV file with the counts of the relevant key determing positions and a numeric output.
#
#********************************************************************
#
# Assumptions
# ===========
#       1. The program requires an input file with a list of antibodies grouped by sequence.
#       2. The antibodies of the group should be separated by commas in the first line and in the next line 
#          the standard deviation and mean should appear as: "SD: 0.2344455 MEAN: 2.48".
#
#       Usage: perl relevantkeySD.pl < ../RedundantChothia_MeanSD.txt
#
#********************************************************************   
# Strategy
# ========
#       1. Open and read each line of the input file.
#       2. Get the sd and the interface residues for the file chosen
#       3. Add the counts of the amino acids
#       4. Print the results in an output file
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

# Call this script as: perl keyres.pl < ../RedundantChothia_MeanSD.txt 

#why #!/usr/bin/perl -s

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
my $name = "relevantKeyresSD";
my $aa = "I,V,L,F,C,M,A,G,T,S,W,Y,P,H,E,Q,D,N,K,R";
open(OUT,">$name.csv") or die "Unable to open the output file $!\n";  
print OUT "Protein,SD,$aa\n";

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
                                # If the sd is not undefined, get the residues from the interface and add calculate the counts
                                unless($1 =~ /undef/)
                                {
                                        my $sd = $1;
                                        
                                        my $seq= &GetSeqInterfaceVHVLCode($file);
                                        my $ref=totalProt::CountAaSequence($seq);
                                        my %hash =%$ref;
                                        print OUT "$pdb,$sd,$hash{ILE},$hash{VAL},$hash{LEU},$hash{PHE},$hash{CYS},$hash{MET},$hash{ALA},$hash{GLY},$hash{THR},$hash{SER},$hash{TRP},$hash{TYR},$hash{PRO},$hash{HIS},$hash{GLU},$hash{GLN},$hash{ASP},$hash{ASN},$hash{LYS},$hash{ARG}\n";
                                }
                          }
                          
                }
        }
}
      
close(OUT);

#********************************************************************
# Purpose: Get the sequence from the relevant key determining residues
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
sub GetSeqInterfaceVHVLCode
{	
        my ($file)=@_;
        
        open(IN,"$dir/$file") or die "Unable to open the pdb file\n";
        
        # Set the position and chain of key determining residues
        my %interRes = (
                L => ["41","38","87"],
                H => ["105","45","60"],
        );        
        my @sequence=();
        # Get the position and check if it corresponds with the ones in the hash (just take the CA atoms)	
        while (my $line = <IN>)
        {
                if($line =~ /^ATOM\s+\d+\s+(\w+)\s+(\w{3})\s(\w)\s+(\w+)\s+.+\s+.+\s+.+\s+.+\s.+\s+\w\s{2}$/)
                {
                        foreach my $key (keys %interRes)
                        {
                                foreach my $res (@{$interRes{$key}})
                                {
                                        if($1 eq "CA" and $key eq $3 and $res eq $4)
                                        {	
                                                push @sequence, $2;
                                        }
                                }
                        }
                }
        }

        close(IN);
        return (\@sequence);
}

exit;



