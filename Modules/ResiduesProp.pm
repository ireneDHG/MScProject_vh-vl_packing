package ResiduesProp;
#********************************************************************
# MSc Bioinformatics with Systems Biology 2015/16
# 
# Analysis of the VH/VL packing
#
# Author:       Irene del Hierro García
# Module name:  ResiduesProp.pm
#
# Description:  This module contains subroutines that are involved in the
#               calculation of the properties of the redundant groups
#
#********************************************************************

#package to calculate the properties of some protein residues
use strict;
use config;
# To have the paths
my $dir;
if(!defined($dir))
{
        $dir = $config::propdir;
}

#********************************************************************
# Purpose: Subrotuine to get the key determining residues of VH/VL interface
#
# Date: 12/08/16
# Version: V1.1
#
# Arguments:
#       string $_[0]: pdb file
#
# Requirements:
#       1. Precise 1 argument
#       2. $_[0] should be a pdb file
#
# Return:
#       hash of arrays reference with the residues of the interface
#
# Die if unable to open the file
#********************************************************************
sub GetAaInterfaceVHVLCode
{	
        my ($file)=@_;
        
        open(IN,"$dir/$file") or die "Unable to open the pdb file\n";
        
        # Set the position and chain of key determining residues
        my %interRes = (
                L => ["38","40","41","44","46","87"],
                H => ["33","42","45","60","62","91","105"],
        );
        
        my $count = 0;
        
        my %aaCode = ();

        # Get the position and check if it corresponds with the ones in the hash (just take the CA atoms)	
        while (my $line = <IN>)
        {
                if($line =~ /^ATOM\s+\d+\s+(\w+)\s+(\w{3})\s(\w)\s+(\w+)\s+.+\s+.+\s+.+\s+.+\s.+\s+\w\s{2}$/)
                {
                        foreach my $key (keys %interRes)
                        {
                                foreach my $res (@{$interRes{$key}})
                                {
                                        # If they are the same, add to the hash
                                        if($1 eq "CA" and $key eq $3 and $res eq $4)
                                        {
                                                $count++;	
                                                push @{$aaCode{$count}}, uc($2),$3,$4;
                                        }
                                }
                        }
                }
        }

        close(IN);
        return (\%aaCode);
}

#********************************************************************
# Purpose: Subrotuine to get the relevant key determining residues of VH/VL interface
#
# Date: 12/08/16
# Version: V1.1
#
# Arguments:
#       string $_[0]: pdb file
#
# Requirements:
#       1. Precise 1 argument
#       2. $_[0] should be a pdb file
#
# Return:
#       hash of arrays reference with the residues of the interface
#
# Die if unable to open the file
#********************************************************************
sub GetAaRelevantVHVLCode
{	
        my ($file)=@_;
        
        open(IN,"$dir/$file") or die "Unable to open the pdb file\n";
        
        # Set the position and chain of key determining residues
        my %interRes = (
                L => ["41","87","38"],
                H => ["105","45","60"],
        );
        
        my $count = 0;
        
        my %aaCode = ();
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
                                                $count++;	
                                                push @{$aaCode{$count}}, uc($2),$3,$4;
                                        }
                                }
                        }
                }
        }

        close(IN);
        return (\%aaCode);
}

#********************************************************************
# Purpose: Subrotuine to classify the residues by being hydrophobic, hydrophilic, amphipathic or glycines
#
# Date: 12/08/16
# Version: V1.1
#
# Arguments:
#       string $_[0]: reference to a hash of arrays
#       string $_[1]: a pdb file
#       string $_[2]: the name of the first output file
#       string $_[3]: the name of the second output file
#       string $_[4]: the standard deviation of the group
#
# Requirements:
#       1. Precise 5 arguments
#       2. $_[0] should be a hash reference with the sequence (the residues in the first position of the array,
#               the chain in the second and the position in the third)
#       3. $_[2] should be a name for a normal txt file
#       4. $_[3] should be a name for a R file
#
# Return:
#       print in two output files
#
#********************************************************************
sub SetResidueType 
{
        my($hashref,$file,$name,$rname,$sd)=@_;
        
        my %aaCount = %$hashref;
        
        # Set the types of the residues
        my %aaTypes = (
                "Hydrophobic" => ["ALA","VAL","PHE","PRO","LEU","ILE"],
                "Hydrophilic" => ["ARG","ASP","GLU","SER","THR","CYS","ASN","GLN","HIS"],
                "Amphipathic" => ["LYS","TYR","MET","TRP"],
                "Other" => ["GLY"],
       );
        
        my $phobic=0; my $philic=0; my $gly=0; my $amphi=0;
        my $total=0; my @sequenceL=(); my @sequenceH=();
        
        my %finalHash = ();
        
        # For each key of the hash, get the sequence, the position and the chain
        foreach my $key (keys %aaCount)
        {
                my $seq="";
                my $aa = @{$aaCount{$key}}[0];
                $total++;
                my $chain = @{$aaCount{$key}}[1];
                my $pos = @{$aaCount{$key}}[2];
                # Join the chain and position
                $seq = $aa . $pos;
                
                # If the chain is light (L) or heavy (H) add to their respective arrays
                if($chain eq "L")
                {
                        push @sequenceL, $seq;
                }
                elsif($chain eq "H")
                {
                        push @sequenceH, $seq;
                }
                
                # For each type of residues, check the type of the amino acids introduced
                # and update counters
                foreach my $key2 (keys %aaTypes)
                {
                        foreach my $element (@{$aaTypes{$key2}})
                        {
                                if($key2 eq "Hydrophobic" and $aa eq $element)
                                {
                                        $phobic++;
                                }
                                elsif($key2 eq "Hydrophilic" and $aa eq $element)
                                {
                                        $philic++;
                                }
                                elsif($key2 eq "Amphipathic" and $aa eq $element)
                                {
                                        $amphi++;
                                }
                                elsif($key2 eq "Other" and $aa eq $element)
                                {
                                        $gly++;
                                }
                        }
                }
        }
        
        my $percPhobic=0; my $percPhilic=0; my $percGly=0; my $percAmphi=0;
        
        # Calculate the percentage of each type of residues
        if($phobic == 0)
        {
                $percPhobic = 0;
        }
        else
        {
                $percPhobic = ($phobic/$total)*100;
        }

        if($philic == 0)
        {
                $percPhilic = 0;
        }
        else
        {
                $percPhilic = ($philic/$total)*100;
        }
        
        if($gly == 0)
        {
                $percGly = 0;
        }
        else
        {
                $percGly = ($gly/$total)*100;
        }

         if($amphi == 0)
        {
                $percAmphi = 0;
        }
        else
        {
                $percAmphi = ($amphi/$total)*100;
        }
        
        # Print the results in the output files.
        print $name "$file: Hydrophobic $percPhobic% Hydrophilic $percPhilic% Amphipatic $percAmphi% Glycines $percGly%\n";
        print $name "L-CHAIN: @sequenceL\n"; 
        print $name "H-CHAIN: @sequenceH\n";
        print $name "Group SD: $sd\n\n";   

        print $rname "$file\t$sd\t$percPhobic\t$percPhilic\t$percAmphi\t$percGly\n"; 
}

#********************************************************************
# Purpose: Subrotuine to calculate the average hydrophobicity of a given sequence
#
# Date: 12/08/16
# Version: V1.1
#
# Arguments:
#       string $_[0]: reference to a hash of arrays
#
# Requirements:
#       1. Precise 1 argument
#       2. $_[0] should be a hash reference with the sequence (the residues in the first position of the array)
#
# Return:
#       hydrophobicity index
#
#********************************************************************
sub CalculateHydropAVG
{
        my($hashRef) =@_;
        
        # Kyte-Doolittle scale
        my %hydroIndex = (
                "ILE" => "4.5",
                "VAL" => "4.2",
                "LEU" => "3.8",
                "PHE" => "2.8",
                "CYS" => "2.5",
                "MET" => "1.9",
                "ALA" => "1.8",
                "GLY" => "-0.4",
                "THR" => "-0.7",
                "SER" => "-0.8",
                "TRP" => "-0.9",
                "TYR" => "-1.3",
                "PRO" => "-1.6",
                "HIS" => "-3.2",
                "GLU" => "-3.5",
                "GLN" => "-3.5",
                "ASP" => "-3.5",
                "ASN" => "-3.5",
                "LYS" => "-3.9",
                "ARG" => "-4.5",
        );

        my %aaHash = %$hashRef; my @indexes;
        
        # For each key of the hash, get the residue, and push the index in the array
        foreach my $key (keys %aaHash)
        {
                my $elem = @{$aaHash{$key}}[0];

                foreach my $aa (keys %hydroIndex)
                {
                        if($aa eq $elem)
                        {
                                push @indexes, $hydroIndex{$aa};
                        }
                }       
        }
        
        my $sum = 0;
        # For each index, sum all values
        foreach my $index (@indexes)
        {
        	$sum += $index;
        }
        
        my $nElem = scalar(@indexes);
        # Calculate the average
        my $hydrophobicity = $sum/$nElem;
        	
        return ($hydrophobicity);
}

#********************************************************************
# Purpose: Subrotuine to calculate the molecular weight of a given sequence
#
# Date: 12/08/16
# Version: V1.1
#
# Arguments:
#       string $_[0]: reference to a hash of arrays
#
# Requirements:
#       1. Precise 1 argument
#       2. $_[0] should be a hash reference with the sequence (the residues in the first position of the array)
#
# Return:
#       molecular weight
#
#********************************************************************
sub CalculateMW
{
        my ($refhash)=@_;

        # Hash with the values for each residue (obtained from Expasy)
        my %hashMW = (
                "ILE" => "131.00",
                "VAL" => "117.00",
                "LEU" => "131.00",
                "PHE" => "165.00",
               	"CYS" => "121.00",
                "MET" => "149.00",
                "ALA" => "89.00",
                "GLY" => "75.00",
                "THR" => "119.00",
                "SER" => "105.00",
                "TRP" => "204.00",
                "TYR" => "181.00",
                "PRO" => "115.00",
                "HIS" => "155.00",
                "GLU" => "147.00",
                "GLN" => "146.00",
                "ASP" => "133.00",
                "ASN" => "132.00",
                "LYS" => "146.00",
                "ARG" => "174.00",
        );
        
        my %aaHash = %$refhash; my @molws;
        
        # For each key of the hash, get the residue, and push the index in the array
        foreach my $key (keys %aaHash)
        {
                my $elem = @{$aaHash{$key}}[0];
                
                foreach my $aa (keys %hashMW)
                {
                        if($aa eq $elem)
                        {
                                push @molws, $hashMW{$aa};
                        }
                }
        }
        
        my $sum = 0;
        # For each index, sum all values
        foreach my $mol (@molws)
        {
                $sum += $mol;
        }
        
        my $nElem = scalar(@molws);
        my $mw = $sum;
        	
        return ($mw);
}

#********************************************************************
# Purpose: Subrotuine to calculate the average isoelectric point of a given sequence
#
# Date: 12/08/16
# Version: V1.1
#
# Arguments:
#       string $_[0]: reference to a hash of arrays
#
# Requirements:
#       1. Precise 1 argument
#       2. $_[0] should be a hash reference with the sequence (the residues in the first position of the array)
#
# Return:
#       average isoelectric point
#
#********************************************************************
sub CalculatePIAVG
{
        my ($refhash) = @_;
        
        # From www.chem.ucalgary.ca and geneinfinity.org
        my %hashPI = (
                "ILE" => "6.02",
                "VAL" => "5.96",
                "LEU" => "5.98",
                "PHE" => "5.48",
                "CYS" => "5.07",
                "MET" => "5.74",
                "ALA" => "6.00",
                "GLY" => "5.97",
                "THR" => "5.60",
                "SER" => "5.68",
                "TRP" => "5.89",
                "TYR" => "5.66",
                "PRO" => "6.30",
                "HIS" => "7.59",
                "GLU" => "3.22",
                "GLN" => "5.65",
                "ASP" => "2.77",
                "ASN" => "5.41",
                "LYS" => "9.74",
                "ARG" => "10.76",
        );
        
        my %aaHash = %$refhash; my @isoelec;
        
         # For each key of the hash, get the residue, and push the index in the array
        foreach my $key (keys %aaHash)
        {
                my $elem = @{$aaHash{$key}}[0];
        		
                foreach my $aa (keys %hashPI)
                {
                        if($aa eq $elem)
                        {
                                push @isoelec, $hashPI{$aa};
                        }
                }
        }

        my $sum = 0;
        # Sum all elements and calculate the average
        foreach my $iso (@isoelec)
        {
                $sum += $iso;
        }
        	
        my $nElem = scalar(@isoelec);
        my $pi = $sum/$nElem;
	
        return ($pi);
}

#********************************************************************
# Purpose: Subrotuine to calculate the isoelectric point of a given protein
#
# Obtained from http://www-nmr.cabm.rutgers.edu/bioinformatics/ZebaView/help.html
#
#
# Return:
#       isoelectric point
#
#********************************************************************

sub CalculatePIProtein
{
        my $protein= "MKCLLLALALTCGAQALIVTQTMKGLDIQKVAGTWYSLAMAASDISLLDAQSAPLRVYVEELKPT
        +PEGDLEILLQKWENGECAQKKIIAEKTKIPAVFKIDALNENKVLVLDTDYKKYLLFCMENSAEPE
        +QSLACQCLVRTPEVDDEALEKFDKALKALPMHIRLSFNPTQLEEQCHI";

        my %residue;
        my $allowed="ACDEFGHIKLMNPQRSTVWY";
        my @acids=('C','D','E','Y');       # Acidic for pI
        my @bases=('H','K','R');           # Basic for pI

        # pKa/pKbs (< & > are default NT & CT)
        my %pka = (
                'C' => "9.0",  
                'D' => "4.05",
                'E' => "4.45",
                'Y' => "10.0",
                'H' => "5.98", 
                'K' => "10.0",  
                'R' => "12.0",   
                '<' => "7.5", 
                '>' => "3.55",
        );

        # NT pKa        #pI should never be over 14
        my %nt = (
                'A' => "7.59", 
                'E' => "7.70", 
                'M' => "7.00",
                'P' => "8.36", 
                'S' => "6.93", 
                'T' => "6.82", 
                'V' => "7.44",
        );

        # CT pKa
        my %ct = (
                'D' => "4.55", 
                'E' => "4.75",
        );

        # A280s
        my %A280 = (
                'W' => "5500.", 
                'Y' => "1490.",
        );

        study($protein);                        # delete non-aa (for some reason
        $protein =~ s/![$allowed]|[\s\n]//go;   # space and \n get spec treatment)
        my $no_aa = length($protein);

        my $first = substr($protein,0,1);  # get NT aa
        my $last = substr($protein, -1,1); # get CT aa 

        my %base;
        my %acid;
        $base{'<'}=1;                   # One amino terminus
        $acid{'>'}=1;                   # One carboxy terminus
 

        if($nt{$first}) {               # NT different for some aa, 7.50 default
                $pka{'<'} = $nt{$first};
        }

        if($ct{$last}) {                # CT different for some aa, 3.55 default
                $pka{'>'} = $ct{$last};
        }
        # amino acid composition
        foreach my $aa (split(//, $allowed)) {
                $residue{$aa} = ($protein =~ s/$aa//g);
        }

        foreach my $aa (@acids) {
                if($residue{$aa}) {
                        $acid{$aa} = $residue{$aa}; # collect acids
                }
        }
        foreach my $aa (@bases) {
                if($residue{$aa}) {
                        $base{$aa} = $residue{$aa}; # collect bases
                }
        } 
                                        # binary search for pI
        my $hi=14;                         # Theoretical max
        my $lo=0;                          # Theoretical min
        my $pI=7;
        my $old_pI=0;
        my $iterations=0;

        while(abs($pI-$old_pI)>0.001) { # Two correct decimals
                if($iterations++ > 15) {      # 14/0.001 > 2^14, so this shouldn't happen
                        last;
                }
                my $result=0;
                foreach my $aa (keys(%acid)) {   # Left (acid) side
                        $result += $acid{$aa}/(1+10**($pka{$aa}-$pI));

                }
                foreach my $aa (keys(%base)) {   # Right (base) side
                        $result -= $base{$aa}/(1+10**($pI-$pka{$aa}));

                }
                $old_pI = $pI;
                if($result > 0){
                        $pI=($pI+$lo)/2;            # Go lower since charge is neg
                        $hi=$old_pI;

                } else {
                        $pI=($hi+$pI)/2;            # Go higher; charge is pos
                        $lo=$old_pI;
                }
        }   
        return($pI);   
}


1;
