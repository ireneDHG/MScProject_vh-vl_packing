#!usr/bin/perl -w
#********************************************************************
# MSc Bioinformatics with Systems Biology 2015/16
# 
# Analysis of the VH/VL packing
#
# Author:       Irene del Hierro García
# Script name:  add_packAngle.pl
# Version:      V1.1
# Date:         12/08/16
#
#********************************************************************
#
# Purpose
# ========
#         Add the packing angles calculated by a C program (previously installed) to a series of antibodies 
#         PDB files.
#         The antibodies have been numerated by the Chothia numbering scheme.
#         Each line of the list of antibodies (text file) contains -if available- multiple
#         structures of the same antibody.
#         The program is then, going to analyze that list, take the name of the antibody, 
#         extract it from the folder specified (LH_Combined_Chothia) and give it to the C program
#         to compute the angle.
#
#********************************************************************
#
# Assumptions
# ===========
#       1. The program requires an input txt file with the list of antbodies grouped by sequence and
#          the presence of a folder containing the pdb files of that list (the path is defined in config.pm).
#       2. The antibodies of the groups in the txt file should be separated by commas
#       3. The location of the modules with the subroutines used in this script should be located in the appropriate path
#          (change the location of the Modules folder to correctly used this program)
#       4. Unless a screen printing is required, add an output file to obtain the results.
#       5. Installation of the C program is required as well as defining the environment variables required.
#
#       Usage: perl add_packAngle.pl < Redundant_LH_Combined_Chothia.txt > RedundantChothia_Angles.txt
#  
#********************************************************************   
# Strategy
# ========
#       1. Open the folder with the PDB files and stored the files in an array.
#       2. Call the packingAngle module to calculate the angle (using a C program) for each PDB file.
#       3. The hash retrieved is introduced again in another subroutine of the packingAngle module to
#          add the angles below their corresponding antibody names in the input file.
#
#********************************************************************            
#
# Revision history
# ================
#       V1.0:   Original
#       V1.1:   Comments and polish
#
#********************************************************************

# Location of the modules
use lib "/home/irene/Documents/MScProject/Modules"; 

use strict;
use config;             # Configuration of the path of the input folder
use packingAngle;       # Subroutines for angle analysis

# Get the directory of the input folder
my $dir;
if(!defined($dir))
{
        $dir = $config::antibodyDir;
}


# Get the PDB files from the directory
opendir (DIR, $dir) or die "Unable to open dir $dir\n";
my @pdbFiles = grep(/.*\.pdb$/,readdir(DIR));
closedir(DIR);

# Create the hash to store the pdbs and the angles
my $pdbRefHash = packingAngle::CreateHashAnglePdb(\@pdbFiles);	

# Call to a subroutine to add the angles with the pdb codes in the output file	
packingAngle::AddAngleToFile($pdbRefHash);
	
exit;
