#!usr/bin/perl
#********************************************************************
# MSc Bioinformatics with Systems Biology 2015/16
# 
# Analysis of the VH/VL packing
#
# Author:       Irene del Hierro García
# Script name:  nominal.pl
# Version:      V1.1
# Date:         19/08/16
#
#********************************************************************
#
# Purpose
# =======
#        Predict the values from a nominal MultiLayer Perceptron model
#
#********************************************************************
#
# Assumptions
# ===========
#       1. The program requires two input files, one with the nodes of the model
#          and other with the values of the predicted element.
#       2. The hidden layers must be sigmoids and the output layer linear
#
#       Usage: perl nominal.pl < model_MLP.txt
#
#********************************************************************   
# Strategy
# ========
#       1. Open and read each line of the input file.
#       2. Get the values of the nodes
#       3. Compute their weights
#       4. Calculate the linear output node
#       5. Give the final result
#
#********************************************************************
#
# Revision history
# ================ 
#       V1.0:   Original
#       V1.1:   Comments
#
#********************************************************************

use strict;


my $input = $ARGV[0];   # The file with the inputs of the predicted value
my $model = $ARGV[1];   # The model

# Open the input file and get the values 
open(IN,$input) or die "Unable to open the file\n";

my %values=();

while (my $line =<IN>)
{
        chomp $line;
        if($line =~ /(.*)\s{5}(.*)/)
        {
                $values{$1} = $2;
        }
}

close(IN);

# Open the model
open(MOD,$model) or die "Unable to open the model\n";

my %sigmoids =();

# Get the sigmoid nodes and store their weights and threshold in a hash of hashes
while(my $line2 = <MOD>)
{
        chomp $line2;

        if($line2 =~ /Sigmoid\sNode\s(\d*)/)
        {
               
                my $number = $1;
                $line2 =<MOD>; $line2 = <MOD>;
                my %attrib=();
                
                if($line2 =~ /\s{4}Threshold\s{4}(-?.*)/)
                {
                        $attrib{"Thres"} = $1;
                        

                        for(my $i=0;$i<72;$i++)
                        {
                                $line2 = <MOD>;
                                if($line2 =~ /\s{4}Attrib\s(.*=?.*)\s{4}(-?.*)/)
                                {
                                        $attrib{$1} = $2;
                                        
                                }
                        }
                }
                $sigmoids{$number} = \%attrib;
        }
}


close(MOD);

my %nodes = ();
# Calculate the weights of each node and sotre them in a different hash
foreach my $key (keys %sigmoids)
{
        my $sum = 0;

        foreach my $key2 (keys %{$sigmoids{$key}})
        {
                
                foreach my $val (keys %values)
                {
                        my $pos=""; my $res="";
                        if($key2 =~ /(.*)=(.*)/){
                                $pos = $1; $res = $2;
                        }
                        if($val eq $pos and $values{$val} eq $res)
                        {
                                my $x = 1 * $sigmoids{$key}{$key2};
                                $sum += $x;
                        }
                }
        }
        # Sum the total, add the threshold and calculate the activation function
        my $total = $sum + $sigmoids{$key}{"Thres"};       
        my $sig = 1 / (1 + exp(-$total));
        $nodes{$key} = $sig;
}

open(M,$model) or die "Unable to open the model\n";

my %linear=();
my $c=scalar keys %sigmoids;
# Get the linear output node and store its weight and threshold in a hash
while(my $lin =<M>)
{
        chomp $lin;
        if($lin =~ /Linear\sNode\s0/)
        {
                
                $lin =<M>; chomp $lin;$lin = <M>;chomp $lin;
                
                if($lin =~ /\s{4}Threshold\s{4}(-?.*)/)
                {
                        $linear{"Thres"} = $1;
                      

                        for(my $i=0;$i<$c;$i++)
                        {
                                $lin = <M>;
                                chomp $lin;
                                if($lin =~ /^\s{4}Node\s(.*)\s{4}(-?.*)/)
                                {
                                        $linear{$1} = $2;
                                        
                                }
                        }
                }
        }
}
close(M);



# Calculate the total weight for the linear node
my $w=0; my $thres = 0;

foreach my $key (keys %linear)
{
        foreach my $key2 (keys %nodes)
        {
                if($key =~ /$key2/)
                {
                       
                        my $weight = $nodes{$key2} * $linear{$key};
                        
                        $w += $weight;
                }
        }
}

$w += $linear{"Thres"};
my $sigma = 1 / $w;
print "The prediction of sd is $sigma\n";


