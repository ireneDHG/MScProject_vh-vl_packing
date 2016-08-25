#!usr/bin/perl -w
#********************************************************************
# MSc Bioinformatics with Systems Biology 2015/16
# 
# Analysis of the VH/VL packing
#
# Author:       Irene del Hierro García
# Script name:  numeric.pl
# Version:      V1.1
# Date:         19/08/16
#
#********************************************************************
#
# Purpose
# =======
#        Predict the values from a numeric MultiLayer Perceptron model
#
#********************************************************************
#
# Assumptions
# ===========
#       1. The program requires two input files, one with the nodes of the model
#          and other with the values of the predicted element.
#       2. The hidden layers must be sigmoids and the output layer linear
#
#       Usage: perl numeric.pl < model_MLP.txt
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

# The first input has the values of the predicted element, the second the model
my $input = $ARGV[0];
my $model = $ARGV[1];

open(IN,$input) or die "Unable to open the file\n";

my %values=();

# Get the values of the predicted thing
while (my $line =<IN>)
{
        chomp $line;
        if($line =~ /(.*)\s{7}(.*)/)
        {
                $values{$1} = $2;
        }
}

# Normalize values:
my %normVal =();
my $max=0;
$_ > $max and $max = $_ for values %values;
my $min=$max;
$_ < $min and $min = $_ for values %values;

foreach my $key (keys %values)
{
        my $norm = ($values{$key} - $min)/($max-$min);
        
        $normVal{$key} = $norm;
}


close(IN);

# Open the model 
open(MOD,$model) or die "Unable to open the model\n";

my %sigmoids =();
my $c=scalar keys %values;

# Get every sigmoid node or hidden layer and store the information in a hash of hashes
while(my $line2 = <MOD>)
{
        chomp $line2;

        if($line2 =~ /Sigmoid\sNode\s(\d*)/)
        {
                my $number = $1;
                $line2 =<MOD>; $line2 = <MOD>;
                my %attrib=();
                
                # Get the threshold of the node
                if($line2 =~ /\s{4}Threshold\s{4}(-?.*)/)
                {
                        $attrib{"Thres"} = $1;

                        # Get all the weights of the node
                        for(my $i=0;$i<$c;$i++)
                        {
                                $line2 = <MOD>;
                                if($line2 =~ /\s{4}Attrib\s(.*)\s{4}(-?.*)/)
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
# For each sigmoid node, compute their weights taking into account the input of the predicted values
foreach my $key (keys %sigmoids)
{
        my $sum = 0;

        foreach my $key2 (keys %{$sigmoids{$key}})
        {
                
                foreach my $val (keys %normVal)
                {
                        if($val eq $key2)
                        {
                                my $x = $normVal{$val} * $sigmoids{$key}{$key2};
                                $sum += $x;
                        }
                }
        }
        # Calculate the total and the activation function
        my $total = $sum + $sigmoids{$key}{"Thres"};       
        my $sig = 1 / (1 + exp(-$total));
        $nodes{$key} = $sig;
}


open(M,$model) or die "Unable to open the model\n";

my %linear=();
my $s=scalar keys %sigmoids;
# Get the linear output node and store the results in a different hash
while(my $lin =<M>)
{
        chomp $lin;
        if($lin =~ /Linear\sNode\s0/)
        {
                
                $lin =<M>; chomp $lin;$lin = <M>;chomp $lin;
                
                # Add the threshold to the hash
                if($lin =~ /\s{4}Threshold\s{4}(-?.*)/)
                {
                        $linear{"Thres"} = $1;
                       
                        # Get the weights for every sigmoid node
                        for(my $i=0;$i<$s;$i++)
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

# Calculate the weight of the output node
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

# Add threshold and then calculate the inverse of the activation function
$w += $linear{"Thres"};
my $sigma = 1 / $w;

# Denormalize
my $a = ($sigma * ($max - $min))+ $min;
print "The prediction is: $a\n";




