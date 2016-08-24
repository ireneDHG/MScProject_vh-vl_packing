######################
# Script name: graphSD.R
# 
# Requirements: 	4 txt input files in tabular format.
#		    	Group number \t sd value
#
# Description: 	Determine the distribution and represent histograms and kernel 
#			density plots for a general distribution and three group size
#			thresholds distributions
#
###################### 

# General distribution
inputSD<-read.table("sd_values.txt", header=TRUE, sep="\t")
head(inputSD)


# Test the type of the distribution
ks.test(inputSD$SD, "pnorm", mean(inputSD$SD), sd(inputSD$SD))
shapiro.test(inputSD$SD)
ks.test(inputSD$SD,"ppois",mean(inputSD$SD),sd(inputSD$SD))
ks.test(inputSD$SD,"pbeta",mean(inputSD$SD),sd(inputSD$SD))

par(mfrow=c(1,2))
# Create the whole SD histogram:
h<-hist(inputSD$SD,breaks=120,col="red",
xlab="Standard deviation",main="Histogram")

# Kernel Density Plot
d <- density(inputSD$SD)
plot(d, main="Kernel Density Plot", xlab= "Standard deviation")

# Calculate mean and sd
mean(inputSD$SD)
sd(inputSD$SD)

par(mfrow=c(2,3))

# For group size threshold 3:
in3<-read.table("sd_AllvaluesThres_3.txt", header=TRUE, sep="\t")
head(in3)
h<-hist(in3$SD,breaks=120,col="red",xlab="SD",
main="Threshold 3") 

# For group thereshold 4:
in4<-read.table("sd_AllvaluesThres_5.txt", header=TRUE, sep="\t")
head(in4)
h<-hist(in4$SD,breaks=120,col="red",xlab="SD",
main="Threshold 4")

# For group thereshold 5:
in5<-read.table("sd_AllvaluesThres_5.txt", header=TRUE, sep="\t")
head(in5)
h<-hist(in5$SD,breaks=120,col="red",xlab="SD",
main="Threshold 5")

# Kernels
d <- density(in3$SD)
plot(d,main="Density",xlab="SD")
e <- density(in4$SD)
plot(e,main="Density",xlab="SD")
f <- density(in5$SD)
plot(f,main="Density",xlab="SD")