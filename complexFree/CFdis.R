######################
# Script name: CFdis.R
# Author: Irene del Hierro
# 
# Requirements: 	1 txt input files in tabular format.
#		    	SET\tSD
#			The SET column should be classified into Mixed, Free or Complex.
#
# Description: 	Represent a normal and a kernel density plot distribution
#			for each class in the SET column.
#
###################### 
rm(list=ls(all=TRUE))

# INPUT FILE
input<-read.table("CFdistro.txt",header=TRUE,sep="\t")
head(input)

free<-input$SD[input$SET=="Free"]
complex<-input$SD[input$SET=="Complex"]
mixed<-input$SD[input$SET=="Mixed"]

# Means
f.mean<-mean(input$SD[input$SET=="Free"])
c.mean<-mean(input$SD[input$SET=="Complex"])
m.mean<-mean(input$SD[input$SET=="Mixed"])

# Sd
f.sd<-sd(input$SD[input$SET=="Free"])
c.sd<-sd(input$SD[input$SET=="Complex"])
m.sd<-sd(input$SD[input$SET=="Mixed"])

# Get the minimums and maximums of the three classes
minF<-min(input$SD[input$SET=="Free"])
maxF<-max(input$SD[input$SET=="Free"])
minC<-min(input$SD[input$SET=="Complex"])
maxC<-max(input$SD[input$SET=="Complex"])
minFC<-min(input$SD[input$SET=="Mixed"])
maxFC<-max(input$SD[input$SET=="Mixed"])

# Print the distributions (assuming a normal distribution for all of them)
x<-seq(from=0,to=60,length.out=500)

par(mfrow=c(1,2))
plot(x,dnorm(x,mean=f.mean,sd=f.sd),ty="l",ylim=c(0,1),col= "green",ylab="Distribution",main="SD distribution of complex, free and mixed groups")
lines(x,dnorm(x,mean=c.mean,sd=c.sd),ty="l",ylim=c(0,1),col="red")
lines(x,dnorm(x,mean=m.mean,sd=m.sd),ty="l",ylim=c(0,1),col="blue")
legend(x=30.5,y=0.9,c("free","complex","mixed"),lty=c(1,1,1),lwd=c(2.5,2.5,2.5),col=c("green","red","blue"))

# Kernel density plot
source("https://bioconductor.org/biocLite.R")
#biocLite("sm") uncomment to install
library(sm)

type.f<-factor(input$SET) 
type.f

sm.density.compare(input$SD,type.f,xlab="SD distrbution",xlim=c(0,60))
title("SD distribution of the redundant groups")

# add legend via mouse click
colfill<-c(2:(2+length(levels(type.f))))
legend(locator(1), levels(type.f), fill=colfill) 

# Histograms
par(mfrow=c(1,3))
h<-hist(free, breaks=120, col="red", xlab="SD",
   main="Free",ylim=c(0,12),xlim=c(0,maxF))


h<-hist(complex, breaks=120, col="red", xlab="SD",
   main="Complex",ylim=c(0,12),xlim=c(0,maxC))


h<-hist(mixed, breaks=120, col="red", xlab="SD",
   main="Mixed",ylim=c(0,12),xlim=c(0,maxFC))




