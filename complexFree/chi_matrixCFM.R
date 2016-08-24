######################
# Script name: chi_matrixCFM.R
# Author: Irene del Hierro
# 
# Requirements: 	1 txt input file in matrix format
#		    	i.e. 	SD	F	C	F/C
#				Low	216	251	110
#				High	8	2	5
#
# Description: 	Perform chi-squared tests to determine the differences amongst
#			the F, C and F/C redundant groups.
#
###################### 
rm(list=ls(all=TRUE))

#INPUT FILE MATRIX
mat<-read.table("div6matCF.txt",header=TRUE,sep="\t")
mat

# F/C:L vs F/C:H differences:
FC_L<-mat[1,4]
FC_L
FC_H<-mat[2,4]
FC_H
FCtotal<-FC_L+FC_H
FCtotal
expFC<-FCtotal/2;
expFC

obs<-c(FC_L,FC_H)
exp<-c(expFC,expFC)
exp

table<-matrix(data=c(obs,exp),nrow=2,ncol=2,byrow=TRUE)
colnames(table)<-c("Low","High")
rownames(table)<-c("O","E")
table
x.sq<-sum((obs-exp)^2/exp)
x.sq
qchisq(0.95,1)

#this is just to check...
prob<-c(0.5,0.5)
chisq.test(obs,p=prob)

#Check C:L vs C:H
C_L<-mat[1,3]
C_L
C_H<-mat[2,3]
C_H
Ctotal<-C_L+C_H
Ctotal
expC<-Ctotal/2;
expC

obs<-c(C_L,C_H)
exp<-c(expC,expC)

table<-matrix(data=c(obs,exp),nrow=2,ncol=2,byrow=TRUE)
colnames(table)<-c("Low","High")
rownames(table)<-c("O","E")
table
x.sq<-sum((obs-exp)^2/exp)
x.sq
qchisq(0.95,1)

#this is just to check...
prob<-c(0.5,0.5)
chisq.test(obs,p=prob)

#Check F:L vs F:H
F_L<-mat[1,2]
F_L
F_H<-mat[2,2]
F_H
Ftotal<-F_L+F_H
Ftotal
expF<-Ftotal/2;
expF

obs<-c(F_L,F_H)
exp<-c(expF,expF)

table<-matrix(data=c(obs,exp),nrow=2,ncol=2,byrow=TRUE)
colnames(table)<-c("Low","High")
rownames(table)<-c("O","E")
table
x.sq<-sum((obs-exp)^2/exp)
x.sq
qchisq(0.95,1)

#this is just to check...
prob<-c(0.5,0.5)
chisq.test(obs,p=prob)

# Compare complex and free only groups!
F_L
F_H
C_L
C_H
obs<-c(F_L,C_L,F_H,C_H)
totalL<-F_L+C_L
totalH<-F_H+C_H
tot<-totalL + totalH
exp<-c(tot/4,tot/4,tot/4,tot/4)
exp
x.sq<-sum((obs-exp)^2/exp)
x.sq
qchisq(0.95,1)
prob<-c(1/4,1/4,1/4,1/4)
chisq.test(obs,p=prob)

# Full table test
obsL<-c(F_L,C_L,FC_L)
obsH<-c(F_H,C_H,FC_H)
obsL
obsH

totalL<-sum(obsL)
totalH<-sum(obsH)
tot<-totalL + totalH

exp<-c(tot/6,tot/6,tot/6,tot/6,tot/6,tot/6)
exp

obs<-c(obsL,obsH)
x.sq<-sum((obs-exp)^2/exp)
x.sq
qchisq(0.95,2)

prob<-c(1/6,1/6,1/6,1/6,1/6,1/6)
chisq.test(obs,p=prob)