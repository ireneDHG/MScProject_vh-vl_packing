######################
# Script name: script3.R
# Author: Irene del Hierro
# 
# Requirements: 	7 pairs txt input files in tabular format (one for the high set
#			and the other for the low set):
#		    	HIGHpdb\tSD\tHYDRO\tMW\tPI for properties (high set)
#			LOWpdb\tSD\tHYDRO\tMW\tPI (low set)
#			PDB\tSD\tHydrophobic\tHydrophilic\tAmphipathic\tGlycines for types
#
# Description: 	Make t-test to test the differences in the properties between
#			the high and low standard deviation sets for division set 3.
#
###################### 

rm(list=ls(all=TRUE))

## PROPERTIES

par(mfrow=c(2,2))

highSD<-read.table("freediv3propH.txt", header=TRUE, sep="\t")
head(highSD)

lowSD<-read.table("freediv3propL.txt", header=TRUE, sep="\t")
head(lowSD)

highHydro<-highSD$HYDRO
lowHydro<-lowSD$HYDRO

t.test(highHydro,lowHydro,var=T)
boxplot(highHydro,lowHydro,names=c("SD>2","SD<=2"),main="Hydrophobicity")

highMW<-highSD$MW
lowMW<-lowSD$MW

highPI<-highSD$PI
lowPI<-lowSD$PI

t.test(highMW,lowMW,var=T)
boxplot(highMW,lowMW,names=c("SD>2","SD<=2"),main="Molecular weight")

t.test(highPI,lowPI,var=T)
boxplot(highPI,lowPI,names=c("SD>2","SD<=2"),main="Isoelectric point")

# Get the p-values
p1<-t.test(highHydro,lowHydro,var=T)$p.value
p2<-t.test(highMW,lowMW,var=T)$p.value
p3<-t.test(highPI,lowPI,var=T)$p.value


####### TYPES

#Types amino acids
typeH<-read.table("freeRdiv3typeH.txt", header=TRUE, sep="\t")
head(typeH)

typeL<-read.table("freeRdiv3typeL.txt", header=TRUE, sep="\t")
head(typeL)

phobicL<-typeL$Hydrophobic
amphiL<-typeL$Amphipathic
glyL<-typeL$Glycines
philicL<-typeL$Hydrophilic

phobicH<-typeH$Hydrophobic
amphiH<-typeH$Amphipathic
glyH<-typeH$Glycines
philicH<-typeH$Hydrophilic

t.test(phobicH,phobicL,var=T)
t.test(philicH,philicL,var=T)
t.test(amphiH,amphiL,var=T)
t.test(glyL,glyH,var=T)

boxplot(phobicH,phobicL,names=c("SD>2","SD<=2"),main="Hydrophobic residues")
boxplot(philicH,philicL,names=c("SD>2","SD<=2"),main="Hydrophilic residues")
boxplot(amphiH,amphiL,names=c("SD>2","SD<=2"),main="Amphipathic residues")
boxplot(glyH,glyL,names=c("SD>2","SD<=2"),main="Glycines residues")

# Get the p-values
p4<-t.test(phobicH,phobicL,var=T)$p.value
p5<-t.test(philicH,philicL,var=T)$p.value
p6<-t.test(amphiH,amphiL,var=T)$p.value
p7<-t.test(glyL,glyH,var=T)$p.value

# Make the Bonferroni correction:
p.valVector<-c(p1,p2,p3,p4,p5,p6,p7)
p.valVector
p.adjust(p.valVector,method="bonferroni",n=length(p.valVector))

# PLOT THE PROPERTIES VS SD

# Join all sd in one variable (same order in graphics)
sdVector<-c(lowSD$SD,highSD$SD)
sdVector

# Join avg hydrophobicity
hydroVector<-c(lowHydro,highHydro)
hydroVector

# Join avg MW
mwVector<-c(lowMW,highMW)
mwVector

# Join avg pI
pIVector<-c(lowPI,highPI)
pIVector

# Join sd types
sdTypes<-c(typeL$SD,typeH$SD)
sdTypes

# Join hydrophobic aa
phobicTypes<-c(typeL$Hydrophobic,typeH$Hydrophobic)
phobicTypes

# Join hydrophilic aa
philicTypes<-c(typeL$Hydrophilic,typeH$Hydrophilic)
philicTypes

# Join amphipathic aa
amphiTypes<-c(typeL$Amphipathic,typeH$Amphipathic)
amphiTypes

# Join glycines aa
glyTypes<-c(typeL$Glycines,typeH$Glycines)
glyTypes

# Scatter plots:

# Hydrophobicity
plot(sdVector,hydroVector,main="Hydrophobicity vs SD",
xlab = "SD",ylab="Average hydrophobicity",pch=19)
abline(lm(hydroVector~sdVector),col="red") # regression line (y~x)

# Mw
plot(sdVector,mwVector,main="MW vs SD",
xlab = "SD",ylab="Molecular weight",pch=19)
abline(lm(mwVector~sdVector),col="red") # regression line (y~x)

# pI
plot(sdVector,pIVector,main="pI vs SD",
xlab = "SD",ylab="Isoelectric point",pch=19)
abline(lm(pIVector~sdVector),col="red") # regression line (y~x)

# Hydrophobic aa
plot(sdTypes,phobicTypes,main="Hydrophobic abundance vs SD",
xlab = "SD",ylab="Hydrophobic abundance",pch=19)
abline(lm(phobicTypes~sdTypes),col="red") # regression line (y~x)

# Hydrophilic aa
plot(sdTypes,philicTypes,main="Hydrophilic abundance vs SD",
xlab = "SD",ylab="Hydrophilic abundance",pch=19)
abline(lm(philicTypes~sdTypes),col="red") # regression line (y~x)

# Amphipathic aa
plot(sdTypes,amphiTypes,main="Amphipathic abundance vs SD",
xlab = "SD",ylab="Amphipathic abundance",pch=19)
abline(lm(amphiTypes~sdTypes),col="red") # regression line (y~x)

# Glycines aa
plot(sdTypes,glyTypes,main="Glycines abundance vs SD",
xlab = "SD",ylab="Glycines abundance",pch=19)
abline(lm(glyTypes~sdTypes),col="red") # regression line (y~x)


