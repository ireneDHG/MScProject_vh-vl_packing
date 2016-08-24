######################
# Script name: sameDiffeLowSD.R
# 
# Requirements: 	2 txt input files in tabular format.
#		    	First column with a class type "Same" or "Different"
#		    	Second column with sd value
#		    	Headers available as PDBs (1º columnd) and SD (2º column)
#
# Description: 	Performs a t-test over the two files to see the differences 
#			in the sd between the two classes.
#
######################

# Delete all variables
rm(list=ls(all=TRUE))

# Input files
samePDB2<-read.table("samePDB_sd2.txt", header=TRUE, sep="\t")
head(samePDB2)

samePDB1<-read.table("samePDB_sd1.txt", header=TRUE, sep="\t")
head(samePDB1)

samePDB3<-read.table("samePDB_sd3.txt", header=TRUE, sep="\t")
head(samePDB3)

samePDB4<-read.table("samePDB_sd4.txt", header=TRUE, sep="\t")
head(samePDB4)

samePDB5<-read.table("samePDB_sd5.txt", header=TRUE, sep="\t")
head(samePDB5)

# Get the sd values of the two classes for the first file
diff2 <- samePDB2$SD[samePDB2$PDBs == "Different"]
same2 <- samePDB2$SD[samePDB2$PDBs == "Same"]

# Make a t-test 
t.test(diff2,same2,var=T)

# Get the sd values of the two classes for the second file
diff1 <- samePDB1$SD[samePDB1$PDBs == "Different"]
same1 <- samePDB1$SD[samePDB1$PDBs == "Same"]

# Make a t-test 
t.test(diff1,same1,var=T)

# Get the sd values of the two classes for the third file
diff3 <- samePDB3$SD[samePDB3$PDBs == "Different"]
same3 <- samePDB3$SD[samePDB3$PDBs == "Same"]

# Make a t-test 
t.test(diff3,same3,var=T)

# Get the sd values of the two classes for the fourth file
diff4 <- samePDB4$SD[samePDB4$PDBs == "Different"]
same4 <- samePDB4$SD[samePDB4$PDBs == "Same"]

# Make a t-test 
t.test(diff4,same4,var=T)

# Get the sd values of the two classes for the fifth file
diff5 <- samePDB5$SD[samePDB5$PDBs == "Different"]
same5 <- samePDB5$SD[samePDB5$PDBs == "Same"]

# Make a t-test 
t.test(diff5,same5,var=T)
