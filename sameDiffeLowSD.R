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

samePDB5<-read.table("samePDB_sd5.txt", header=TRUE, sep="\t")
head(samePDB5)

samePDB10<-read.table("samePDB_sd10.txt", header=TRUE, sep="\t")
head(samePDB10)

samePDB40<-read.table("samePDB_sd40.txt", header=TRUE, sep="\t")
head(samePDB40)

samePDB60<-read.table("samePDB_sd60.txt", header=TRUE, sep="\t")
head(samePDB60)

# Get the sd values of the two classes for the 1 file
diff2 <- samePDB2$SD[samePDB2$PDBs == "Different"]
same2 <- samePDB2$SD[samePDB2$PDBs == "Same"]

# Make a t-test 
t.test(diff2,same2,var=T)

# Get the sd values of the two classes for the 2 file
diff1 <- samePDB1$SD[samePDB1$PDBs == "Different"]
same1 <- samePDB1$SD[samePDB1$PDBs == "Same"]

# Make a t-test 
t.test(diff1,same1,var=T)


# Get the sd values of the two classes for the 3 file
diff5 <- samePDB5$SD[samePDB5$PDBs == "Different"]
same5 <- samePDB5$SD[samePDB5$PDBs == "Same"]

# Make a t-test 
t.test(diff5,same5,var=T)

# Get the sd values of the two classes for the 4 file
diff10 <- samePDB10$SD[samePDB10$PDBs == "Different"]
same10 <- samePDB10$SD[samePDB10$PDBs == "Same"]

# Make a t-test 
t.test(diff10,same10,var=T)

# Get the sd values of the two classes for the 5 file
diff40 <- samePDB40$SD[samePDB40$PDBs == "Different"]
same40 <- samePDB40$SD[samePDB40$PDBs == "Same"]

# Make a t-test 
t.test(diff40,same40,var=T)

# Get the sd values of the two classes for the 6 file
diff60 <- samePDB60$SD[samePDB60$PDBs == "Different"]
same60 <- samePDB60$SD[samePDB60$PDBs == "Same"]

# Make a t-test 
t.test(diff60,same60,var=T)
