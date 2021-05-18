#set the working directory

#set repeat type
#-dn	Search for direct non-complementary repeats
#-dc	Search for direct complementary repeats
#-in	Search for inverse non-complementary repeats
#-ic	Search for inverse complementary repeats
repeat_type = '\'dn\''

#retrieve repeats data from ibmdb2 database
#set up driver
#user can use other type of connections/drivers as well
#user can use RData file with already prepared data for the analysis
#install.packages("RODBC")
library(RODBC)

#set database name
dsn.name <- "MITOHOND"
#dsn.name <- "RIPITIV1"

con1 <- odbcConnect(dsn=dsn.name)

#pull details about repeats frequencies
queryPatternFreq <- gsub("[\r\n\t]", " ", paste("WITH LOCATIONS AS
					(
					SELECT SEQID, ID_FRAGMENT1, LOCATION1 AS LOCATION
					FROM DB2ADMIN.MATCH JOIN DB2ADMIN.FRAGMENT ON DB2ADMIN.MATCH.ID_FRAGMENT1 = DB2ADMIN.FRAGMENT.ID
  				    JOIN  DB2ADMIN.SEQUENCE ON SEQID = DB2ADMIN.SEQUENCE.ID
  				    WHERE DB2ADMIN.SEQUENCE.TYPE = ", repeat_type , "
					UNION 
					SELECT SEQID, ID_FRAGMENT1,LOCATION2 AS LOCATION
					FROM DB2ADMIN.MATCH JOIN DB2ADMIN.FRAGMENT ON DB2ADMIN.MATCH.ID_FRAGMENT1 = DB2ADMIN.FRAGMENT.ID
					JOIN  DB2ADMIN.SEQUENCE ON SEQID = DB2ADMIN.SEQUENCE.ID
  				    WHERE DB2ADMIN.SEQUENCE.TYPE = ", repeat_type , "  
            		)
						SELECT DISTINCT SEQID, ID_FRAGMENT1, COUNT(*) AS PF
						FROM LOCATIONS
						GROUP BY SEQID, ID_FRAGMENT1;", sep = " "))
resPatternFreq <- sqlQuery(con1, queryPatternFreq, errors=FALSE)

#pull details about repeats position
queryPatternPositions <- gsub("[\r\n\t]", " ", paste("WITH LOCATIONS AS
					(
					SELECT SEQID, ID_FRAGMENT1, LOCATION1 AS LOCATION
					FROM DB2ADMIN.MATCH JOIN DB2ADMIN.FRAGMENT ON DB2ADMIN.MATCH.ID_FRAGMENT1 = DB2ADMIN.FRAGMENT.ID
           			JOIN  DB2ADMIN.SEQUENCE ON SEQID = DB2ADMIN.SEQUENCE.ID
  				    WHERE DB2ADMIN.SEQUENCE.TYPE = ", repeat_type , "
					UNION 
					SELECT SEQID, ID_FRAGMENT1,LOCATION2 AS LOCATION
					FROM DB2ADMIN.MATCH JOIN DB2ADMIN.FRAGMENT ON DB2ADMIN.MATCH.ID_FRAGMENT1 = DB2ADMIN.FRAGMENT.ID
           			JOIN  DB2ADMIN.SEQUENCE ON SEQID = DB2ADMIN.SEQUENCE.ID
  		    		WHERE DB2ADMIN.SEQUENCE.TYPE = ", repeat_type , "
					)
						SELECT DISTINCT SEQID, ID_FRAGMENT1, LOCATION 
						FROM LOCATIONS;", sep = " "))
resPatternPositions <- sqlQuery(con1, queryPatternPositions, errors=FALSE)

#pull details about sequence length
querySeqLen <- gsub("[\r\n\t]", " ", paste("SELECT ID,  SEQLENGTH
					FROM DB2ADMIN.SEQUENCE
					WHERE DB2ADMIN.SEQUENCE.TYPE = ", repeat_type , ";", sep = " "))
resSeqLen <- sqlQuery(con1, querySeqLen, errors=FALSE)

#create Vocabular - the set pf repeats
queryVocabular <- gsub("[\r\n\t]", " ", paste(" SELECT DISTINCT DB2ADMIN.FRAGMENT.ID, TEXT, TEXT_LENGTH
				    FROM DB2ADMIN.FRAGMENT LEFT JOIN DB2ADMIN.MATCH ON DB2ADMIN.MATCH.ID_FRAGMENT1 = DB2ADMIN.FRAGMENT.ID
				    LEFT JOIN  DB2ADMIN.SEQUENCE ON SEQID = DB2ADMIN.SEQUENCE.ID
  				    WHERE DB2ADMIN.SEQUENCE.TYPE = ", repeat_type , "
				    ORDER BY TEXT_LENGTH, TEXT", sep = " "))
resVocabular <- sqlQuery(con1, queryVocabular, errors=FALSE)

#pull details about sequences
querySeq <- gsub("[\r\n\t]", " ", paste(" SELECT ID
					FROM DB2ADMIN.SEQUENCE  WHERE DB2ADMIN.SEQUENCE.TYPE = ", repeat_type , "", sep = " "))
resSeq <- sqlQuery(con1, querySeq, errors=FALSE)

#pull details about host mith
queryHost <- gsub("[\r\n\t]", " ", paste("SELECT ID, SUBSTRING(NAME, 0,LOCATE('.', NAME)+2) as Name,  HOST
					FROM DB2ADMIN.SEQUENCE WHERE DB2ADMIN.SEQUENCE.TYPE = ", repeat_type , "", sep = " "))
resHost <- sqlQuery(con1, queryHost, errors=FALSE)

#pull details about sequences
queryNames <- gsub("[\r\n\t]", " ", paste("  SELECT ID, SUBSTRING(NAME, 0,LOCATE('.', NAME)+2)
					FROM DB2ADMIN.SEQUENCE WHERE DB2ADMIN.SEQUENCE.TYPE = ", repeat_type , "", sep = " "))
resNames <- sqlQuery(con1, queryNames, errors=FALSE)


#pull details about host viruses
#queryHost <- gsub("[\r\n\t]", " ", paste("SELECT ID, LOCATE('|', NAME)-4,  SUBSTRING(NAME, LOCATE(':', NAME)+1,LOCATE('|', NAME)-LOCATE(':', NAME)-1) as Name,  VIRUS_SPECIES AS HOST
#  FROM DB2ADMIN.SEQUENCE WHERE DB2ADMIN.SEQUENCE.TYPE = ", repeat_type , "", sep = " "))

#resHost <- sqlQuery(con1, queryHost, errors=FALSE)

#pull details about sequences
#queryNames <- gsub("[\r\n\t]", " ", paste("  SELECT ID, LOCATE('|', NAME)-4,  SUBSTRING(NAME, LOCATE(':', NAME)+1,LOCATE('|', NAME)-LOCATE(':', NAME)-1)
#  FROM DB2ADMIN.SEQUENCE WHERE DB2ADMIN.SEQUENCE.TYPE = ", repeat_type , "", sep = " "))

#resNames <- sqlQuery(con1, queryNames, errors=FALSE)


# Close connections
odbcCloseAll()
cat("Database connections are closed.\n")

seqVectorsF_P <- data.frame(matrix(0, ncol = nrow(resSeq), nrow = nrow(resVocabular))  ) 
rownames(seqVectorsF_P) <- resVocabular[]$ID
colnames(seqVectorsF_P) <- resSeq[]$ID
dim(seqVectorsF_P)

totalSeq <- nrow(resSeq)
totalVoc <- nrow(resVocabular)

save.image("RData")