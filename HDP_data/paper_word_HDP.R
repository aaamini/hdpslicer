setwd("~/Dropbox/niloo/m/HDP_data")
library(data.table)
library(tidyr)
library(splitstackshape)
library(igraph)


ML_data<- read.csv("ML.csv",stringsAsFactors = F)
ML_data$authors <- gsub("\\[", "", ML_data$authors)
ML_data$authors <- gsub("\\]", "", ML_data$authors)
ML_data11<- cSplit(ML_data, 'authors', ',', 'long')
data11 <- cbind(ML_data11,"ML")
colnames(data11)[3] <- 'label'
paper_author_1<- unique(cbind(data11[,2],data11[,3]))[1:300]
#####################################################################
security_data<- read.csv("security.csv",stringsAsFactors = F)
security_data$authors <- gsub("\\[", "", security_data$authors)
security_data$authors <- gsub("\\]", "", security_data$authors)
security_data11<- cSplit(security_data, 'authors', ',', 'long')
data22 <- cbind(security_data11,"security")
colnames(data22)[3] <- 'label'
paper_author_2<- unique(cbind(data22[,2],data22[,3]))[1:300]
###################################################################
multimedia_data<- read.csv("multimedia.csv",stringsAsFactors = F)
multimedia_data$authors <- gsub("\\[", "", multimedia_data$authors)
multimedia_data$authors <- gsub("\\]", "", multimedia_data$authors)
multimedia_data33<- cSplit(multimedia_data, 'authors', ',', 'long')
data33 <- cbind(multimedia_data33,"multimedia")
colnames(data33)[3] <- 'label'
paper_author_3<- unique(cbind(data33[,2],data33[,3]))[1:300]

##############################################################


##################################################
papers_and_labels <- as.data.frame(unique(rbind(paper_author_1,paper_author_2,paper_author_3)))


library(tm)
library(stringi)
library(proxy)

# load data ---------------------------------------------------------------

papers <-papers_and_labels$title
docs <- Corpus(VectorSource(papers))

docs2 <- tm_map(docs, PlainTextDocument)
docs2 <- tm_map(docs2, stripWhitespace)
#docs2 <- tm_map(docs2, removePunctuation)
docs2 <- tm_map(docs2, removeNumbers)
docs2 <- tm_map(docs2, content_transformer(tolower))
docs2 <- tm_map(docs2, stemDocument, language="english")
docs2 <- tm_map(docs2, removeWords,stopwords("english"))



tdm <- TermDocumentMatrix(docs2)
tdm4 <- removeSparseTerms(tdm, 0.993)
tdm.mat4 <- as.matrix(tdm4)
dim(tdm.mat4)
id4 <- which(colSums(tdm.mat4)==0)
tdm.mat5<- tdm.mat4[,-id4]
papers_label <-papers_and_labels$label[-id4]


setwd("~/Dropbox/niloo/m/HDP_data/")
save(tdm.mat5,papers_label, file="paper_words_HDP.RData")
