library(Matrix)
library(purrr)
library(extraDistr)
library(tidyverse)

source("helpers.R")
source("net_common.R")
examp_vocab <- read_csv("data/vocab.nips.txt", col_names = F) %>% pull(X1)

# K = 3 # number of clusters
gen_doc_data = function(J, K = 3, n = 100, W=15, alpha0=3) {
  tru_doc_label = sample(K, J, replace=T) # true document labels
  # word_label_dist = 0.7*label_vec2mat(tru_doc_label) + 0.3*Matrix(rdirichlet(J,rep(1,K))) # sample word distributions among labels in each document
  word_label_dist = 0.7*label_vec2mat(tru_doc_label) + 0.3*rdirichlet(J,rep(1,K)) # sample word distributions among labels in each document  
  # image(Matrix(word_label_dist))

    # z = lapply(1:J, function(j) label_mat2vec(t(rmultinom(n, 1, prob=word_label_dist[j,]))))
  z = lapply(1:J, function(j) sample(K, n, replace=T,  prob=word_label_dist[j,])) # sample word label in each document
  shifts = seq(1,0.9*W,length.out = K)
  shapes = purrr::map(1:K, ~dnorm(1:W,shifts[.],W/6))
  shapes = purrr::map(1:K, ~(shapes[[.]])/sum(shapes[[.]]))
  shape_mat = do.call(rbind, shapes)
  # plot(shape_mat[3,])
  # phi = rdirichlet(K, rep(1/W,W))
  phi = rdirichlet(K, alpha0*shape_mat)
  # image(Matrix(phi))
  y = lapply(1:J, function(j) rcat(n, phi[z[[j]],]))
  y = word2num(y)$corpus
  vocab_subset = sample(examp_vocab, W)
  return( list(corpus1 = y, corpus2 = lapply(y, function(elem) vocab_subset[elem]),
               word_label = z, doc_label = tru_doc_label, phi=phi, vocab = vocab_subset) )
}

# library(reticulate)

# out = gen_doc_data(3, n=100, W=20)
# write_wordlist_corp_to_file(out$corpus2)
# docs = out$corpus1
# saveRDS(docs,"temp.rds")

# reticulate::source_python("dai/dai_hdp.py")
# docs = readRDS("temp.rds")
# doc_topic_freq = dai_warpper(docs, 5L)
# lapply(doc_topic_freq, function(freq_table) max.col(freq_table))

# out = gen_doc_data(3, n= 20)
# vocab = out$vocab
# library(topicmodels)
# library(tm)
# library(lda)
# 
# lexicalize(out$corpus2)


# corp = VCorpus( VectorSource(docs_from_words(out$corpus2)) )
# dtm = DocumentTermMatrix(corp)
# dtm2ldaformat(dtm)
