require(stm)
require(purrr)
# Test compute_aggregate_nmi(word2num(y)$corpus,y)
# map words in each document in the list to a sequence of numbers
word2num <- function(y) {
  temp = factor(unlist(y))
  lvs = levels(temp)
  return(list(corpus = lapply(y, function(elm) as.integer(factor(elm,lvs))),
              nwords = length(lvs)))
}

# compute_aggregate_nmi <- function (z_estim, z, method = "NMI"){
#   # require(clue)
#   # tru_labs <- as.cl_hard_partition( unlist(z) )
#   # estim_labs <- as.cl_hard_partition( unlist(z_estim) )
#   # 
#   # temp = as.numeric(cl_agreement(estim_labs, tru_labs, method = method))
#   # # temp[is.nan(temp)] = 0
#   # if (is.nan(temp))
#   #   return(0)
#   # else
#   #   return(temp)
#   adjustedRandIndex(unlist(z), unlist(z_estim))
# }


compute_doc_labels = function(curr_z) {
  K0 <- max(unlist(curr_z))
  return( sapply(1:length(curr_z), function(j) which.max(tabulate(curr_z[[j]], K0))) )
}

create_docword_freq = function(curr_y) {
  W = max(unlist(curr_y))
  do.call(rbind, lapply(1:length(curr_y), function(j) tabulate(curr_y[[j]], nbins = W)))
}

intcorp_to_ldac <- function(intcorp, fname="test.ldac") {
  y_table <- create_docword_freq(intcorp)
  nW <- dim(y_table)[2]
  stm::writeLdac(plyr::alply(y_table, 1, function(row) rbind(1:nW,row)), fname)  
}

# Turn a list of list of words to a list of sentences
docs_from_words <- function(wordlist_corp) {
  purrr::map(wordlist_corp, ~paste(., collapse = " "))
}


write_wordlist_corp_to_file = function(wordlist_corp, fname = 'dai/test.txt') {
  docs = docs_from_words(wordlist_corp)
  file.create(fname)
  lapply(docs, function(doc) write.table( doc, fname, append=T, row.names = F, col.names = F, quote=F))
}
