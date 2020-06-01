require(stm)
require(purrr)


# Change point detection --------------------------------------------------

# Try it with: 
# breaks = c(0,1,2,69,100)
intervals_from_breaks = function(breaks) {
  m = length(breaks)
  idx = c(1, rep(2:(m-1), each=2), m)
  Intervals( matrix(breaks[idx], ncol=2, byrow = T) )
}

detect_change = function(y) {
  # require(intervals)
  # require(purrr)
  n = length(y)
  # interp = approx(1:n, y, n=5*length(y))
  out = trendsegment(y)
  
  breaks = c(0, out$cpt, n)
  ints = intervals_from_breaks(breaks) # split into intervals according to break points
  # new_breaks = c(0, as.matrix(contract(ints, delta=1))[,1], n) # merge short intervals
  new_breaks = c(as.matrix(contract(ints, delta=1))[,1]-1, n) # merge short intervals
  # breaks = breaks[c(T, diff(breaks) > 2)]  # remove short intervals
  
  interval_sds = unlist( purrr::map( split(y, cut(seq_along(y), new_breaks)), sd) )
  return( new_breaks[ which.max(interval_sds) + 1 ] )
}



# Documment manipulation --------------------------------------------------

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
which_most_frequent = function(vec) {
  as.integer(names(which.max(table(vec))))
}

# input should be a list otherwise the behavior is unintended
compute_doc_labels = function(label_vec_list) {
  # label_vec_list: a list of label vectors -> we return which label in each vector is most frequent
  # K0 <- max(unlist(curr_z))
  # return( sapply(1:length(curr_z), function(j) which.max(tabulate(curr_z[[j]], K0))) )
  return( sapply(label_vec_list, which_most_frequent) )
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
