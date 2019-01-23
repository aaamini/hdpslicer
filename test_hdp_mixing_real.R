library(plotly)
plot1 <-function(y) plot_ly(x=seq(1,length(y)),y=~y)

source('hdp_module.R')

library(MASS)
library(mvtnorm)
library(clue)

compute_aggregate_nmi <- function (z_estim, z){
  require(clue)
  tru_labs <- as.cl_hard_partition( unlist(z) )
  estim_labs <- as.cl_hard_partition( unlist(z_estim) )
  
  cl_agreement(estim_labs, tru_labs, method = 'NMI') 
}


library(data.table)
data_filename <- "data/docword.nips.txt"
temp <- unlist(read.csv(data_filename, nrows=3, header = F))
(D <- temp[1])
(W <- temp[2])
(N <- temp[3])

vocab <- fread("data/vocab.nips.txt")
#dat <- read.csv(file = "data/docword.nips.txt", skip = 3)
data_full <- fread(file = data_filename, skip = 3, col.names = c("docID","wordID","count"))
summary(dat_full)

top_words <- data_full[, sum(count), by=wordID]

top_wID <- top_words[order(-V1)]$wordID[1:50]



J <- 5
dat <- data_full[wordID %in% top_wID, ][docID %in% 1:J,]
dat
#dat <- data_full[docID %in% 1:J,]
unique(dat[,2])



library(splitstackshape)
dat.ex <- expandRows(dat,"count")
y <- lapply( split(dat.ex, dat.ex$docID, flatten=T), function(x) x$wordID ) 
str(y)

ITRmax <- 1000
zh <- hdp_slice_sampler(y, beta0=3, gam0=1, 
                        ITRmax=ITRmax, Kcap=50, Tcap=50,
                        categorical=T, W=W)

  
nmi <- sapply(3:ITRmax, function(itr) compute_aggregate_nmi(zh[[itr]], zh[itr-1]))

p <- plot_ly(x=3:ITRmax, width = 600, height = 400) %>% 
  #add_markers(y = nmi[[1]], name = 'n = 30', marker=list(symbol="diamond-open-dot")) %>%
  add_markers(y = nmi, name = 'n = 100') %>%
  #add_markers(y = nmi[[3]], name = 'n = 300',marker=list(symbol="x")) %>%
  #add_trace(y = nmi[[1]], name = 'trace 1', mode = 'lines+markers') %>%
  layout(xaxis = list(title="iteration"), yaxis=list(title="NMI")) %>%
  layout(legend = list(x = 0.7, y = 0.1)) 

p

export(p, file=paste("nmi", fname, ".pdf"))
export(p, file=paste("nmi", fname, ".png"))

p <- plot_ly(nmi_df, x=2:ITRmax, y=nmi) %>% layout(xaxis = list(title="iteration"), yaxis=list(title="NMI"))
p

#pdf("nmi_1.png")
#print(p)
#dev.off()

library(hdp)
example_data_hdp

set.seed(10)
quick_chain <- hdp_posterior(out$y, burnin=1, n=50, space=1, seed=1234)
temp <- hdp_extract_components(quick_chain)
str(temp)
plot_comp_size(temp, bty="L", lab=c(3, 5, 7))
base(quick_chain)

#plot_lik(quick_chain, bty="L")
# cRand <- cl_agreement(estim_labs, tru_labs, method = "cRand")
# acc <- cl_agreement(estim_labs, tru_labs, method = "diag")
# cat("NMI     = ", round(nmi,3), "\nacc.    = ", round(acc,3), "\nadjRand = ", round(cRand,3), "\n")

