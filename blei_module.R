require(lpSolve)
# create or empty directory
recreate_dir = function(dir) {
  unlink(dir, recursive = T) # remove the directory
  dir.create(dir)
}

blei_splitmerge_wrapper = function(curr_y, max_iter=20, split_merge=TRUE, remove_dir=TRUE) {
  bleid = list(indir = "blei_in", outdir = "blei_out")
  lapply(bleid, recreate_dir)  
  
  ## Create ldac data file
  for (i in 1:length(curr_y)) {
    dat = paste(curr_y[[i]], collapse = " ")
    write(dat, file.path(bleid$indir, paste0("doc", i, ".txt")))
  }
  system(paste("python2 text2ldac.py",bleid$indir))          ## for Linux
  datfile_path = file.path(bleid$indir, paste0(bleid$indir, ".dat"))
  
  ## Run Blei's
  # max_iter = 20
  # split_merge = TRUE
  split_merge_flag = ifelse(split_merge, " --split_merge yes", "")
  cmd = sprintf("blei_hdp/hdp/hdp --algorithm train --data %s --max_iter %s --save_lag 1 --directory %s%s", 
                datfile_path, max_iter, bleid$outdir, split_merge_flag)
  system(cmd)
  
  ## process the output 
  flist = list.files(path=bleid$outdir, pattern="word")
  flist_pieces = strsplit(flist, "-")
  
  estim_doc_labels = vector("list",max_iter)
  for (i in 1:length(flist)) {
    itr_num = flist_pieces[[i]][1]  # something like "00000" "00001" or "mode"
    if (is.na(suppressWarnings(as.numeric(itr_num))))  { # if == "mode" ignore
      # print(as.integer(itr_num))
      next
    }
    itr_num = as.integer(itr_num)+1
    
    y_table <- create_docword_freq(curr_y)  # create the doc_word table table
    nwords = dim(y_table)[2]
    
    # temp = readr::read_delim(file.path(bleid$outdir, flist[[i]]), delim=" ")
    temp = read.delim(file.path(bleid$outdir, flist[[i]]), sep=" ")
    temp = temp + 1 # 0- to 1-based indexing
    
    doc_list =  temp %>% group_by(d) %>% group_split()
    out_table = t( sapply(doc_list, function(x) tabulate(x$w, nwords)) )
    
    zlabels = purrr::map(doc_list, ~.$z)
    temp_doc_labels = compute_doc_labels(zlabels)
    
    # Matching dociments 
    # They output document and words in random order
    # Solve linear assignement problem to recover document permutation
    b = t(apply(y_table, 1, sort))
    a = t(apply(out_table, 1, sort))
    costs = -a %*% t(b)
    Pmat = lp.assign(costs)$solution
    # Pmat = doc i in the output corresponds to document j in the input where Pmat[i,j] = 1
    # Pmat %*% 1:J moves document i to position j where Pmat[i,j] = 1
    
    # permute documents to original order
    estim_doc_labels[[itr_num]] <- as.vector(t(Pmat) %*% temp_doc_labels)  # t(Pmat) is the inverse permutation matrix
  }
  
  if (remove_dir) lapply(bleid, function(dir) unlink(dir, recursive = T))
  return(estim_doc_labels)
}






##############  running blei package  ############## 
 
##  performs posterior sampling
# for (i in 1:n_samps) {
#   ##  create folders to store sampling outputs
#   train_dir_s = paste0("comp2/split/train_dir_", i)      ##  folder for sampling with split-merge
#   train_dir_ns = paste0("comp2/no_split/train_dir_", i)    ##  folder for sampling without split-merge
#   if (!dir.exists(train_dir_s)) {
#     dir.create(paste0(train_dir_s), recursive = T)
#   }
#   if (!dir.exists(train_dir_ns)) {
#     dir.create(paste0(train_dir_ns), recursive = T)
#   }
#   
#   # ## For Windows
#   # command_s = paste0("wsl ./hdp --algorithm train --data text2ldac/blei_input/blei_input.dat --max_iter ", n_iters," --save_lag 1 --directory ", train_dir_s, " --split_merge yes")
#   # command_ns = paste0("wsl ./hdp --algorithm train --data text2ldac/blei_input/blei_input.dat --max_iter ", n_iters," --save_lag 1 --directory ", train_dir_ns)
#   # shell(command_s)
#   # shell(command_ns)
# 
#   ##  For Linux (?)
#   command_s = paste0("./hdp --algorithm train --data text2ldac/blei_input/blei_input.dat --max_iter ", n_iters," --save_lag 1 --directory ", train_dir_s, " --split_merge yes")
#   command_ns = paste0("./hdp --algorithm train --data text2ldac/blei_input/blei_input.dat --max_iter ", n_iters," --save_lag 1 --directory ", train_dir_ns)
#   system(command_s)
#   system(command_ns)
# }




# test
# P = diag(5)[sample(5),]
# P %*% y_table
# b = t(apply(y_table, 1, sort))
# a = t(apply(P %*% y_table, 1, sort))
# costs = -a %*% t(b)
# lp.assign(costs)$solution








# ##############  rematching documents  ############## 
# 
# ##  match_d is a data frame with doc.in as input document id and doc.out as output document id
# # match_d = quick_rematch(J, paste0(hdp_dir, "/comp2"), paste0(text2ldac_dir, "/blei_input"), W = W)
# match_d = quick_rematch(J, bleid$indir, blied$outdir, W = W)
# 
# ##############  finding cRand  ############## 
# 
# ##  compute the cRand for each iteration
# nmi_mat_s = matrix(0, nrow = n_samps, ncol = n_iters)
# nmi_mat_ns = matrix(0, nrow = n_samps, ncol = n_iters)
# for (j in 1:n_samps) {
#   ## returns the topic for each document
#   out_s = find_topics(paste0(hdp_dir, "/comp2/split"), set = j, doc_labels, match_d, n_iters, 1)      ##  document topic for split-merge output
#   out_ns = find_topics(paste0(hdp_dir, "/comp2/no_split"), set = j, doc_labels, match_d, n_iters, 1)   ##  document topic for non split-merge output
#   #nmi_mat_s[j,] = sapply(1:n_iters, function(itr) compute_aggregate_nmi(compute_doc_labels(out_s[[1]][[itr]]), doc_labels, method="cRand")) 
#   #nmi_mat_ns[j,] = sapply(1:n_iters, function(itr) compute_aggregate_nmi(compute_doc_labels(out_ns[[1]][[itr]]), doc_labels, method="cRand")) 
#   nmi_mat_s[j,] = sapply(1:n_iters, function(itr) compute_mutual_info(compute_doc_labels(out_s[[itr]]), doc_labels)) 
#   nmi_mat_ns[j,] = sapply(1:n_iters, function(itr) compute_mutual_info(compute_doc_labels(out_ns[[itr]]), doc_labels)) 
# }
# 
# 
# plot(apply(nmi_mat_s, 2, mean), ylim = c(0.5, 1))
# plot(apply(nmi_mat_ns, 2, mean), ylim = c(0.5, 1))