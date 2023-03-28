### Figure S5
file1 <- '04_model_raw_RData/single_LODO/02lasso/'
load('./01_data/cohort_info.RData')

# data was in 04_model_raw_RData/single_LODO/02lasso

for (i_t in 1:nrow(cross)){
  m <- cross$markers[i_t]
  print(cross$disease[i_t])
  print(m)
  load(paste0(file1,m))
  
  
  marker_4E <- subset(marker.adj[["marker_data"]],class=='adjust')
  marker_4E$LDA <- as.numeric(marker_4E$LDA)
  l_t<- unique(marker_4E$scientific_name)
  marker_4E_factor <- data.frame(p=rep(0,length(l_t)),n=rep(0,length(l_t)),row.names =l_t)
  
  for (i in unique(marker_4E$scientific_name)){
    marker_4E_factor[i,'p'] <- nrow(subset(marker_4E,(scientific_name%in%i)&(LDA>0)))
    marker_4E_factor[i,'n'] <- nrow(subset(marker_4E,(scientific_name%in%i)&(LDA<0)))
  }
  
  
  marker.adj[["marker_data"]]$LDA <- as.numeric(marker.adj[["marker_data"]]$LDA)
  marker.adj[["marker_data"]]$nrproj <- as.numeric(marker.adj[["marker_data"]]$nrproj)
  marker.adj[["marker_data"]] <- arrange(marker.adj[["marker_data"]],LDA)
  my_marker.plot(marker.adj[["marker_data"]],feat_list=NULL,meta_list=NULL,lda_cutoff=2,nproj_cutoff=1,level='genus')
  
  # exclude different
  marker_4E_factor <- rownames(marker_4E_factor)[marker_4E_factor$p!=0&marker_4E_factor$n!=0]
  marker.adj[["marker_data"]] <- subset(marker_4E,!scientific_name %in% marker_4E_factor)
  marker.adj[["marker_data"]]$LDA <- as.numeric(marker.adj[["marker_data"]]$LDA)
  marker.adj[["marker_data"]]$nrproj <- as.numeric(marker.adj[["marker_data"]]$nrproj)
  marker.adj[["marker_data"]] <- arrange(marker.adj[["marker_data"]],LDA)
  my_marker.plot(marker.adj[["marker_data"]],feat_list=NULL,meta_list=NULL,lda_cutoff=2,nproj_cutoff=1,level='genus')
}


