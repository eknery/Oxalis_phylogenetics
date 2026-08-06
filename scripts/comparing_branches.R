### libraries
if(!require("tidyverse")) install.packages("tidyverse"); library("tidyverse")
if(!require("ggplot2")) install.packages("ggplot2"); library("ggplot2")
if(!require("phangorn")) install.packages("phangorn"); library("phangorn")
if(!require("ape")) install.packages("ape"); library("ape")
if(!require("seqinr")) install.packages("seqinr"); library("seqinr")
if(!require("gamlss")) install.packages("gamlss"); library("gamlss")
if(!require("ggeffects")) install.packages("ggeffects"); library("ggeffects")

### file names
dir_input = "7_ml_phylogenies/bootstrap/"
### files in directory
file_names = list.files(dir_input)
### loading data
tree_list = list()
for(i in 1:length(file_names) ){
  tree_name = file_names[i]
  tree_list[[i]] = read.tree(file = paste0(dir_input, tree_name))
  names(tree_list)[i] =  str_remove(string = tree_name, 
                                    pattern = ".tre")
}

############################### PROCESSING DATA ################################

### getting species per locus
all_names = c()
for(i in 1:length(tree_list)){
  some_names = tree_list[[i]][[1]]$tip.label
  all_names = c(all_names, some_names)
}
### into one dataframe
all_names = sort(unique(all_names))

### get names per locus
names_loci = all_names
for(i in 1:length(tree_list)){
  boll_names = all_names %in% tree_list[[i]][[1]]$tip.label
  names_loci = cbind(names_loci, boll_names)
}
### transform to tibble
names_loci = as_tibble(names_loci)

### get species with all loci sequenced 
common_names = names_loci %>% 
  filter_at(vars(-names_loci), all_vars(. == TRUE) ) %>% 
  select(names_loci) %>% 
  pull()

### only Oxalis names
common_names = common_names[grepl("Oxalis_", common_names)]

### pruning trees to sampled species
pruned_tree_list = tree_list
for (i in 1:length(tree_list) ){
  pruned_trees = tree_list[[i]]
  for(j in 1:length(tree_list[[i]]) ){
    pruned_trees[[j]] = keep.tip(phy = tree_list[[i]][[j]], 
                                 tip = common_names)
    
  }
  pruned_tree_list[[i]] = pruned_trees
}

########################### CALCULATING BRANCH LENGTHS #########################

### vectors for results
locus = c()
branch_len = c()
### loop over prunned trees
for (i in 1:length(pruned_tree_list) ){
  for(j in 1:length(pruned_tree_list[[i]]) ){
    locus = c(locus, names(pruned_tree_list[i]) )
    branch_len = c(branch_len, sum(pruned_tree_list[[i]][[j]]$edge.length) )
  }
}

### dataframe 
branch_df = data.frame(
  "loucs" = locus,
  "branch_len" = branch_len
)

###
branch_locus = ggplot(
  data = branch_df, 
  aes(x= locus, y= as.numeric(branch_len), fill= locus)
) +
  geom_boxplot(
    width = 0.4, 
    outlier.shape = NA,
    alpha = 0.5
  )+
  scale_fill_manual(
    values = c(
      "nuclear"= "darkorange",
      "plastidial" = "darkgreen"
    )
  ) +
  labs(
    x = "Partition",
    y = "Tree length",
  ) +
  theme(panel.background=element_rect(fill="white"), 
        panel.grid=element_line(colour=NULL),
        panel.border=element_rect(fill=NA,colour="black"),
        axis.title=element_text(size=12,face="bold"),
        axis.text.x=element_text(size=10),
        legend.position = "none")

tiff("figures/branch_partition.tiff",
     units="cm", 
     width= 9.5, 
     height= 8.5,
     res=1200)
  branch_locus
dev.off()
