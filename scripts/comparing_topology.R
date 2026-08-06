### load libraries
if(!require("tidyverse")) install.packages("tidyverse"); library("tidyverse")
if(!require("ggplot2")) install.packages("ggplot2"); library("ggplot2")
if(!require("phangorn")) install.packages("phangorn"); library("phangorn")
if(!require("ape")) install.packages("ape"); library("ape")
if(!require("seqinr")) install.packages("seqinr"); library("seqinr")

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

################################# FIND COMMON SPECIES ###########################

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
### column names
colnames(names_loci) = c("names_loci", names(tree_list))
### transform to tibble
names_loci = as_tibble(names_loci)

### get species with all loci sequenced 
common_names = names_loci %>% 
  filter_at(vars(-names_loci), all_vars(. == TRUE) ) %>% 
  select(names_loci) %>% 
  pull()

### only Oxalis names
common_names = common_names[grepl("Oxalis_", common_names)]

################################ PROCESSING TREES #############################

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

### transforming in vector
pruned_trees_vec = c()
for(i in 1:length(pruned_tree_list)){
  pruned_trees_vec = c(pruned_trees_vec, pruned_tree_list[[i]] )
}


################################ CONGRUENCE ANALYSES ##########################

### distance
dist = RF.dist(
  c(pruned_tree_list[[1]],  
    pruned_tree_list[[2]]
    ),
  tree2 = NULL, 
  normalize = T
)

### getting vector with names for each distance
locus = c()
for(i in 1:length(pruned_tree_list)){
  locus_name = names(pruned_tree_list)[i]
  locus_rep = rep(locus_name, length.out= length(pruned_tree_list[[i]]) )
  locus = c(locus, locus_rep)
}

### PCOA
pcoa = pcoa(dist, correction="none", rn=NULL)
pcoa_df = as.data.frame(cbind(locus, pcoa$vectors))
### get % var 
pc_rel_var = pcoa$values$Relative_eig
### axis names
pc_axis_1 = paste0("PCoA (", round(pc_rel_var[1]*100, 2), "%)" )
pc_axis_2 = paste0("PCoA (", round(pc_rel_var[2]*100, 2), "%)" )

### plot pcoa
pcoa_plot = ggplot(
  data = pcoa_df,
  aes(
    x=as.numeric(Axis.1),
    y=as.numeric(Axis.2),
    color=locus)
  ) +
  geom_point(size = 1, alpha = 0.5) +
  scale_colour_manual(values=c(
    "nuclear"= "darkorange",
    "plastidial" = "darkgreen"
    )
  )+
  labs(
    x= pc_axis_1, 
    y= pc_axis_2
  )+
  theme(panel.background=element_rect(fill="white"), 
        panel.grid=element_line(colour=NULL),
        panel.border=element_rect(fill=NA,colour="black"),
        axis.title=element_text(size=12,face="bold"),
        axis.text.x=element_text(size=10),
        legend.position = "none")

### export plot
tiff("figures/pcoa_ml_trees.tiff", 
     units="cm", width=10, height=9, res=600)
 pcoa_plot
dev.off()
