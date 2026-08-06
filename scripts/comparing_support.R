### libraries
if(!require("tidyverse")) install.packages("tidyverse"); library("tidyverse")
if(!require("ggplot2")) install.packages("ggplot2"); library("ggplot2")
if(!require("phangorn")) install.packages("phangorn"); library("phangorn")
if(!require("ape")) install.packages("ape"); library("ape")
if(!require("seqinr")) install.packages("seqinr"); library("seqinr")
if(!require("gamlss")) install.packages("gamlss"); library("gamlss")
if(!require("ggeffects")) install.packages("ggeffects"); library("ggeffects")

### load Ml phylogenetic tree
mcc = read.tree(file = "8_bayesian_phylogenies/mcc_median.tre")

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

### pruning MCC tree
mcc_pruned = keep.tip(mcc, common_names)

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

############################## CALCULATING SUPPORT #############################

### support from nuclear
mcc_nuc = plotBS(
  tree = mcc_pruned, 
  trees =  pruned_tree_list[[1]], 
  type = "phylogram", 
  method = "TBE"
  )
### get nuclear support 
nuc_supp = mcc_nuc$node.label

### support from plastidial
mcc_pla = plotBS(
  tree = mcc_pruned, 
  trees =  pruned_tree_list[[2]], 
  type = "phylogram", 
  method = "TBE"
)
### get nuclear support 
pla_supp = mcc_pla$node.label

### locus vector
locus = c(rep("nuclear", length(nuc_supp) ), 
          rep("plastidial", length(pla_supp) )
          )
### support vector
support = c(nuc_supp, pla_supp)

### node ages
age = branching.times(mcc_pruned)
age = c(age, age)

### into dataframe
supp_df = data.frame(locus, support, age)

### removing NAs
supp_df = supp_df[!is.na(supp_df$support),]

################################# ANALYZING SUPPORT ############################

### https://www.rdocumentation.org/packages/gamlss.dist/versions/6.1-1/topics/gamlss.family
beta1 = gamlss(
  support ~locus, 
  family = BEINF, 
  data = supp_df
)
### check
summary(beta1)

###
supp_locus = ggplot(
  data = supp_df, 
  aes(x= locus, y= as.numeric(support), fill= locus)
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
    y = "Bootstrap value",
  ) +
  theme(panel.background=element_rect(fill="white"), 
        panel.grid=element_line(colour=NULL),
        panel.border=element_rect(fill=NA,colour="black"),
        axis.title=element_text(size=12,face="bold"),
        axis.text.x=element_text(size=10),
        legend.position = "none")

tiff("figures/boots_partition.tiff",
     units="cm", 
     width= 9.5, 
     height= 8.5,
     res=1200)
  supp_locus
dev.off()

####### AGE EFFECT

### https://www.rdocumentation.org/packages/gamlss.dist/versions/6.1-1/topics/gamlss.family
beta2 = gamlss(
  support ~ age*locus, 
  family = BEINF, 
  data = supp_df
)
### check
summary(beta2)

### predictions
pred2 = ggpredict(beta2, terms = c("age [all]", "locus"))

### plot
supp_age = ggplot() +
  geom_point(
    data = supp_df, 
    aes(x=  as.numeric(age), y= as.numeric(support), color= locus),
    size = 2.5,
    alpha = 0.10
  ) +
  geom_line(
    data = pred2, 
    aes(x = x, y = predicted, color = group), 
    size = 1.2
  ) +
  scale_color_manual(
    values = c(
      "nuclear"= "darkorange",
      "plastidial" = "darkgreen"
    )
  ) +
  labs(
    x = "Node age",
    y = "Bootstrap value",
  ) +
  theme(panel.background=element_rect(fill="white"), 
        panel.grid=element_line(colour=NULL),
        panel.border=element_rect(fill=NA,colour="black"),
        axis.title=element_text(size=12,face="bold"),
        axis.text.x=element_text(size=10),
        legend.position = "none")

tiff("figures/boots_age.tiff",
     units="cm", 
     width= 9.5, 
     height= 8.5,
     res=1200)
  supp_age
dev.off()
