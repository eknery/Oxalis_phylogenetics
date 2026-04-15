### library
if(!require("ips")) install.packages("ips"); library("ips")
if(!require("seqinr")) install.packages("seqinr"); library("seqinr")
if(!require("ape")) install.packages("ape"); library("ape")

### chose input directory
dir_input = "5_concatenated_sequences/"

### partition finder 
partitionfinder(
  alignment = paste0(dir_input, "align.phy"),
  # branchlengths = "linked",
  # models = "raxml",
  # model.selection = "BIC",
  # search = "greedy",
  exec = "C:/Users/eduar/OneDrive/Documents/2_Programas/partitionfinder-2.1.1"
)
