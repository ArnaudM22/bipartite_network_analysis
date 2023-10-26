
# Setting the working directory.
setwd('C:/Users/arnau/Desktop/bipartite_network_analysis/Code')

# Importing the libraries.
library(devtools)
library(RVAideMemoire)
library(netZooR)
library(dplyr)
library(data.table)

# Setting the seed.
set.seed(65)
# Importing the edgelist
elist = read.table('../datasets/phase_2/tcx_lostAD_link.csv', sep = ',', header = TRUE)
elist2 <- elist[, c('teName','geneName', 'coef')]
elist2[,'coef'] = abs(elist2[,'coef'])
# Constructing the condor object
condor.object <- createCondorObject(elist2)
# Runing the bipartite Community detection
condor.object <- condorCluster(condor.object, cs.method = "LEC", project =F)


# Family enrichment analysis
#Defining the family enrichment function : 
#return as a list Fisher exact score, FIsher pairwise comparison and count table.
tefamily_comu_test <- function(fam_of_interest, count_table){
  colname = colnames(count_table)
  colname = colname[colname != fam_of_interest]
  other_family_sum = as.list(rowSums(count_table[,colname]))
  count_table = as.data.frame.matrix(cbind(count_table[,fam_of_interest], other_family_sum))
  colnames(count_table)[1]= fam_of_interest
  count_table[,1] = unlist(count_table[,1])
  count_table[,2] = unlist(count_table[,2])
  mosaicplot(count_table,
             main = "Mosaic plot",
             color = TRUE
  )
  exp = chisq.test(count_table)$expected
  results = fisher.test(count_table, workspace = 2e9, simulate.p.value = TRUE, B = 1000000)
  #pairwise analysis :
  count_table_for_fish = as.matrix(count_table)
  multicomp = fisher.multcomp(count_table_for_fish)
  return(as.list(results, multicomp, count_table))
}
tlist =  condor.object$red.memb
names(tlist)[1] ="teName"
famlist = elist[c('teName', 'repFamily')]
tlist = merge(tlist, famlist, by = "teName")
count_table = as.data.frame.matrix(table(tlist[c('com', 'repFamily')]))
#Computing the family enrichment
result_list_L1 = tefamily_comu_test(fam_of_interest = 'L1', count_table = count_table)
result_list_Alu = tefamily_comu_test(fam_of_interest ='Alu' , count_table = count_table)
result_list_SVA = tefamily_comu_test(fam_of_interest ='SVAs' , count_table = count_table)

# TE-ZNF age independance analysis
age_table = elist[,c('kznf_age', 'te_age')]
age_table = as.data.frame.matrix(table(age_table))
exp = chisq.test(age_table)$expected
fisher.test(age_table, workspace = 2e9, simulate.p.value = TRUE, B = 1000000)
chisq.test(age_table)

# ZNF age enrichment in modules analysis.
zlist = condor.object$blue.memb
names(zlist)[1] ="geneName"
agelist = elist[c('geneName', 'kznf_age')]
zlist = merge(zlist, agelist, by = "geneName")
count_table = as.data.frame.matrix(table(zlist[c('com', 'kznf_age')]))
mosaicplot(count_table,
           main = "Mosaic plot",
           color = TRUE
)
exp = chisq.test(count_table)$expected
fisher.test(count_table, workspace = 2e9, simulate.p.value = TRUE, B = 1000000)
count_table_for_multicomp = as.matrix(count_table)
multi_comp_table = fisher.multcomp(count_table_for_multicomp)
multi_comp_table = multi_comp_table$p.value
write.csv(multi_comp_table,file="../results/age_enrichment_munticount.csv")

 
