library(limma)
library(reshape2)
library(tidyr)
library(gtools)
library(ggplot2)
library(utils)
## working directory for analysis 
dire = '/home/b3053674/Documents/LargeStudy/Limma'
baseline_analysis_dir = r'/home/b3053674/Documents/LargeStudy/Limma/Contrasts/BaselineAnalysis'
setwd(dire)

## file containing delta Ct values
data_file = file.path(dire, 'DataFromWaferGen2.csv')

## read the data
data = read.csv(data_file)

## create unique id per sample
data['id'] = paste(data$cell_line, data$treatment, data$time, data$replicate, sep='_')

## remove the Ct column as we only want to work with delta ct values for now
data = data[, !(names(data) %in% c('Ct'))]

## look at the data
head(data)

## isolate pheno data from data 
pd = data[!(names(data) %in% c('Norm2ref', 'Assay' ))]

pd 

## by removing the Norm2ref and Assay columns
pd = pd[!duplicated(pd), ]


## and setting id to index (and removing the id column)
rownames(pd) = pd$id
pd$id = NULL

## Convert replicates and time points to factors
pd$replicate = factor(pd$replicate)
pd$time = factor(pd$time)

## get genes as rows and samples as columns
## First remove pheno data
data$cell_line = NULL
data$treatment = NULL
data$replicate = NULL
data$time = NULL

## reorganize df so id is along the rows and genes the columns
data = spread(data, key='Assay', value='Norm2ref')

## assign id's to rownames and delete as column
rownames(data) = data$id
data$id = NULL

## transpose to get samples along the columns and genes down the rows
data = t(data)

##construct a constrasts matrix for limma
pd$treat_cell_time = factor(paste(pd$treatment, pd$cell_line, pd$time, sep='_'))

design_matrix = model.matrix(~0 + pd$treat_cell_time)
dim(design_matrix)
colnames(design_matrix) = levels(pd$treat_cell_time)
rownames(design_matrix) = pd$treat_cell_time

## view design as csv to validate accuracy
write.csv(design_matrix, file='design_matrix.csv')


fit = lmFit(data, design=design_matrix)

##First set of contrasts I want to make are between the 0 and 96h time 
## point for each cell line. Reassuringly, there are no statistiaclly 
## significant differential expression between these groups


## make contrasts matrix
contrasts_matrix = makeContrasts(A=Baseline_A_0 - Baseline_A_96, 
                                 B=Baseline_B_0 - Baseline_B_96, 
                                 C=Baseline_C_0 - Baseline_C_96, 
                                 D=Baseline_D_0 - Baseline_D_96, 
                                 E=Baseline_E_0 - Baseline_E_96, 
                                 F=Baseline_F_0 - Baseline_F_96, 
                                 G=Baseline_G_0 - Baseline_G_96, 
                                 H=Baseline_H_0 - Baseline_H_96, 
                                 I=Baseline_I_0 - Baseline_I_96, 
                                 levels=design_matrix)

## fit contrasts and do empirical bayes
fit2 = contrasts.fit(fit, contrasts_matrix)
fit2 = eBayes(fit2)
class(fit2)

## look at the results
topTable(fit2, coef=1)
results = decideTests(fit2)
write.csv(results, file='baseline_0_96_contrasts.csv')



contrasts_matrix = makeContrasts(
  'Baseline_A_0 - Baseline_C_0',
  levels=design_matrix
)

fit2 = contrasts.fit(fit, contrasts_matrix)
fit2 = eBayes(fit2)

topTable(fit2, coef=1)
results = decideTests(fit2)
vennDiagram(results)



### Create all combinations of contrast
### for limma  
create_baseline_contrasts = function(){
  l = vector()
  for (cell in levels(pd$cell_line)){
    for (h in c('0', '96')){
      l= c(l, paste0('Baseline_', cell, '_', h))
    }
  }
  
  df = data.frame(combinations(n=length(l), r=2, v=l, repeats.allowed = F))
  df$Contrasts = paste(df$X1 , '-', df$X2)
  df = df[,'Contrasts']
  return (df)
}

baseline_contrasts = create_baseline_contrasts()


## function to remove nan rows
delete.na <- function(DF, n=0) {
  DF[rowSums(is.na(DF)) <= n,]
}


do_contrast = function(contrast, fit, design){
  # print(contrast)
  contrasts_matrix = makeContrasts(contrasts=contrast, levels=design)
  fit2 = contrasts.fit(fit, contrasts_matrix)
  fit2 = eBayes(fit2)

  topTable(fit2, coef=1)
  results = decideTests(fit2)
  # vennDiagram(results)
  
  return (delete.na(results))
}

## Baseline Analysis
setwd(baseline_analysis_dir)
## do all combinations of contrast for baseline
bc = do_contrast(baseline_contrasts, fit, design_matrix)

## summarize  
sample_sums = sort(colSums(bc))
gene_sums = sort(rowSums(bc))

sample_sums = sample_sums[!sample_sums == 0]
gene_sums = gene_sums[!gene_sums == 0]

length(sample_sums)
length(gene_sums)

write.csv(bc, 'baseline_contrasts.csv')
write.csv(sort(sample_sums), 'sum_of_sample_contrasts.csv')
write.csv(sort(gene_sums), 'sum_of_gene_contrasts.csv')


gene_sums
gene_sums = data.frame(gene_sums)
sample_sums = data.frame(sample_sums)


## plot summary of numbers of contrasts each gene is differentially regualted in 
barplot(sort(abs(gene_sums)/153*100, decreasing=TRUE), names.arg=rownames(gene_sums), las=2, col='black',
        main='Barchart summarizing number of baseline contrasts each gene is \ndifferentially regulated in',
        ylab='Percentage of Contrasts (n=153)',
        ylim=c(0, 100),
        cex.axis=1.4,
        )       


sample_sums

## plot summary of numbers of contrasts each gene is differentially regualted in 

barplot(sort(abs(sample_sums$sample_sums), decreasing=TRUE), col='black',
        main='Barchart summarizing number genes differentially regulated \nin each contrast',
        ylab='Percentage of Genes in Contrast (n_genes=66)',
        cex.axis=1.4,
)       


bc_sort = bc[,rownames(sample_sums)]

dim(bc_sort)

bc_sort
vennDiagram(bc)
bc_sort[,1:4]
vennDiagram(bc_sort[,1:3], names=c('C_0 - H_96', 'D_0 - H_96', 'B_0 - H_96', 'D_0 - I_96'))


write.csv(bc_sort, file='baseline_contrasts_sorted.csv')

bc_sort[,1:3]#][!bc_sort[,1:3] == 0]


get_common_genes = function(...){
  args = list(...)
  #Reduce(function(x, y) paste(x, '&' y), args)
  print (args)
}

common = bc_sort[, 'Baseline_H_96 - Baseline_I_96'] & bc_sort[, 'Baseline_H_96 - Baseline_I_0']

topTable(fit2[common,])

##############################################################
##                                                          ##
##            Questions                                     ##
##############################################################
##
## 1) Which genes are differentially regulated between neonatal, adult and senescent cells?  
##     
## 2) Are any of these genes per cell line DE as a result of time only?
##
## 3) Are any of these genes per cell line DE as a result of media change?
##
## 4) Are any of these genes per cell line DE as a result of adding TGFb?
##
## after Q3 and Q4 are answered then dig into the time element i.e. at which times are the genes DE? 
##
##
## 5) are there are genes in the second WG experiment that were also DE in the first?
##
## 6) Which genes are differently regulated
##
##
##


## 1) Which genes are differentially regulated between neonatal, adult and senescent cells?  




















































