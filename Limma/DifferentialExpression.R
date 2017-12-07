library(limma)
library(reshape2)
library(tidyr)
library(gtools)
library(ggplot2)
library(utils)
library(stringi)
library(stringr)
library(Rmisc)
################################################################
##                                                            ##   
##            Data Preprocessing                              ##
##                                                            ##
################################################################

## working directory for analysis 
dire = '/home/b3053674/Documents/LargeStudy/Limma'
baseline_analysis_dir = '/home/b3053674/Documents/LargeStudy/Limma/Contrasts/BaselineAnalysis'
setwd(dire)

## file containing delta Ct values
data_file = file.path(dire, 'DataFromWaferGen2.csv')

## read the data
data = read.csv(data_file)

## Remove Cell line F, time 8h, rep 6, TGFb treated as an anomoly (consult PCA)
data = data[!(data$cell_line == 'F' & data$treatment == 'TGFb' & data$replicate == 6 & data$time == 8),]

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

################################################################
##                                                            ##   
##        Plot COL1A1 against CTGF. Correlated?               ##
##                                                            ##
################################################################

##isolate CTGF and COL1A1
col_vs_ctgf_data = t(data[c('CTGF', 'COL1A1'),])

## Look at data
head(col_vs_ctgf_data)

## split id into separate columns in dataframe
split = data.frame(str_split_fixed(rownames(col_vs_ctgf_data),'_', 4))

## row and col names
rownames(split) = rownames(col_vs_ctgf_data)
head(split)
colnames(split) = c('cell_line', 'treatment', 'time', 'replicate')

## merge split frame into col_vs_ctgf_data
col_vs_ctgf_data = merge(col_vs_ctgf_data, split, by=0, all=T) 

dim(col_vs_ctgf_data)
head(col_vs_ctgf_data)


## create function for calculating r2
reg = function(df){
  reg_fun = lm(formula=df$COL1A1 ~ df$CTGF)
  slope = round(coef(reg_fun)[2], 3)
  r2 = round(as.numeric(summary(reg_fun)[8]), 3)
  r2.adj = round(as.numeric(summary(reg_fun)[9]), 3)
  intercept = round(coef(reg_fun)[1], 3)
  c(slope, intercept, r2, r2.adj)
}

## calculate regression data and relabel cols
reg_data = ddply(col_vs_ctgf_data, c("treatment", 'cell_line'), reg)
reg_data
colnames(reg_data) = c('treatment', 'cell_line', 'slope', 'intercept', 'r2', 'r2.adj')

## plot data 
ggplot(col_vs_ctgf_data, 
       aes(x=CTGF, y=COL1A1, group=cell_line, color=treatment)) +
  geom_point() +
  
  ## use black and white theme as base
  theme_bw() + 
  theme(strip.text.x = element_text(size=16), 
        strip.text.y = element_text(size=16), 
        plot.title = element_text(hjust = 0.5, size=22),
        axis.text = element_text(size=10),
        axis.text.x = element_text(size=14), 
        axis.text.y = element_text(size=14)) + 

  ## add regression line with standard error
  geom_smooth(method=lm, se=T) +
  
  ## add title
  ggtitle('Linear Regression. x=CTGF, y=COL1A1') +
  
  ## convert to facet grid to plot cell lines against treatments
  facet_grid(cell_line ~ treatment, scales = 'fixed') + 
  
  ## add equations for labels
  geom_label(data=reg_data, 
             inherit.aes = F, 
             aes(x=4, y=45, 
                 label=paste('y=', 
                             intercept,"+", 
                             slope, "x, r2=", r2))) +
  labs(x='CTGF (DeltaCt)', y='COL1A1 (DeltaCt)') 
  # guide_legend(name='Treatment')



################################################################
##                                                            ##   
##            Based on correlation of CTGF and COL1A1         ##
##            make general function to do same with any gene  ##
##            pair                                            ##
################################################################
##isolate CTGF and COL1A1

data = t(data)
data = data.frame(data)


## create function for calculating r2
reg = function(df){
  reg_fun = lm(formula=df$COL1A1 ~ df$CTGF)
  slope = round(coef(reg_fun)[2], 3)
  r2 = round(as.numeric(summary(reg_fun)[8]), 3)
  r2.adj = round(as.numeric(summary(reg_fun)[9]), 3)
  intercept = round(coef(reg_fun)[1], 3)
  c(slope, intercept, r2, r2.adj)
  
}


## create function for calculating r2
reg = function(df, x, y){
  reg_fun = lm(formula=df[[y]] ~ df[[x]])
  slope = round(coef(reg_fun)[2], 3)
  r2 = round(as.numeric(summary(reg_fun)[8]), 3)
  r2.adj = round(as.numeric(summary(reg_fun)[9]), 3)
  intercept = round(coef(reg_fun)[1], 3)
  c(slope, intercept, r2, r2.adj)
  
}


reg(data, 'CTGF', 'COL1A1')

class(data)



head(data)

corr_plot = function(data, x, y){
  ##isolate CTGF and COL1A1
  data = t(data[c(x, y),])
  
  ## Look at data
  head(data)

  ## split id into separate columns in dataframe
  split = data.frame(str_split_fixed(rownames(data),'_', 4))

  ## row and col names
  rownames(split) = rownames(data)
  head(split)
  colnames(split) = c('cell_line', 'treatment', 'time', 'replicate')

  ## merge split frame into col_vs_ctgf_data
  data = merge(data, split, by=0, all=T)

  ## calculate regression data and relabel cols
  reg_data = ddply(data, c("treatment", 'cell_line'), reg)
  reg_data
  # colnames(reg_data) = c('treatment', 'cell_line', 'slope', 'intercept', 'r2', 'r2.adj')
  # 
  # ## plot data 
  # ggplot(data, 
  #        aes(x=x, y=y, group=cell_line, color=treatment)) +
  #   geom_point() +
  #   
  #   ## use black and white theme as base
  #   theme_bw() + 
  #   theme(strip.text.x = element_text(size=16), 
  #         strip.text.y = element_text(size=16), 
  #         plot.title = element_text(hjust = 0.5, size=22),
  #         axis.text = element_text(size=10),
  #         axis.text.x = element_text(size=14), 
  #         axis.text.y = element_text(size=14)) + 
  #   
  #   ## add regression line with standard error
  #   geom_smooth(method=lm, se=T) +
  #   
  #   ## add title
  #   ggtitle(paste('Linear Regression. x=', x, 'y=', y)) +
  #   
  #   ## convert to facet grid to plot cell lines against treatments
  #   facet_grid(cell_line ~ treatment, scales = 'fixed') + 
  #   
  #   ## add equations for labels
  #   geom_label(data=reg_data, 
  #              inherit.aes = F, 
  #              aes(x=4, y=45, 
  #                  label=paste('y=', 
  #                              intercept,"+", 
  #                              slope, "x, r2=", r2))) +
  #   labs(x='CTGF (DeltaCt)', y='COL1A1 (DeltaCt)') 
}

reg(data, 'COL1A1', 'COL1A2')

################################################################
##                                                            ##   
##            LIMMA Stuff                                     ##
##                                                            ##
################################################################




##construct a constrasts matrix for limma
pd$treat_cell_time = factor(paste(pd$treatment, pd$cell_line, pd$time, sep='_'))

design_matrix = model.matrix(~0 + pd$treat_cell_time)



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


contrasts_matrix = makeContrasts(A=Baseline_A_96 - Baseline_A_0, 
                                 B=Baseline_B_96 - Baseline_B_0, 
                                 C=Baseline_C_96 - Baseline_C_0, 
                                 D=Baseline_D_96 - Baseline_D_0, 
                                 E=Baseline_E_96 - Baseline_E_0, 
                                 F=Baseline_F_96 - Baseline_F_0, 
                                 G=Baseline_G_96 - Baseline_G_0, 
                                 H=Baseline_H_96 - Baseline_H_0, 
                                 I=Baseline_I_96 - Baseline_I_0, 
                                 levels=design_matrix) 
## fit contrasts and do empirical bayes
fit2 = contrasts.fit(fit, contrasts_matrix)
fit2 = eBayes(fit2)
class(fit2)

## look at the results
topTable(fit2, coef=1, sort.by = 'P', adjust.method = 'holm')
results = decideTests(fit2)
write.csv(results, file='baseline_0_96_contrasts.csv')



contrasts_matrix = makeContrasts(
  'Baseline_A_0 - Baseline_D_0',
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






## 2) Are any of these genes per cell line DE as a result of time only?


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
''' FN1      0.8352343 0.4555378 0.70056509 0.7066410 0.5994690'''


contrasts_matrix = makeContrasts(A=Baseline_A_0 - Baseline_A_96, 
                                 levels=design_matrix)
                                 
                                 
## fit contrasts and do empirical bayes
fit2 = contrasts.fit(fit, contrasts_matrix)
fit2 = eBayes(fit2)
class(fit2)
fit2$p.value

## look at the results
topTable(fit2)
results = decideTests(fit2)

results
write.csv(results, file='baseline_0_96_contrasts.csv')


## trying new design matrix
## Sinc baseline is kind of separate from the rest of the data, shold I only include the baseline 
## data in the statistical model? 

## I think it is worth extracting the baseline data prior to fitting the first model
## Analyze the Baseline condition in isolation from the control and treated
design_matrix_interaction = model.matrix(~0 + pd$cell_line*pd$treatment*pd$time)


write.csv(design_matrix_interaction, file='design_matrix_interaction.csv')
dim(design_matrix_interaction)









































