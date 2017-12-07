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
## Same as in the other scripts

## working directory for analysis 
dire = '/home/b3053674/Documents/LargeStudy/Limma'
regression_dire = '/home/b3053674/Documents/LargeStudy/RegressionAnalysis'
setwd(regression_dire)

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
data = data.frame(data)
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

## create function for calculating r2
## df = dataframe. genes on columns samples in rows
## x = predictor variable
## y = response variable
reg = function(df, x, y){
  reg_fun = lm(formula=df[[y]] ~ df[[x]])
  slope = round(coef(reg_fun)[2], 3)
  r2 = round(as.numeric(summary(reg_fun)[8]), 3)
  r2.adj = round(as.numeric(summary(reg_fun)[9]), 3)
  intercept = round(coef(reg_fun)[1], 3)
  c(slope, intercept, r2, r2.adj)
  
}

## test reg
reg(data, 'CTGF', 'COL1A1')

## Function to plot correlations between x and y using data. 

corr_plot = function(data, x, y, eq_loc=c(5, 20)){
  
  ##isolate CTGF and COL1A1
  data = data[, c(x, y)]
  
  ## split id into separate columns in dataframe
  split = data.frame(str_split_fixed(rownames(data),'_', 4))

  ## row and col names
  rownames(split) = rownames(data)
  head(split)
  colnames(split) = c('cell_line', 'treatment', 'time', 'replicate')

  ## merge split frame into col_vs_ctgf_data
  data = merge(data, split, by=0, all=T)
  ## calculate regression data and relabel cols
  
  print (data)
  reg_data = ddply(data, c("treatment", 'cell_line'), reg, x, y)
  colnames(reg_data) = c('treatment', 'cell_line', 'slope', 'intercept', 'r2', 'r2.adj')

  # ## plot data
  ggplot(data,
         aes_string(x=x, y=y, group='cell_line', color='treatment')) +
    geom_point()  +
    ## convert to facet grid to plot cell lines against treatments
    facet_grid(cell_line ~ treatment, scales = 'fixed') +

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
    ggtitle(paste('Linear Regression. x=', x, 'y=', y)) +

    ## add equations for labels
    geom_label(data=reg_data,
               inherit.aes = F,
               aes(x=eq_loc[1], y=eq_loc[2],
                   label=paste('y=',
                               intercept,"+",
                               slope, "x, r2=", r2),
                   size=6)) +
    labs(x=paste(x, '(DeltaCt)'), y=paste(y, '(DeltaCt)'))
}


## Genes available
# [1] "ACTA2"    "ADAMTS1"  "ATP6AP1"  "B2M"      "BGN"      "BHLHE40"  "CAV1"     "CDKN2A"   "CDKN2B"  
# [10] "COL1A1"   "COL1A2"   "COL2A1"   "COL4A1"   "COL5A1"   "CTGF"     "DCN"      "EGR2"     "ELN"     
# [19] "ENG"      "ETS1"     "FBLN1"    "FBN1"     "FGF2"     "FN1"      "FOSB"     "GADD45B"  "GUSB"    
# [28] "HAS2"     "ID1"      "IL1A"     "IL1B"     "IL6"      "ITGA1"    "ITGA2"    "JUN"      "JUNB"    
# [37] "LARP6"    "LOXL1"    "LOXL2"    "LTBP2"    "MMP1"     "MMP13"    "MMP14"    "MMP2"     "NOX4"    
# [46] "PDGFA"    "PMEPA1"   "PPIA"     "PPP3CA"   "PSMD14"   "PTEN"     "RARA"     "RARG"     "RHOB"    
# [55] "SERPINE1" "SERPINE2" "SKI"      "SKIL"     "SMAD3"    "SMAD7"    "SPARC"    "TGFB1"    "TGFBR1"  
# [64] "TGFBR2"   "THBS1"    "THBS2"    "TIMP1"    "TIMP3"    "TP53BP1"  "VCAN"     "VEGFA"    "VIM" 


corr_plot(data, 'ELN', 'COL1A1', eq_loc=c(1, 40))






