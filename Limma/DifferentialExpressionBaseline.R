library(limma)
library(reshape2)
library(tidyr)
library(gtools)
library(ggplot2)
library(utils)
## working directory for analysis 
dire = '/home/b3053674/Documents/LargeStudy/Limma'
baseline_analysis_dir = '/home/b3053674/Documents/LargeStudy/Limma/Contrasts/BaselineAnalysis'
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


pd_baseline = pd[pd$treatment == 'Baseline',]
baseline_ids = rownames(pd_baseline)
baseline_ids

baseline_data = data[,baseline_ids]



































