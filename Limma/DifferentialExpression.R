library(limma)
library(reshape2)
library(tidyr)
## working directory for analysis 
dire = '/home/b3053674/Documents/LargeStudy/Limma'
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


## First I'll analyse the baseline expression levels. This is simpler than 
## analysing the whole experiment and will help with getting used to limma again
head(pd)

time_vec= unique(pd$time)
replicates = unique(pd$replicate)


# treat = c('Baseline', 'Control', 'TGFb')
# cell_line = c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I')
# time = c('t0', 't0.5', 't1', 't2', 't3', 't4', 't8', 't12', 't24', 't48', 't72', 't96')

head(pd)

pd$cell_treat_time = factor(paste(pd$cell_line, pd$treatment, pd$time, sep='_'))

design_matrix = model.matrix(~0 + pd$cell_treat_time)

colnames(design_matrix) = pd$cell_treat_time


# 
# ### view design as csv to validate that it looks reasonable 
write.csv(design_matrix, file='design_matrix.csv')
# 
# ## A note on design matrices. Some variables look like they are missing. 
# ## This is because the missing column is the default (i.e. when all 
# ## columns are 0 the design take on the \'missing\' value by default)
# 
# 
# design_matrix
# 
# fit = lmFit(data, design=design_matrix)
# 
# contrasts_matrix = makeContrasts(Baseline+C, levels=design_matrix)
# 
# 
# contrasts_matrix
# 
# fit2 = contrasts.fit(fit, contrasts_matrix)
# fit2 = eBayes(fit2)
# topTable(fit2, coef=1, adjust='BH')
# 
# 
# results = decideTests(fit2)
# 
# vennDiagram(results)
# 
# 
# 






















