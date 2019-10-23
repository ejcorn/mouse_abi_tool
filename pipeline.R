rm(list=ls())
basedir <- '~/Dropbox/Cornblath_Bassett_Projects/MouseGeneExpressionABA/mouse_abi_tool/'
setwd(basedir)
params <- list(basedir=basedir)
source('code/misc/miscfxns.R')
params$source.save <- source.save
params$opdir <- paste('ABA_expression',params$c.max,'/',sep='')
dir.create(params$opdir,recursive = T)

#################################################
### Load packages & create output directories ###
#################################################

source('code/misc/packages.R')
source('code/misc/directories.R')

############################
### Load connectome data ###
############################

source('code/process/process.R')

#########################################
### Save all ABI gene expression data ###
#########################################

source('code/aba/process_ontology.R')
fname <- 'AllGeneExpression'
source('code/aba/download_gene_expression.R')
source('code/aba/probeselection_normalization.R')
source('code/aba/probenormalization_selection.R')
