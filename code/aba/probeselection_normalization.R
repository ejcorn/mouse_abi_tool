rm(list=setdiff(ls(),c('params')))
basedir <- params$basedir # make this whatever you like, end with /
setwd(basedir)
source('code/aba/aba_fxns.R')
abadir <- 'data/aba/'
dir.create(paste0(abadir,'expression'),recursive = T)
load(file=paste0(abadir,'ontologies/keys.RData'))
load(paste(params$opdir,'processed/connectome.RData',sep=''))  # load path data and ROI names
fname <- 'OpioidGeneExpression'

datasets <- read.csv('data/aba/mouse_expression_data_sets.csv',stringsAsFactors = F)
df.expression.probes <- read.csv(file = paste0(abadir,'expression/',fname,'.csv'))

#######################
### Probe selection ###
#######################

# see https://www.sciencedirect.com/science/article/pii/S1053811919300114?via%3Dihub
# mean across probes is most consistent when benchmarked against RNA seq data, though other methods are possible

df.expression <- as.data.frame(matrix(ncol=0,nrow=length(region.names.hemi)))
rownames(df.expression) <- region.names.hemi
genes <- unique(datasets$gene_symbol) # select gene names
genes <- unique(datasets$gene_symbol[grep('Opr',datasets$gene_symbol)]) # select gene names
for(gene in genes){
  gene.mask <- grep(gene,names(df.expression.probes)) # find all probes/datasets for a given gene
  df.expression[,gene] <- rowMeans(df.expression.probes[,gene.mask],na.rm = T) # store the mean across those datasets, ignoring NaNs
}

#################################
### Normalize gene expression ###
#################################

sigmoid.norm <- function(x){
  # https://www.pnas.org/content/113/5/1435
  # robust to outliers, puts each gene on same scale relative to its own expression across brain
  return((1+exp(-scale(x,center = T)))^-1)
}

for(gene in genes){
  df.expression[,gene] <- sigmoid.norm(df.expression[,gene])
}

write.csv(x=df.expression,file = paste0(abadir,'expression/',fname,'_Normalized.csv'))