rm(list=setdiff(ls(),c('params','fname')))
basedir <- params$basedir # make this whatever you like, end with /
setwd(basedir)
source('code/aba/aba_fxns.R')
abadir <- 'data/aba/'
dir.create(paste0(abadir,'expression'),recursive = T)
load(file=paste0(abadir,'ontologies/keys.RData'))
load(paste(params$opdir,'processed/connectome.RData',sep=''))  # load path data and ROI names

datasets <- read.csv('data/aba/mouse_expression_data_sets.csv',stringsAsFactors = F)
load(file = paste0(abadir,'expression/',fname,'.RData'))
df.expression.probes <- df.expression 

#################################
### Normalize gene expression ###
#################################

# in this file, I'm normalizing before averaging across probes because it intuitively seems better

sigmoid.norm <- function(x){
  # https://www.pnas.org/content/113/5/1435
  # robust to outliers, puts each gene on same scale relative to its own expression across brain
  return((1+exp(-scale(x,center = T)))^-1)
}

for(probe in names(df.expression.probes)){
  df.expression.probes[,probe] <- sigmoid.norm(df.expression.probes[,probe])
}

#######################
### Probe selection ###
#######################

# see https://www.sciencedirect.com/science/article/pii/S1053811919300114?via%3Dihub
# mean across probes is most consistent when benchmarked against RNA seq data, though other methods are possible

df.expression <- as.data.frame(matrix(ncol=0,nrow=length(region.names.hemi)))
rownames(df.expression) <- region.names.hemi
genes <- unique(datasets$gene_symbol) # select gene names
for(gene in genes){
  geneprobeid.names <- names(original.geneprobeid)[grep(gene,names(original.geneprobeid))] # find all probes/datasets for a given gene based on original names
  geneprobeid.names <- unname(original.geneprobeid[geneprobeid.names]) # get dataframe column names
  df.expression[,gene] <- rowMeans(df.expression.probes[,geneprobeid.names,drop = FALSE],na.rm = T) # store the mean across those datasets, ignoring NaNs
}

write.csv(x=df.expression,file = paste0(abadir,'expression/',fname,'_NormalizedFirst.csv'))
