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

gene.probe.idx <- names(df.expression.probes)
probe.names <- gsub('\\.I\\..*','',gsub('.*\\.P\\.','',gene.probe.idx))
identical(probe.names,datasets$probe_name)

#######################################################################################
### per Fulcher et al. 2019, discard genes to enrich for brain genes w quality data ### 
#######################################################################################
# https://www.pnas.org/content/pnas/suppl/2019/02/18/1814144116.DCSupplemental/pnas.1814144116.sapp.pdf
# they retain genes that are either:
# coronal OR
# r >= 0.5 between multiple probes

# 09/17/20: in the future can use this data: https://figshare.com/articles/dataset/Processed_voxel-level_gene_expression_data_processed_from_the_Allen_Developing_Mouse_Brain_Atlas/12763604
# tabula muris consortium: cell composition for entire mouse but not too many brain areas
# allen has single cell data


coronal.mask <- datasets$plane_of_section == 'coronal'
coronal.genes <- datasets$gene_symbol[coronal.mask]
sagittal.genes <- datasets$gene_symbol[!coronal.mask] # at this point only retain sagittal sections with more than one probe
sagittal.genes <- unique(sagittal.genes[duplicated(sagittal.genes)])
genes.retained <- unique(c(coronal.genes,sagittal.genes))

#########################################
### Probe selection and normalization ###
#########################################

norm.fun <- function(x){
  # https://www.pnas.org/content/113/5/1435
  # robust to outliers, puts each gene on same scale relative to its own expression across brain
  return((1+exp(-scale(x,center = T)))^-1)
}

# see https://www.sciencedirect.com/science/article/pii/S1053811919300114?via%3Dihub
# mean across probes is most consistent when benchmarked against RNA seq data, though other methods are possible

agreement.thresh <- 0.5
df.expression <- as.data.frame(matrix(ncol=0,nrow=length(region.names.hemi)))
rownames(df.expression) <- region.names.hemi
for(gene in genes.retained){
  print(which(genes.retained %in% gene))
  geneprobeid.names <- names(original.geneprobeid)[grep(gene,names(original.geneprobeid),fixed = T)] # find all probes/datasets for a given gene based on original names
  geneprobeid.names <- unname(original.geneprobeid[geneprobeid.names]) # get dataframe column names
  G <- df.expression.probes[,geneprobeid.names,drop = FALSE] # get gene expression in matrix
  if(gene %in% coronal.genes){
    df.expression[,gene] <- rowMeans(apply(G,2,norm.fun),na.rm = T) # store the mean across those datasets, ignoring NaNs
  } else{
    r.mat <- cor(G,use='pairwise.complete.obs') # correlation between multiple section datasets
    r.mean <- mean(r.mat[!diag(TRUE,length(geneprobeid.names))]) # agreement between datasets averaged
    if(r.mean >= agreement.thresh){
      df.expression[,gene] <- rowMeans(apply(G,2,norm.fun),na.rm = T) # store the mean across those datasets, ignoring NaNs
    }
  }
}

write.csv(x=df.expression,file = paste0(abadir,'expression/',fname,'_Fulcher2019.csv'))
