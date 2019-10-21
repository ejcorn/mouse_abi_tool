rm(list=setdiff(ls(),c('params')))
basedir <- params$basedir # make this whatever you like, end with /
setwd(basedir)
source('code/aba/aba_fxns.R')
abadir <- 'data/aba/'
dir.create(paste0(abadir,'expression'),recursive = T)
load(file=paste0(abadir,'ontologies/keys.RData'))
load(paste(params$opdir,'processed/connectome.RData',sep=''))  # load path data and ROI names

dfile <- paste0(basedir,'data/aba/mouse_expression_data_sets.csv')
if(!file.exists(dfile)){
  download.file(url = 'http://download.alleninstitute.org/informatics-archive/october-2014/mouse_expression/mouse_expression_data_sets.csv',
  destfile = dfile)}
datasets <- read.csv('data/aba/mouse_expression_data_sets.csv',stringsAsFactors = F)
# see http://download.alleninstitute.org/informatics-archive/october-2014/mouse_expression/Accessing_October_2014_expression_data.pdf
# obtained from http://download.alleninstitute.org/informatics-archive/october-2014/mouse_expression/mouse_expression_data_sets.csv
# https://www.pnas.org/content/113/5/1435 says it's okay to use both sagittal and coronal sections

df.expression <- as.data.frame(matrix(ncol=0,nrow=length(region.names.hemi)))
rownames(df.expression) <- region.names.hemi
opioid.receptors <- grep('Opr',datasets$gene_symbol)
goi.idx <- 1:nrow(datasets)
for(dataset in opioid.receptors){ # iterate through every gene and every probe
  goi <- datasets$gene_symbol[dataset] # get gene name
  probe <- datasets$probe_name[dataset] # get probe name
  id <- datasets$data_set_id[dataset]
  url <- datasets$structure_unionizes_file_url[dataset] # get url
  dfile <- paste0(abadir,'expression/',goi,probe,'.csv') # download file
  download.file(url = url,destfile = dfile)
  gene.exp <- read.csv(dfile)
  gene.exp <- gene.exp[,c('structure_id','expression_energy')] # remove extraneous info
  #write.csv(x = gene.exp[gene.exp$structure_id %in% unlist(region.names.id.key),], file = dfile,row.names = F) # rewrite full gene expression file
  df.expression[,paste0(goi,'.P.',probe,'.I.',id)] <- ontology.assign.expression(region.names.id.key,gene.exp) # default is expression energy
  file.remove(dfile)
}

write.csv(x = df.expression,file = paste0(abadir,'expression/OpioidGeneExpression.csv'))
