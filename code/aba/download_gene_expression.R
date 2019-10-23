rm(list=setdiff(ls(),c('params','fname')))
basedir <- params$basedir # make this whatever you like, end with /
setwd(basedir)
source('code/aba/aba_fxns.R')
source('code/misc/miscfxns.R')
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

fname <- 'OpioidGeneExpression' # name for output file
fname <- 'AllGeneExpression'

df.expression <- as.data.frame(matrix(ncol=0,nrow=length(region.names.hemi)))
rownames(df.expression) <- region.names.hemi
#opioid.receptors <- grep('Opr',datasets$gene_symbol)
goi.idx <- 1:nrow(datasets) 
for(dataset in goi.idx){ # iterate through every gene and every probe
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

# check all genes were downloaded ... process gene-probe-id naming convention in same way that dataframe naming automatically processes strings
original.geneprobeid <- rep(NA,nrow(datasets))
for(dataset in goi.idx){
  goi <- datasets$gene_symbol[dataset] # get gene name
  if(!is.na(as.numeric(substr(goi,1,1)))){goi <- paste0('X',goi)} # if first character is a number add X
  goi <- gsub('-','.',goi) # replace - with .
  goi <- gsub('\\*','.',goi) # replace * with .
  goi <- gsub('\\(','.',goi) # replace () with .
  goi <- gsub('\\)','.',goi) # replace () with .
  probe <- datasets$probe_name[dataset] # get probe name
  id <- datasets$data_set_id[dataset]
  original.geneprobeid[dataset] <- paste0(goi,'.P.',probe,'.I.',id)
}
# name each processed column name for its original unprocessed name
names(original.geneprobeid) <- sapply(goi.idx, function(dataset) paste0(datasets$gene_symbol[dataset],
                                                                '.P.',datasets$probe_name[dataset],'.I.',datasets$data_set_id[dataset]))
unit.test(identical(unname(original.geneprobeid),names(df.expression)),'all probes downloaded and correctly named','ERROR: probe names not identical')
unit.test(nrow(datasets) == ncol(df.expression),'all probes downloaded','ERROR: missing some probes')

save(df.expression,original.geneprobeid,file = paste0(abadir,'expression/',fname,'.RData'))
write.csv(x = df.expression,file = paste0(abadir,'expression/',fname,'.csv'))

