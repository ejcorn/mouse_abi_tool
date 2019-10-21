rm(list=setdiff(ls(),c('params')))
basedir <- params$basedir # make this whatever you like, end with /
setwd(basedir)
source('code/aba/aba_fxns.R')
source('code/misc/miscfxns.R')
abadir <- paste0(basedir,'data/aba/')
dir.create(paste0(abadir,'ontologies'),recursive = T)

# now we need to get a list of relevant structure ids to retrieve the expression data for the regions we measured in our experiments
# this link has info on the ABA hierarchical structural ontology:
# http://help.brain-map.org/display/api/Atlas+Drawings+and+Ontologies#AtlasDrawingsandOntologies-StructuresAndOntologies
# http://help.brain-map.org/display/api/Downloading+an+Ontology%27s+Structure+Graph
# the Ontology id of adult mouse is 1, so you plug in this link to get data
# http://api.brain-map.org/api/v2/structure_graph_download/1.json
# need to download R json package to work with this data

ontology.url <- 'http://api.brain-map.org/api/v2/structure_graph_download/1.json' # ontology of structures for adult mouse
ontology.destfile <- paste0(abadir,'ontologies/adult_mouse.json')
download.file(url = ontology.url, destfile = ontology.destfile) # save ontology file

ontology <- fromJSON(file = ontology.destfile)$msg[[1]] # remove initial fluff, file starts at $msg[[1]]

# make a simplified ontology that only contains acronym and id
unique.identifier <- 'ElIcOrNbLaTh' # use any relatively unique string just to make sure you extract the right list elements in next step
ontology <- ontology.name.children(ontology)
ontology <- ontology.remove.fluff(ontology,unique.identifier)

# now make a data frame that links every single structure name with its id, that corresponds to expression files
ont.un <- unlist(ontology) # unlist ontology
# extract ids and acronyms as a vector, searching for unique identifier to ensure specificity
structure_ids <- ont.un[grep(paste0(".id",unique.identifier),names(ont.un))]
mode(structure_ids) <- 'numeric' # convert to numeric while preserving names
acronyms <- ont.un[grep(paste0(".acronym",unique.identifier),names(ont.un))]
# remove .id* or .acronym* and make sure the stem (i.e. grey.children.CTX.children.CH. etc... matches up)
# this ensures that the id definitely matches the acronym
names(structure_ids) <- str_replace(string = names(structure_ids), pattern =paste0('.id',unique.identifier),replacement = '')
names(acronyms) <- str_replace(string = names(acronyms), pattern =paste0('.acronym',unique.identifier),replacement = '')
# now the names of these two vectors should be identical if they match
unit.test(identical(names(acronyms),names(structure_ids)),'ids and structure names matched correctly','ERROR: ids and structure names matched incorrectly')
# cool, now we can join them into the same data frame
df.key <- data.frame(structure_id = structure_ids, 
                 acronym = acronyms,stringsAsFactors = F)

# next we can look through the region names from our experiment and get the corresponding structure ids, for ABA and CNDR names
load(paste(params$opdir,'processed/connectome.RData',sep=''))  # load path data and ROI names

region.names.id.key <- ontology.get.name.id.key(region.names.hemi,df.key)
# check for regions that weren't able to be labeled
unlist(region.names.id.key)[is.na(unlist(region.names.id.key))]

save(region.names.id.key,file=paste0(abadir,'ontologies/keys.RData'))
