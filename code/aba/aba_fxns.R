ontology.name.children <- function(ontology){ 
  # INPUTS:
  # ontology: nested list from ABA structure ontology
  
  # OUTPUTS:
  # ontology: input ontology list, with children recursively renamed for their acronyms
  
  names(ontology$children) <- sapply(ontology$children, function(x) x$acronym)
  ontology$children <- lapply(ontology$children, function(x) ontology.name.children(x))
  
  return(ontology)
}

ontology.remove.fluff <- function(ontology,unique.identifier = 'ElIcOrNbLaTh'){
  # INPUTS:
  # ontology: nested list from ABA structure ontology
  
  # OUTPUTS:
  # ontology: input ontology list, with all fields removed except for id, children and acronym
  
  ontology[!names(ontology) %in% c('acronym','id','children')] <- NULL
  names(ontology)[names(ontology) == 'id'] <- paste0('id',unique.identifier) # give it a unique identifier so no chance of ambiguity later
  names(ontology)[names(ontology) == 'acronym'] <- paste0('acronym',unique.identifier) # give it a unique identifier so no chance of ambiguity later
  ontology$children <- lapply(ontology$children, function(x) ontology.remove.fluff(x))
  return(ontology)
}

ontology.get.name.id.key <- function(name.vec,df.key){
  # INPUTS
  # name.vec: vector of region names
  # df.key: data frame with structure_id column containing numeric ABA structure_ids and acronym column containing character structure name
  # OUTPUTS:
  # names.id.key: list whose element names are the names in name.vec, and whose elements contain corresponding structure ids to that name
  
  names.id.key <- as.list(rep(NA,length(name.vec))) # default id is NA
  names(names.id.key) <- name.vec
  for(i.name in name.vec){ # get structure_id corresponding to each name and store in the corresponding list element
    i.name.id <- df.key$structure_id[df.key$acronym == i.name]
    if(length(i.name.id) > 0){names.id.key[i.name] <- i.name.id} # leave as NA if name doesn't exist in acronym
  }
  return(names.id.key)
}

ontology.assign.expression <- function(names.id.key,goi.exp,expression.metric = 'expression_energy'){
  # INPUTS:
  # names.id.key: list whose element names are the names in name.vec, and whose elements contain corresponding structure ids to that name
  # goi.exp: dataframe containing expression values for each structure id
  # expression.metric: character corresponding to column of goi.exp that names a metric of expression to extract
  # ^ i.e. energy, intensity, density
  
  # OUTPUTS:
  # expression: named vector of numeric gene expression values for given expression metric and each region in names.id.key
  # NA means names.id.key had no match to a structure_id to begin with
  # or there was no match for structure_id in expression data
  # -2 means the structure_id matched multiple times
  
  expression <- rep(NA,length(names.id.key))
  names.id.key <- unlist(names.id.key)
  names(expression) <- names(names.id.key)
  for(region in names(names.id.key[!is.na(names.id.key)])){
    region.expression <- goi.exp[goi.exp$structure_id == names.id.key[[region]],expression.metric]
    if(length(region.expression) == 1){expression[region] <- region.expression}
    if(length(region.expression) == 0){expression[region] <- NA}
    if(length(region.expression) > 1){expression[region] <- -2}
  }
  return(expression)
}

tryDownload <- function(url,destfile){
  # INPUTS:
  # url: character of url to file
  # destfile: character of path for where to save file
  
  # OUTPUTS:
  # 'good' if it worked, 'bad' if it didn't
  out <- tryCatch({
    download.file(url = url,destfile = destfile)
  }, # try to download expression for gene of interest
  error= function(cond) {return('bad')})
  return(out)
}
