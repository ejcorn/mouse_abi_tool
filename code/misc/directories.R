rm(list=setdiff(ls(),'params'))
basedir <- params$basedir
setwd(basedir)
# create directories for each set of analyses in advance

dir.create(paste(params$opdir,'diffmodel/',sep=''),recursive=TRUE)