#################
### Load data ###
#################

rm(list=setdiff(ls(),'params'))
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'processed/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/miscfxns.R')

connectivity.ipsi <- read.csv('data/Connectome_Ipsi.csv',row.names = 1,header = TRUE,check.names = FALSE)
connectivity.contra <- read.csv('data/Connectome_Contra.csv',row.names = 1,header = TRUE, check.names = FALSE)

###################################
### Process connectivity matrix ###
###################################

conn.names.ipsi <- colnames(connectivity.ipsi)
conn.names.contra <- colnames(connectivity.contra)

# checks 
if(identical(colnames(connectivity.contra),rownames(connectivity.contra))){
  print('contra connectivity colnames and rownames equal')}
if(identical(colnames(connectivity.ipsi),rownames(connectivity.ipsi))){
  print('ipsi connectivity colnames and rownames equal')}
region.names.hemi <- colnames(connectivity.ipsi)

W <- rbind(cbind(connectivity.ipsi,connectivity.contra),cbind(connectivity.contra,connectivity.ipsi))
rownames(W) <- c(paste('i',rownames(connectivity.ipsi),sep=''), # add i to ipsilateral regions
                 paste('c',rownames(connectivity.contra),sep='')) # add c to contralateral regions
colnames(W) <- c(paste('i',rownames(connectivity.ipsi),sep=''), # add i to ipsilateral regions
                 paste('c',rownames(connectivity.contra),sep='')) # add c to contralateral regions

n.regions.ABA <- nrow(W)
n.regions.ABA.hemi <- n.regions.ABA/2
# check if connectivity was tiled into alternating blocks correctly

unit.test(all(W[(n.regions.ABA.hemi+1):n.regions.ABA,(n.regions.ABA.hemi+1):n.regions.ABA] == W[1:n.regions.ABA.hemi,1:n.regions.ABA.hemi]),
          'tiling on-diagonal blocks worked','tiling on-diagonal blocks failed') 
unit.test(all(W[1:n.regions.ABA.hemi,(n.regions.ABA.hemi+1):n.regions.ABA] == W[(n.regions.ABA.hemi+1):n.regions.ABA,1:n.regions.ABA.hemi]),
          'tiling off-diagonal blocks worked','tiling off-diagonal blocks failed')

unit.test(all(colnames(W) == rownames(W)),'row and column names of conn mat are same','ERROR with conn mat names')
region.names <- colnames(W) # store connectivity matrix names

save(region.names, region.names.hemi, n.regions.ABA, n.regions.ABA.hemi,file = paste(savedir,'connectome.RData',sep=''))
write.csv(x=W,file=paste0(basedir,'data/W.csv'))
writeMat(paste(savedir,'W.mat',sep=''),W=as.matrix(W))