

# This script parses gdx files from gtapinGAMS (http://www.mpsge.org/gtap6/) and builds multi-regional 
# input-output tables based on Peters, G. P., Andrew, R., & Lennox, J. (2011). Constructing an 
# environmentally-extended multi-regional input-output table using the GTAP database. Economic Systems 
# Research, 23(2), 131-152.

# Version: 0.3
# Author: David Font Vivanco
# Date: 22/11/2018

##################################################
# loading packages, files, defining variables, etc.

# load required packages
list.of.packages <- c('gdxtools')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if (length(new.packages))  install.packages(new.packages, dependencies = T)
lapply(list.of.packages, require, character.only = TRUE)

# set working directory
wd = "mywd"
setwd(wd)

# define a gdx
mygdx <- gdx('mygdx.gdx')

# parameters, variables, and sets in gdx file
parameters = mygdx$parameters
variables = mygdx$variables
sets = mygdx$sets

# list of required parameters to construct MRIOT

MRIOT.par = c('vdfm','vifm','vxmd','vxmd','vdpm','vdgm','vdim','vipm','vigm','viim','vst','vtw','vtwr')

# extracting necessary parameters

pars = list()

for (i in MRIOT.par){
  
  name = paste(toupper(i),'T',sep='') # parameter name as reported
  pars[[i]] = mygdx[name]
  
}

# defining dimensions and names

regions = unlist(as.list(mygdx['r'])) # region names
nr = length(regions) # number of regions

intermediates = unique(mygdx["vdfm"]$i) # intermediate names
nintermediates = length(intermediates) # number of intermediates

demandnames = c('PRIV','GOVT','CGDS') # PRIV:private households; GOVT: government; CGDS:investment goods
ndemand = length(demandnames) # n of final demand types

transportnames = unique(mygdx["vst"]$i) # transport services names
ntransport = length(transportnames)


# Z matrix and Z labels

Z = matrix(0,nintermediates*nr,nintermediates*nr)
Zlabels = data.frame(i = rep(intermediates,times=nr))
Zlabels$r = rep(regions,each=nintermediates)

# Y matrix and Y labels

Y = matrix(0,nintermediates*nr,ndemand*nr)
Ylabels = data.frame(i = rep(demandnames,nr))
Ylabels$r = rep(regions,each=ndemand)


#######################################################
# implementing model with endogenised transport margins

year = 2011 # this applies only if thre is a time dimension in the gdx data

# domestic Z (diagonal elements of Z)

for (r in regions){ # for each region r

  # converting vdfm into matrices
  # vdfm: domestic firm purchases of i by j in region r
  vdfm = pars[['vdfm']]
  vdfmr = vdfm[(vdfm$t == year & vdfm$r == r),]

  vdfmr.mat = reshape(vdfmr[,-c(3:4)], idvar = "i", timevar = "j", direction = "wide")[,-1]
  vdfmr.mat[is.na(vdfmr.mat)] = 0
  
  # filling Z
  cZ = which(Zlabels$r == r)
  rZ = cZ

  Z[rZ,cZ] = as.matrix(vdfmr.mat)
  
}

# domestic Y (diagonal elements of Y)

for (r in regions){ # for every region r
  
  # vd*m = domestic purchases
  
  # private households
  vdpm = pars[['vdpm']]
  vdpmr = vdpm[(vdpm$t == year & vdpm$r == r),]

  vdpmr.mat = as.matrix(vdpmr$value)
  vdpmr.mat[is.na(vdpmr.mat)] = 0
  
  # filling Y
  cY = intersect(which(Ylabels$r == r),which(Ylabels$i == 'PRIV'))
  rY = intersect(which(Zlabels$r == r),which(Zlabels[,1] %in% vdpmr$j))
  Y[rY,cY] = vdpmr.mat
  
  # government
  vdgm = pars[['vdgm']]
  vdgmr = vdgm[(vdgm$t == year & vdgm$r == r),]

  vdgmr.mat = as.matrix(vdgmr$value)
  vdgmr.mat[is.na(vdgmr.mat)] = 0
  
  # filling Y
  cY = intersect(which(Ylabels$r == r),which(Ylabels$i == 'GOVT'))
  rY = intersect(which(Zlabels$r == r),which(Zlabels[,1] %in% vdgmr$j))
  Y[rY,cY] = vdgmr.mat
  
  # investment
  vdim = pars[['vdim']]
  vdimr = vdim[(vdim$t == year & vdim$r == r),]
  vdimr.mat = as.matrix(vdimr$value)
  vdimr.mat[is.na(vdimr.mat)] = 0
  
  # filling Y
  cY = intersect(which(Ylabels$r == r),which(Ylabels$i == 'CGDS'))
  rY = intersect(which(Zlabels$r == r),which(Zlabels[,1] %in% vdimr$j))
  Y[rY,cY] = vdimr.mat
 
}


# filling imports in Z and Y (off-diagonal elements Zrs and Yrs) and endogenising transport margins

endogenous = T # endogenise transport margins?

for ( r in regions){ # for each importing region r
  
  # we construct the off-diagonal terms by distributing the row sum vxmd over the columns 
  # using vifm and this retains the various balances when margins are taken into account
  
  # extracting vifm
  
  vifm = pars[['vifm']]
  vifmr = vifm[(vifm$t == year & vifm$r == r),]

  # vifm as distribution matrix
  
  vifmr.mat = reshape(vifmr[,-c(3:4)], idvar = "i", timevar = "j", direction = "wide",
                     new.row.names = unique(vifmr$i))[,-1]
  names(vifmr.mat) = gsub("value.", "", names(vifmr.mat))
  vifmr.mat2 = matrix(0,nintermediates,nintermediates) # zero matrix will full dimensions
  
  # add missing rows and cols if necessary to retain full dimensions
  if (nrow(vifmr.mat) != nintermediates){
    
    vifmr.mat2[intermediates %in% rownames(vifmr.mat),] = 
      as.matrix(vifmr.mat[intermediates %in% rownames(vifmr.mat),])
  } 
  
  if (ncol(vifmr.mat) != nintermediates){
    
    vifm.mat2[intermediates %in% rownames(vifm.mat),] = 
      as.matrix(vifm.mat[,intermediates %in% colnames(vifm.mat)])
  }
  
  vifmr.mat2[is.na(vifmr.mat2)] = 0
  
  # private households
  # extracting vipm and converting to matrix
  
  vipm = pars[['vipm']]
  vipmr = vipm[(vipm$t == year & vipm$r == r),]
  vipmr.mat = as.matrix(vipmr$value)
  vipmr.mat[is.na(vipmr.mat)] = 0
  
  # zero vector will full dimensions
  vipmr.mat2 = matrix(0,nintermediates,1)
  
  # vipm with full dimensions
  vipmr.mat2[intermediates %in% unique(vipmr$j)] = vipmr.mat
  
  # government
  # extracting vigm and converting to matrix
  
  vigm = pars[['vigm']]
  vigmr = vigm[(vigm$t == year &vigm$r == r),]
  vigmr.mat = as.matrix(vigmr$value)
  vigmr.mat[is.na(vigmr.mat)] = 0
  
  # zero vector will full dimensions
  vigmr.mat2 = matrix(0,nintermediates,1)
  
  # vigm with full dimensions
  vigmr.mat2[intermediates %in% unique(vigmr$j)] = vigmr.mat
  
  # investment
  # extracting viim and converting to matrix
  
  viim = pars[['viim']]
  viimr = viim[(viim$t == year & viim$r == r),]
  viimr.mat = as.matrix(viimr$value)
  viimr.mat[is.na(viimr.mat)] = 0
  
  # zero vector will full dimensions
  viimr.mat2 = matrix(0,nintermediates,1)
  
  # viim with full dimensions
  viimr.mat2[intermediates %in% unique(viimr$j)] = viimr.mat
  
  
  # vim = vifm + vipm + vigm + viim (dim = nintemerdiates:(nintermediates+ndemand))
  
  vimr = cbind(vifmr.mat2,vipmr.mat2,vigmr.mat2,viimr.mat2)
  
  # vim as distribution matrix: vim.dist = (diag( rowsum(vim) ))^-1 * vim
  
  vimr.dist = sweep(x = vimr,MARGIN = 1,STATS = rowSums(vimr),FUN = '/')
  vimr.dist[is.nan(vimr.dist)] = 0  

  # extracting vxmd
  vxmd = pars[['vxmd']]

  for (s in regions){ # for each exporting region s

      # extracting vxmd for region s
      vxmds = vxmd[(vxmd$t == year & vxmd$s == s & vxmd$r == r),]
      
      # vxmd as matrix
      vxmds.mat = as.matrix(vxmds$value)
      vxmds.mat[is.na(vxmds.mat)] = 0
      
      # zero vector will full dimensions
      vxmds.mat2 = matrix(0,nintermediates,1)
      
      # vxmd with full dimensions
      vxmds.mat2[intermediates %in% unique(vxmds$i),1] = vxmds.mat
      
      # Zrs and Yrs: distributing vxmd using the structure of vim
      Zrs = diag(as.vector(vxmds.mat2)) %*% vimr.dist[,1:nintermediates]
      Yrs = diag(as.vector(vxmds.mat2)) %*% vimr.dist[,(nintermediates+1):ncol(vimr.dist)]
      
      # Inserting Zrs and Yrs (diagonal elements)
      cZ = which(Zlabels$r == r)
      rZ = which(Zlabels$r == s)
      Z[rZ,cZ] = Z[rZ,cZ] + Zrs
      
      cY = which(Ylabels$r == r)
      rY = which(Zlabels$r == s)
      Y[rY,cY] = Y[rY,cY] + Yrs
      
      if (endogenous == T){
        
        # transport margins (vtwr) are endogenised
        
        # extracting vtwr from region s to region r and converting to matrix
        vtwr = pars[['vtwr']]
        vtwrsr = vtwr[(vtwr$t == year & vtwr$s == s & vtwr$r == r),]
        vtwrsr = aggregate(value ~ i,data=vtwrsr ,'sum') # transport types are aggregated
        
        vtwrsr.mat = as.matrix(vtwrsr$value)
        
        print(paste('vtwrsr',s,r,sum(vtwrsr.mat)))
        
        # zero vector with full dimensions
        vtwrsr.mat2 = matrix(0,nintermediates,1)
        
        # vtwr with full dimensions
        vtwrsr.mat2[intermediates %in% unique(vtwrsr$i),1] = vtwrsr.mat
        
        # distribution of transport margins (vtwr) over users (only firms)
        
        # vifm as distribution matrix
        
        vifmr.dist = sweep(x = vifmr.mat2,MARGIN = 1,STATS = rowSums(vifmr.mat2),FUN = '/')
        vifmr.dist[is.nan(vifmr.dist)] = 0 
        
        # distributing vtwr according to the structure of vifm
        U.Zrs = diag(as.vector(vtwrsr.mat2)) %*% vifmr.dist

        # transport pool (vtw)
        vtw = pars[['vtw']]
        vtw = vtw[vtw$t == year,]
        vtw.mat = as.matrix(vtw$value)
        
        # zero vector with full dimensions
        vtw.mat2 = matrix(0,nintermediates,1)
        
        # vtw with full dimensions
        vtw.mat2[intermediates %in% unique(vtw$j),1] = vtw.mat
        
        # transport margins are distributed over suppliers
        
        for (ss in regions){ # for each region ss that supplies transport for the exports of region s to region r
          
          # transport supply (vst)
          vst = pars[['vst']]
          vstss = vst[(vst$t == year & vst$r == ss),]
          vstss.mat = as.matrix(vstss$value)
          
          # zero vector with full dimensions
          vstss.mat2 = matrix(0,nintermediates,1)
          
          # vst with full dimensions
          vstss.mat2[intermediates %in% unique(vstss$j),1] = vstss.mat
          
          # transport pool, vtw, is distributed evenly in proportion to supply, vst
          vtw.dist = vstss.mat2/sum(vtw.mat2)
          vtw.dist[is.na(vtw.dist)] = 0

          # transport supply over users
          T.Zrs = vtw.dist %*% colSums(U.Zrs)
          
          # inserting T.Zrs
          cZ = which(Zlabels$r == r)
          rZ = which(Zlabels$r == ss)
          Z[rZ,cZ] = Z[rZ,cZ] + T.Zrs

        }
      }
  }
}


# calculating A and x

x = rowSums(Z) + rowSums(Y)
A = sweep(x = Z,MARGIN = 2,STATS = x,FUN = '/')

# checks

vdfm = pars[['vdfm']]
vdfm = vdfm[(vdfm$t == year),]

vxmd = pars[['vxmd']]
vxmd = vxmd[(vxmd$t == year),]

vdpm = pars[['vdpm']]
vdpm = vdpm[(vdpm$t == year),]

vdgm = pars[['vdgm']]
vdgm = vdgm[(vdgm$t == year),]

vdim = pars[['vdim']]
vdim = vdim[(vdim$t == year),]

vtwr = pars[['vtwr']]
vtwr = vtwr[(vtwr$t == year),]

# x = vdfm + vxmd + vd*m + vtwr
x1 = sum(x)
x2 = sum(vdfm$value) + sum(vxmd$value) + sum(vdpm$value) + sum(vdgm$value) + sum(vdim$value) + sum(vtwr$value)
all.equal(x1,x2)

# end

