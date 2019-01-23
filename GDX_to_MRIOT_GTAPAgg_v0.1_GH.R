

# This script parses gdx files created by 'GTAPAgg' and builds multi-regional 
# input-output tables based on Peters, G. P., Andrew, R., & Lennox, J. (2011). Constructing an 
# environmentally-extended multi-regional input-output table using the GTAP database. Economic Systems 
# Research, 23(2), 131-152.

# Version: 0.1
# Author: David Font Vivanco
# Date: 23/01/2018

# Notes: 
# Example based on version 'GTAP9a_2011'.
# VDIM (investment domestic demand) and VIIM (investment imports demand) are endogenous as
# investments are considered a sector and not a demand category in the original data.

##################################################
# loading packages, files, defining variables, etc.

# load required packages
list.of.packages <- c('gdxtools')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if (length(new.packages))  install.packages(new.packages, dependencies = TRUE)
lapply(list.of.packages, require, character.only = TRUE)

# set working directories
wd_input = ".../input/GTAP9a_2011"
wd_output = ".../output"


setwd(wd_input)

# define a gdx
mygdx <- gdx('basedata.gdx')

# parameters, variables, and sets in gdx file
parameters = mygdx$parameters
sets = mygdx$sets

# defining dimensions and names

regions = unlist(mygdx['REG']) # region names
nr = length(regions) # number of regions

intermediates.CGDS = unlist(mygdx["PROD_COMM"]) # intermediate names, with CGDS
intermediates = intermediates.CGDS[!intermediates.CGDS %in% 'CGDS'] # removing investment
nintermediates = length(intermediates) # number of intermediates

demandnames = c('PRIV','GOVT','CGDS') # PRIV:private households; GOVT: government; CGDS: investment
ndemand = length(demandnames) # n of final demand types

transportnames = unlist(mygdx["MARG_COMM"]) # transport services names
ntransport = length(transportnames)

# list of required parameters to construct MRIOT

MRIOT.par = c('VDFM','VIFM','VXMD','VDPM','VDGM','VIPM','VIGM','VST','VTWR')

# extracting necessary parameters and completing values

pars = list()

for (i in MRIOT.par){
  
  data = as.data.frame(mygdx[i])

  # adding dimensions to parameters
  
  text = parameters$text[grepl(pattern = i,x = parameters$text)] # find complete string
  dims = strsplit(x = gsub(".*[:]([^.]+)[]]].*", "\\1", text),split = '\\*')[[1]] # list of dimensions

  colnames(data) <- c(dims,'value') # changing column names
  
  # adding missing entries to obtain full dimensions later on
  # it is not very elegant but it works...
  
  if (i == "VDFM" | i == "VIFM"){
    
    all.comb = expand.grid(TRAD_COMM=intermediates, PROD_COMM=intermediates.CGDS,REG=regions)
    data2 = merge(all.comb,data,all=TRUE)
    
  } else if (i == "VXMD"){
    
    colnames(data) = c('TRAD_COMM','REGr','REGs','value')
    all.comb = expand.grid(TRAD_COMM=intermediates,REGr=regions,REGs=regions)
    data2 = merge(all.comb,data,all=TRUE)
    
  } else if (i == "VDPM" | i == "VDGM" | i == 'VIPM' | i == 'VIGM'){
      
    all.comb = expand.grid(TRAD_COMM=intermediates,REG=regions)
    data2 = merge(all.comb,data,all=TRUE)
      
  } else if (i == 'VST'){
    
    all.comb = expand.grid(MARG_COMM=transportnames,REG=regions)
    data2 = merge(all.comb,data,all=TRUE)
    
  } else if (i == 'VTWR'){
    
    colnames(data) = c('MARG_COMM','TRAD_COMM','REGr','REGs','value')
    all.comb = expand.grid(MARG_COMM=transportnames,TRAD_COMM=intermediates,REGr=regions,REGs=regions)
    data2 = merge(all.comb,data,all=TRUE)
    
  }
    
  data2[is.na(data2)] = 0
  
  pars[[i]] = data2

}


# calculating VTW (aggregate international transportation services) for posterior use
pars[['VTW']] = aggregate(value ~  MARG_COMM, pars[['VTWR']], sum)


# Z matrix and Z labels

Z = matrix(0,nintermediates*nr,nintermediates*nr)
Zlabels = data.frame('PROD_COMM' = rep(intermediates,times=nr))
Zlabels$REG = rep(regions,each=nintermediates)

# Y matrix and Y labels

Y = matrix(0,nintermediates*nr,ndemand*nr)
Ylabels = data.frame('DEMAND' = rep(demandnames,nr))
Ylabels$REG = rep(regions,each=ndemand)


######################################################
# implementing model with endogenous transport margins

# domestic Z (diagonal elements of Z)

for (r in regions){ # for each region r
  
  # converting VDFM into matrices
  # VDFM: domestic firm purchases of i by j in region r
  VDFM = pars[['VDFM']]
  VDFMr = VDFM[ VDFM$REG == r,]
  
  VDFMr.mat = reshape(VDFMr[,-3], idvar = "TRAD_COMM", timevar = "PROD_COMM", direction = "wide")[,-1]
  VDFMr.mat[is.na(VDFMr.mat)] = 0
  
  # moving investments to final demand
  VDFMr.CGDS = VDFMr.mat[,which(colnames(VDFMr.mat) %in% 'value.CGDS')]
  VDFMr.mat = VDFMr.mat[,-which(colnames(VDFMr.mat) %in% 'value.CGDS')]
  
  # filling CGDS in Y
  cY = intersect( which(Ylabels$REG == r),which(Ylabels$DEMAND == 'CGDS') )
  rY = which(Zlabels$REG == r)
  Y[rY,cY] = VDFMr.CGDS
  
  # filling Z
  cZ = which(Zlabels$REG == r)
  rZ = cZ
  
  Z[rZ,cZ] = as.matrix(VDFMr.mat)
  
}

# domestic Y (diagonal elements of Y)

for (r in regions){ # for every region r
  
  # VD*M = domestic purchases
  
  # private households
  VDPM = pars[['VDPM']]
  VDPMr = VDPM[VDPM$REG == r,]
  
  VDPMr.mat = as.matrix(VDPMr$value)
  VDPMr.mat[is.na(VDPMr.mat)] = 0
  
  # filling Y
  cY = intersect( which(Ylabels$REG == r),which(Ylabels$DEMAND == 'PRIV') )
  rY = which(Zlabels$REG == r)
  Y[rY,cY] = VDPMr.mat
  
  # government
  VDGM = pars[['VDGM']]
  VDGMr = VDGM[VDGM$REG == r,]
  
  VDGMr.mat = as.matrix(VDGMr$value)
  VDGMr.mat[is.na(VDGMr.mat)] = 0
  
  # filling Y
  cY = intersect(which(Ylabels$REG == r),which(Ylabels$DEMAND == 'GOVT'))
  rY = which(Zlabels$REG == r)
  Y[rY,cY] = VDGMr.mat

}


# filling imports in Z and Y (off-diagonal elements Zrs and Yrs) and endogenising transport margins

endogenous = TRUE # endogenise transport margins?

for ( r in regions ){ # for each importing region r
  
  # we construct the off-diagonal terms by distributing the row sum VXMD over the columns 
  # using VIFM and this retains the various balances when margins are taken into account
  
  # extracting VIFM
  
  VIFM = pars[['VIFM']]
  VIFMr = VIFM[VIFM$REG == r,]
  
  # VIFM as distribution matrix
  
  VIFMr.mat = reshape(VIFMr[,-3], idvar = "TRAD_COMM", timevar = "PROD_COMM", direction = "wide",
                      new.row.names = unique(VIFMr$TRAD_COMM))[,-1]
  VIFMr.mat[is.na(VIFMr.mat)] = 0
  VIFMr.CGDS = VIFMr.mat[,intermediates.CGDS %in% 'CGDS']
  VIFMr.mat = VIFMr.mat[,!intermediates.CGDS %in% 'CGDS']
  
  # private households
  # extracting VIPM and converting to matrix
  
  VIPM = pars[['VIPM']]
  VIPMr = VIPM[VIPM$REG == r,]
  VIPMr.mat = as.matrix(VIPMr$value)
  VIPMr.mat[is.na(VIPMr.mat)] = 0
  
  # government
  # extracting VIGM and converting to matrix
  
  VIGM = pars[['VIGM']]
  VIGMr = VIGM[VIGM$REG == r,]
  VIGMr.mat = as.matrix(VIGMr$value)
  VIGMr.mat[is.na(VIGMr.mat)] = 0

  # VIM = VIFM + VIPM + VIGM + VIFM.CGDS (dim = nintemerdiates:(nintermediates+ndemand))
  
  VIMr = cbind(VIFMr.mat,VIPMr.mat,VIGMr.mat,VIFMr.CGDS)
  
  # VIM as distribution matrix: VIM.dist = (diag( rowsum(VIM) ))^-1 * VIM
  
  VIMr.dist = as.matrix(sweep(x = VIMr,MARGIN = 1,STATS = rowSums(VIMr),FUN = '/'))
  VIMr.dist[is.na(VIMr.dist)] = 0  
  
  # extracting VXMD
  VXMD = pars[['VXMD']]
  
  for (s in regions){ # for each exporting region s
    
    # extracting VXMD for region s
    VXMDs = VXMD[(VXMD$REGs == s & VXMD$REGr == r),]
    
    # VXMD as matrix
    VXMDs.mat = as.matrix(VXMDs$value)
    VXMDs.mat[is.na(VXMDs.mat)] = 0
    
    # Zrs and Yrs: distributing VXMD using the structure of VIM
    # Yrs includes CGDS, so imported CGDS is moved to final demand
    Zrs = diag(as.vector(VXMDs.mat)) %*% VIMr.dist[,1:nintermediates]
    Yrs = diag(as.vector(VXMDs.mat)) %*% VIMr.dist[,(nintermediates+1):ncol(VIMr.dist)]
    
    # Inserting Zrs and Yrs (off-diagonal elements)
    cZ = which(Zlabels$REG == r)
    rZ = which(Zlabels$REG == s)
    Z[rZ,cZ] = Z[rZ,cZ] + Zrs
    
    cY = which(Ylabels$REG == r)
    rY = which(Zlabels$REG == s)
    Y[rY,cY] = Y[rY,cY] + Yrs
    
    if (endogenous == T){
      
      # transport margins (VTWR) are endogenised
      
      # extracting VTWR from region s to region r and converting to matrix
      VTWR = pars[['VTWR']]
      VTWRsr = VTWR[(VTWR$REGs == s & VTWR$REGr == r),]
      VTWRsr = aggregate(value ~ TRAD_COMM,data=VTWRsr ,'sum') # transport types are aggregated
      VTWRsr.mat = as.matrix(VTWRsr$value)
      
      # distribution of transport margins (VTWR) over users (only firms)
      
      # VIFM as distribution matrix
      
      VIFMr.dist = as.matrix(sweep(x = VIFMr.mat,MARGIN = 1,STATS = rowSums(VIFMr.mat),FUN = '/'))
      VIFMr.dist[is.na(VIFMr.dist)] = 0 
      
      # distributing VTWR according to the structure of VIFM
      U.Zrs = diag(as.vector(VTWRsr.mat)) %*% VIFMr.dist
      
      # transport pool (VTW)
      VTW = pars[['VTW']]
      VTW.mat = as.matrix(VTW$value)

      # transport margins are distributed over suppliers
      
      for (ss in regions){ # for each region ss that supplies transport for the exports of region s to region r
        
        # transport supply (VST)
        VST = pars[['VST']]
        VSTss = VST[VST$REG == ss,]
        VSTss.mat = as.matrix(VSTss$value)

        # transport pool, VTW, is distributed evenly in proportion to supply, VST
        VTW.dist = VSTss.mat/sum(VTW.mat)
        VTW.dist[is.na(VTW.dist)] = 0
        
        # transport supply over users
        T.Zrs = VTW.dist %*% colSums(U.Zrs)
        
        # inserting T.Zrs
        cZ = which(Zlabels$REG == r)
        rZ = which(Zlabels$REG == ss & Zlabels$PROD_COMM == transportnames)
        Z[rZ,cZ] = Z[rZ,cZ] + T.Zrs
        
      }
    }
  }
  
  print(r)
}


# calculating A and x

x = rowSums(Z) + rowSums(Y)
A = sweep(x = Z,MARGIN = 2,STATS = x,FUN = '/')

# saving objects

setwd(wd_output)
saveRDS(A,'A.rds')
saveRDS(x,'x.rds')


# checking that all parameters are fully represented in the input-output system

VDFM = pars[['VDFM']]

VXMD = pars[['VXMD']]

VDPM = pars[['VDPM']]

VDGM = pars[['VDGM']]

VTWR = pars[['VTWR']]

# x = VDFM + VXMD + VD*M + VTWR
x1 = sum(x)
x2 = sum(VDFM$value) + sum(VXMD$value) + sum(VDPM$value) + sum(VDGM$value) + sum(VTWR$value)
all.equal(x1,x2)

# end