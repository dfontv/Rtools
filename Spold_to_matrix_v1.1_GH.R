

# This script parses spold files from an Ecoinvent system and creates calculation-ready matrices
# It is largely based on the presentation 'Building matrices from spold files' by Guillaume Bourgault

# Author: David Font Vivanco
# Date: 04/12/2018
# Version: 1.1

# Fixes respect to previous version:
# 1. both child and parent datasets are correctly parsed


## install & require libraries ##
list.of.packages <- c('XML','tictoc',"data.table",'Matrix','Rcpp','RcppArmadillo')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if (length(new.packages))  install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)


# Working directories (based on the Ecoinvent v3 files downloaded from the web)
wd_datasets = ".../datasets"
wd_MasterData = ".../MasterData"
wd_output = ".../output"


##########################
# Step 1: building indexes
##########################

# A labels
setwd(wd_datasets)
files = list.files()

A.labels = data.frame( do.call('rbind',strsplit(x = gsub(pattern = '.spold',replacement = '',x = files),
                                                split = '_')) )
A.labels$index = 1:nrow(A.labels)
colnames(A.labels)[1:2] = c('activity id','Reference product intermediateExchangeId')


# B labels
setwd(wd_MasterData)

elementary.filename = 'ElementaryExchanges.xml'
elementary.parsed <- xmlTreeParse(elementary.filename,useInternalNodes = T) # parsing xml file
elementary.parsed.list = xmlToList(elementary.parsed) # xml to list
B.list = list()

# Exchanges to list
for ( i in 1:length(elementary.parsed.list) ){
  
  exchange = elementary.parsed.list[i]
  elementaryExchangeId = exchange$elementaryExchange$.attrs[1]
  subcompartmentId = exchange$elementaryExchange$compartment$.attrs
  index = i
  names(index) = 'index'
  if( length(c(elementaryExchangeId,subcompartmentId,index)) == 3 ){
    B.list[[i]] = c(elementaryExchangeId,subcompartmentId,index)}
  
  
}

B.labels = data.frame(do.call('rbind',B.list))


run.step.1.2 = TRUE # step 1.2 should be only run once as it may take up to a few hours

if (run.step.1.2 == T){
  
  ######################################################
  # Step 1.2: extracting and storing data (run only once)
  ######################################################  
  
  # Instead of extracting matrix values in each file, target variables are first extracted into a single dataset
  # This allows to backtrack information from A and B matrices and store additional information (e.g.
  # prices and production volume as in the example)
  
  setwd(wd_datasets)
  
  activity.description = list() # additional activity description for each file in a single list
  flow.Data = list() # all flow data for each file in a single list
  
  nflow = 1 # flow counter 
  nactivity = 1 # activity counter
  elapsed = list() # time counter
  elapsed.total = 0  # time counter
  
  # This loop extracts additional activity description and all flow data plus price and production 
  # volume data (can be adapted to extract other desired variables)
  
  for ( f in 1:length(files) ){ # for every spold file
    
    start = Sys.time()
    
    # resetting targeted variables each file
    classificationId = NA
    geographyId = NA
    
    filename = gsub(pattern = '.spold',replacement = '',x = files[f])
    xml = files[f]
    Parsed.xml <- xmlTreeParse(xml,useInternalNodes = T) # parsing xml file
    xml.list = xmlToList(Parsed.xml) # xml to list
    
    if ( is.null(xml.list$childActivityDataset)==F ){ # if child dataset
      
      classificationId = xml.list$childActivityDataset$activityDescription$classification$.attrs
      geographyId = xml.list$childActivityDataset$activityDescription$geography$.attrs
      
    } else if ( is.null(xml.list$activityDataset)==F ){ # if parent dataset
      
      classificationId = xml.list$activityDataset$activityDescription$classification$.attrs
      geographyId = xml.list$activityDataset$activityDescription$geography$.attrs
      
    }
    
    activity.description[[nactivity]] = c(classificationId, geographyId,f)
    nactivity = nactivity + 1
    
    for (ff in 1:1000000){ # for as many flows as there are
      
      if ( is.null(xml.list$childActivityDataset)==F ){ # if child dataset
        
        flow.data = xml.list$childActivityDataset$flowData[ff]
        
      } else if ( is.null(xml.list$activityDataset)==F ){ # if parent dataset
        
        flow.data = xml.list$activityDataset$flowData[ff]
        
      }
      
      
      # resetting targeted variables each flow
      Exchange = NA
      Group = NA
      Groupnumber = NA
      intermediateExchangeId = NA
      activityLinkId = NA
      amount = NA
      productionVolumeAmount = NA
      price_EUR2005 = NA
      
      if ( is.null(unlist(flow.data)) == T ){ break } # if no more flows, then break
      
      if ( is.null(flow.data$intermediateExchange) == F ){ # if intermediate Exchange
        
        names.att = names(flow.data$intermediateExchange$.attrs) # names of attributes
        
        Exchange = 'intermediate'
        ExchangeId = flow.data$intermediateExchange$.attrs[names.att == 'intermediateExchangeId']
        activityLinkId = flow.data$intermediateExchange$.attrs[names.att == 'activityLinkId']
        amount = flow.data$intermediateExchange$.attrs[names.att == 'amount']
        productionVolumeAmount = flow.data$intermediateExchange$.attrs[names.att == 'productionVolumeAmount']
        if ( length(activityLinkId) == 0){activityLinkId = NA}
        if ( length(productionVolumeAmount) == 0 ){ productionVolumeAmount = NA}
        
        if ( is.null(flow.data$intermediateExchange$inputGroup) == F ){ # if input
          
          Group = 'input'
          Groupnumber = flow.data$intermediateExchange$inputGroup
          
        } else if ( is.null(flow.data$intermediateExchange$outputGroup) == F ){ # if output
          
          Group = 'output'
          Groupnumber = flow.data$intermediateExchange$outputGroup
          
        }
        
        # finding price
        
        for ( i in 1:1000000 ){ # for as many items as there are in flow ff
          
          if ( is.null(unlist(flow.data$intermediateExchange[i])) == T ){break} # if no more items then break
          
          if ( is.null(flow.data$intermediateExchange[i]$property) == F ){ # if it is a property
            
            names.att = names(flow.data$intermediateExchange[i]$property$.attrs)
            propertyId = flow.data$intermediateExchange[i]$property$.attrs[names.att == 'propertyId']
            
            if ( propertyId == '38f94dd1-d5aa-41b8-b182-c0c42985d9dc' ){ # if property is 'price_EUR2005'
              
              price_EUR2005 = flow.data$intermediateExchange[i]$property$.attrs[names.att == 'amount']
              
            }
            
          }
          
        }
        
        flow.Data[[nflow]] = c(filename, Exchange, Group, Groupnumber, ExchangeId, activityLinkId, amount,
                               price_EUR2005, productionVolumeAmount) 
        nflow = nflow + 1 
        
      } else if ( is.null(flow.data$elementaryExchange) == F ){ # if elementary Exchange
        
        
        names.att = names(flow.data$elementaryExchange$.attrs) # names of attributes
        
        Exchange = 'elementary'
        ExchangeId = flow.data$elementaryExchange$.attrs[names.att == 'elementaryExchangeId']
        activityLinkId = NA
        amount = flow.data$elementaryExchange$.attrs[names.att == 'amount']
        productionVolumeAmount = NA
        
        if ( is.null(flow.data$elementaryExchange$inputGroup) == F ){ # if input
          
          Group = 'input'
          Groupnumber = flow.data$elementaryExchange$inputGroup
          
        } else if ( is.null(flow.data$elementaryExchange$outputGroup) == F ){ # if output
          
          Group = 'output'
          Groupnumber = flow.data$elementaryExchange$outputGroup
          
        }
        
        flow.Data[[nflow]] = c(filename, Exchange, Group, Groupnumber, ExchangeId, activityLinkId, amount,
                               price_EUR2005, productionVolumeAmount) 
        nflow = nflow + 1 
      }
      
    }
    
    end. = Sys.time()
    elapsed.f = as.numeric(difftime(time1 = end., time2 = start, units = "secs"))
    elapsed[[f]] = elapsed.f
    average.time = mean( unlist(elapsed) )
    elapsed.total = elapsed.total + (elapsed.f/(60*60))
    
    estimated =  round( ((mean( unlist(elapsed) )/(60*60) ) * length(files) - elapsed.total),2 )
    
    print( paste( round(f/length(files) * 100,digits = 2),'% completed',' | ',round(elapsed.total,2), ' hours elapsed, ',
                  estimated,' hours remaining',sep='' ) )
  }
  
  # converting lists to data frames
  flow.Data.df = data.frame(do.call('rbind',flow.Data))
  colnames(flow.Data.df) = c('filename', 'Exchange', 'Group', 'Groupnumber', 'ExchangeId', 'activityLinkId',
                             'amount', 'price_EUR2005', 'productionVolumeAmount')
  
  activity.description.df = data.frame(do.call('rbind',activity.description))
  colnames(activity.description.df) = c('classificationId', 'geographyId')
  
  # saving outputs as R objects
  setwd(wd_output)
  
  saveRDS(activity.description.df, 'activity.description.df.rds')
  saveRDS(flow.Data.df, 'flow.Data.df.rds')
  
}

#################################
# Step 2: populating the matrices
#################################

# Loading R objects
setwd(wd_output)

activity.description.df = readRDS('activity.description.df.rds')
flow.Data.df = readRDS('flow.Data.df.rds')

# Preparing to build A and B matrices
# First, we calculate row and column indices for each value, and modify sign of value if needed

row.col.value = list()

tic()
for ( f in 1:nrow(flow.Data.df) ){ # for each intermediate/elementary flow
  
  # resetting target variables
  col.index = NA
  row.index = NA
  value = NA
  
  activity_id = strsplit(as.character(flow.Data.df$filename[f]),split = '_')[[1]][1]
  Reference_product_intermediateExchangeId = strsplit(as.character(flow.Data.df$filename[f]),split = '_')[[1]][2]
  
  col.index = intersect( which( A.labels$`activity id` == activity_id ),
                         which( A.labels$`Reference product intermediateExchangeId` == Reference_product_intermediateExchangeId ))
  
  # if reference product (outputGroup = 0)
  if ( flow.Data.df$Exchange[f] == 'intermediate' && flow.Data.df$Group[f] == 'output' &&
       flow.Data.df$Groupnumber[f] == '0' ){
    
    value = as.numeric(as.character(flow.Data.df$amount[f]))
    row.index = col.index # diagonal element
    
    # if intermediate input (inputGroup = 5)
  } else if ( flow.Data.df$Exchange[f] == 'intermediate' && flow.Data.df$Group[f] == 'input' &&
              flow.Data.df$Groupnumber[f] == '5' ){
    
    # The amount sign SHOULD be changed, because in the matrix A, the sign convention is positive 
    # for production, and negative for consumption. A negative amount should be changed to positive.
    value = -as.numeric( as.character(flow.Data.df$amount[f]) )
    activityLinkId = as.character(flow.Data.df$activityLinkId[f])
    ExchangeId = as.character(flow.Data.df$ExchangeId[f])
    
    if ( !is.na(activityLinkId) && value != 0 ){
      
      row.index = which( A.labels$`activity id` == activityLinkId & 
                           A.labels$`Reference product intermediateExchangeId` == ExchangeId )
    }
    
    
    # if elementary exchange  (inputGroup = 4 and outputGroup = 4)
  } else if ( flow.Data.df$Exchange[f] == 'elementary' && flow.Data.df$Groupnumber[f] == '4' ){
    
    value = as.numeric( as.character(flow.Data.df$amount[f]) )
    elementaryExchangeId = as.character(flow.Data.df$ExchangeId[f])
    row.index = which(B.labels$id == elementaryExchangeId)
    
  }
  
  row.col.value[[f]] = c(row.index,col.index,value)
  
  print( paste( round(f/nrow(flow.Data.df) * 100,digits = 2),'% completed',sep='' ) )
  
}  
toc() 

# adding additional variables
flow.Data.df = cbind( flow.Data.df, data.frame(do.call('rbind',row.col.value)) )
colnames(flow.Data.df)[10:12] = c('row.index','col.index','value')

# subsetting valid intermediate and elementary exchanges  
intermediate.exchanges = flow.Data.df[ (flow.Data.df$Exchange == 'intermediate' & 
                                          is.na(flow.Data.df$value) == F & 
                                          is.na(flow.Data.df$row.index) == F),]
elementary.exchanges = flow.Data.df[ (flow.Data.df$Exchange == 'elementary' & 
                                        is.na(flow.Data.df$value) == F & 
                                        is.na(flow.Data.df$row.index) == F),]

# Using row and column indexes, we build compressed, column-oriented, sparse matrices of class "dgCMatrix"
# Sparse matrices have optimized size and can be used in computations in most cases
# Sparse matrices can be converted to (dense) numeric with 'as.matrix()' if needed
A.dgC <- sparseMatrix(i=intermediate.exchanges$row.index, j=intermediate.exchanges$col.index, 
                      x = intermediate.exchanges$value, dims = c(nrow(A.labels),nrow(A.labels)))
B.dgC <- sparseMatrix(i=elementary.exchanges$row.index, j=elementary.exchanges$col.index, 
                      x = elementary.exchanges$value, dims = c(nrow(B.labels),nrow(A.labels)))

# adding additional information to A.labels (classification and geography)
A.labels = cbind(A.labels, activity.description.df)

# Saving R objects
setwd(wd_output)

saveRDS(flow.Data.df, 'flow.Data.df2.rds')

saveRDS(A.dgC, 'A.dgC.rds')
saveRDS(A.labels, 'A.labels.rds')

saveRDS(B.dgC, 'B.dgC.rds')
saveRDS(B.labels, 'B.labels.rds')


###################
# Step 3: Calculate
###################

setwd(wd_output)

A.dgC = readRDS('A.dgC.rds')
A.labels = readRDS('A.labels.rds')

B.dgC = readRDS('B.dgC.rds')
B.labels = readRDS('B.labels.rds')

# Defining Rcpp function to speed up matrix algebra
Rcpp::cppFunction("arma::mat armaInv(arma::mat x) { return arma::inv(x); }", depends="RcppArmadillo")

# (dummy) demand vector
y = matrix(0,nrow(A.dgC),1)
y[1] = 1

# Leontief inverse
# If many calculations are planned, it is recommended to solve the system through the Leontief inverse
# Otherwise, it will be quicker to solve the system Ax=y, and even quicker through iterative methods (e.g. power series)
tic()
L = armaInv( as.matrix(A.dgC) )
toc()

# total output
tic()
x = L %*% y
toc()

# inventory
tic()
g = as.matrix(B.dgC) %*% x
toc()


