# Amanda Leonti
# 07 Aug, 2019
# Purpose: Take in a list of gene symbols, and identify gene aliases that are found in our counts data
##################################

identify_geneAlias <- function(genes){

  
  # This is OLD and I found a faster way to do it, only keeping this around in case I incorporate parts of it into the new function someday...
  get_aliases_old <- function(geneSymbol){
    
    # If the gene symbol fed into this function is NOT in the expression data already...
    if(!geneSymbol %in% allGenes$geneSymbol){
      print(paste0("Can't find ", geneSymbol, " in expression data, looking up aliases..."))
      
      # Identifying the row of the alias spreadsheet where the gene symbol IS found
      aliases <- aliasDF[apply(aliasDF, 1, function(i) any(grep(paste0("^\\b", geneSymbol, "\\b$"), i))),]
      
      # Only runs the steps below if at least 1 alias is returned
      if(nrow(aliases) >= 1){
        
        # Checking every alias in the row to identify the one that IS in our expression data
        if(any(apply(aliases, 2, function(x) x %in% allGenes$geneSymbol))){
          
          new_alias <- unname(unlist(aliases[apply(aliases, 2, function(x) x %in% allGenes$geneSymbol)]))
          
          # Checking to make sure it's in the expression data
          if(any(new_alias %in% allGenes$geneSymbol)){
            
            # If multiple aliases are identified, this will print a message for all aliases (instead of just the first one identified)
            sapply(new_alias, function(x) print(paste0("Alias found! ", geneSymbol, " is referred to as ", x, " in our expression data.")))
            
            results <- data.frame(Old.symbol = geneSymbol, 
                                  Final = new_alias)
            
            # Returns all identified aliases that are also in our expression data
            return(results)
            
          }else{
          
            # Also searching for the alias using the geneSynonym() package (our alias sheet isn't up-to-date and is sometimes missing aliases)
            synonyms <- unname(unlist(humanSyno(geneSymbol)))
            
            if(any(geneSymbol %in% synonyms)){
              
              new_alias <- synonyms[synonyms %in% allGenes$geneSymbol]
              
              print(paste0("Alias found! ", geneSymbol, " is referred to as ", new_alias, " in our expression data.")) 
              
              results <- data.frame(Old.symbol = geneSymbol, 
                                    Final = new_alias)
              return(results)
            }
            
          }
        
      }else{
        
        # Same deal as above
        print(paste0(geneSymbol, " is not in our expression data, nor does it have an alias that is."))
        
        results <- data.frame(Old.symbol = geneSymbol, 
                              Final = NA)
        
        return(results)
      }
      
    }else{   
      
      # If the gene symbol is already in our allGenes...
      print(paste0(geneSymbol, " is already in our expression data!")) # Skip the whole process above and print this message
      
      results <- data.frame(Old.symbol = geneSymbol, 
                            Final = geneSymbol)
      
      return(results)
    }
  
    }
    
  }
  
  
  # Defining a function to look up aliases for each gene in the gene list
  get_aliases_new <- function(geneSymbol){
    
    require(geneSynonym)
    
    # If the gene symbol is already in our expression data...
    if(geneSymbol %in% allGenes$geneSymbol){
      
    print(paste0(geneSymbol, " is already in our expression data!")) # Skip the whole process above and print this message
    
    results <- data.frame(Original = geneSymbol, 
                          Final = geneSymbol)
    
    return(results)
    }
    
    
    # If the gene symbol fed into this function is NOT in the expression data already...
    if(!geneSymbol %in% allGenes$geneSymbol){
      
      print(paste0("Can't find ", geneSymbol, " in expression data, looking up aliases..."))
      
      # Searching for a gene synonym with the geneSynonyms package
      synonyms <- unname(unlist(humanSyno(geneSymbol)))
      
      # The length requirement is needed because humanSyno() will return the original gene symbol if it can't find any aliases (and we already know the symbol isn't in the exp data)
      if(length(synonyms) > 1 & any(geneSymbol %in% synonyms)){
        
        new_alias <- synonyms[synonyms %in% allGenes$geneSymbol][1] # Sometimes this will return multiple (?), will only keep 1 for now!
        
        print(paste0("Alias found! ", geneSymbol, " is referred to as ", new_alias, " in our expression data.")) 
        
        results <- data.frame(Original = geneSymbol, 
                              Final = new_alias)
        return(results)
        
      }else{
        
        # Message printed when none of the synonyms identified by the geneSynonym package are in our expression data
        print(paste0(geneSymbol, " is not in our expression data, nor does it have an alias that is."))
        
        results <- data.frame(Original = geneSymbol, 
                              Final = geneSymbol)
        
        return(results)
      }
      
    }
    
  }
  
  # Reading in a matrix of gene symbols & their associated aliases - only needed for the old version of the code
  # aliasDF <- read.delim("/Volumes/homes/meshinchi_s/workingDir/TARGET/AML_TARGET/RNA/mRNAseq/analysis/0000.00.02_Reference_GeneInfo/Hs_GeneInfo.Symb.Alias.txt", header = F, stringsAsFactors = F)
  
  # Reading a list of ~51k gene symbols that are found in our expression data
  allGenes <- read.csv("All_GeneSymbols_in_our_Ribodepleted_RNAseq_Data_08.07.2019.csv", stringsAsFactors = F, strip.white = T)
  
  # Removing whitespaces
  genes <- trimws(genes, which = "both") # This doesn't always work, for whatever reason
  genes <- gsub("\\s$|^\\s", "", genes)
  
  # Identifying aliases for all genes in the "genes" list that are found in our expression data (the "allGenes" object)
  all_aliases <- lapply(genes, get_aliases_new)
  
  # Collapsing list of dataframes into 1 big dataframe
  final_table <- suppressWarnings(bind_rows(all_aliases)) # Will produce warnings about coercing a factor to character vector, which is harmless
  
  return(final_table)
  
}
