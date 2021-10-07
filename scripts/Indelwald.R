# THE CODE BELOW IS A COPY OF THE SCRIPT 3_Indelwald_spectrum.R,
# PART OF INDELWALD, WRITTEN BY MAXIMILIAN STAMMNITZ:
# https://github.com/MaximilianStammnitz/Indelwald

# Input:
#  'x' - character matrix of variants, with columns named CHROM, POS, REF, ALT
#  'reference' - DNAStringSet object containing reference genome
#                (from function Biostrings::readDNAStringSet)

# Modifications to the original code are indicated by the tag:
# ### MODIFIED (ao7)


############ INDELWALD - INDEL SPECTRUM ##############

## Last Update - 08/07/2021 ##
## mrs72 / Max Stammnitz ##
## maxrupsta@gmail.com ##
## University of Cambridge  ##

### MODIFIED (ao7) - avoid vcfR, stringr
# library(vcfR)
# library(Biostrings)
# library(stringr)
library(Biostrings)
###

### MODIFIED (ao7) - replace 'str_split_fixed' with 'strsplit' alternative
str_split_fixed <- function(string, sep, n) {
    t(sapply(strsplit(string, sep), function(x) {
        if (length(x) < n) x <- c(x, rep("", n-length(x)))
        if (length(x) > n) x <- c(x[1:(n-1)], paste(x[n:length(x)], collapse=sep))
        x
    }))
}
###

## Function to produce an indel spectrum
indel.spectrum <- function(x, reference){
  
  ### MODIFIED (ao7) - force to matrix, add TRIPLET column, filter problematic indels
  stopifnot(all.equal(colnames(x), c('CHROM', 'POS', 'REF', 'ALT')))
  x <- as.matrix(cbind(x, 'TRIPLET'=NA))
  pos <- as.integer(x[,'POS'])
  names(reference) <- sapply(strsplit(names(reference), " "), `[`, 1)
  idx1 <- (pos < 1000) | (pos > width(reference)[match(x[,'CHROM'], names(reference))] - 1000)
  idx2 <- nchar(x[,'REF']) > 50 | nchar(x[,'ALT']) > 50
  if (any(idx1 | idx2)) {
    x <- x[!(idx1 | idx2),]
    cat(' Indelwald discarded', sum(idx1 | idx2), 'indels\n')
  }
  ###
  
  ## 1. split VCF indels into:
  # (i) 1-bp deletions
  dels.1bp <- x[nchar(as.character(x[,'ALT'])) == nchar(as.character(x[,'REF']))-1,,drop=F]
  
  # (ii) 1-bp insertions
  ins.1bp <- x[nchar(as.character(x[,'REF'])) == nchar(as.character(x[,'ALT']))-1,,drop=F]
  
  # (iii) >1 bp deletions
  
  ## 2
  dels.2bp <- x[nchar(as.character(x[,'ALT'])) == nchar(as.character(x[,'REF']))-2,,drop=F]
  
  ## 3
  dels.3bp <- x[nchar(as.character(x[,'ALT'])) == nchar(as.character(x[,'REF']))-3,,drop=F]
  
  ## 4
  dels.4bp <- x[nchar(as.character(x[,'ALT'])) == nchar(as.character(x[,'REF']))-4,,drop=F]
  
  ## 5+
  dels.5bp <- x[nchar(as.character(x[,'ALT'])) <= nchar(as.character(x[,'REF']))-5,,drop=F]
  
  # (iv) >1 bp insertions
  
  ## 2
  ins.2bp <- x[nchar(as.character(x[,'REF'])) == nchar(as.character(x[,'ALT']))-2,,drop=F]
  
  ## 3
  ins.3bp <- x[nchar(as.character(x[,'REF'])) == nchar(as.character(x[,'ALT']))-3,,drop=F]
  
  ## 4
  ins.4bp <- x[nchar(as.character(x[,'REF'])) == nchar(as.character(x[,'ALT']))-4,,drop=F]
  
  ## 5+
  ins.5bp <- x[nchar(as.character(x[,'REF'])) <= nchar(as.character(x[,'ALT']))-5,,drop=F]
  
  ## 2. classify 1 bp events into:
  
  # (i) 1 bp deletions at homopolymers (length 1 == "no neighbouring homopolymers")
  if(nrow(dels.1bp) > 0){
    
    ## extract 10 bp upstream/downstream sequence context from reference
    dels.1bp.context <- as.character(subseq(x = reference[as.character(dels.1bp[,'CHROM'])], 
                                            start = as.numeric(dels.1bp[,'POS']) - 9, 
                                            end = as.numeric(dels.1bp[,'POS']) + 11))
    dels.1bp.context.middle <- paste0('[', str_split_fixed(dels.1bp.context, '', 21)[,11,drop=F], ']')
    dels.1bp.context.start <- str_split_fixed(dels.1bp.context, '', 11)[,1:10,drop=F]
    dels.1bp.context.start <- paste(dels.1bp.context.start[,1], dels.1bp.context.start[,2], dels.1bp.context.start[,3],
                                    dels.1bp.context.start[,4], dels.1bp.context.start[,5], dels.1bp.context.start[,6],
                                    dels.1bp.context.start[,7], dels.1bp.context.start[,8], dels.1bp.context.start[,9],
                                    dels.1bp.context.start[,10], sep = '')
    dels.1bp.context.end <- str_split_fixed(dels.1bp.context, '', 12)[,12,drop=F]
    dels.1bp[,'TRIPLET'] <- paste(dels.1bp.context.start, dels.1bp.context.middle, dels.1bp.context.end, sep = '')
    colnames(dels.1bp)[5] <- 'CONTEXT FW' 
    
    ## need the central base to be pyrimidine-centred, i.e. C or T, hence also run reverseComplements on full contexts
    dels.1bp.context.rc <- as.character(reverseComplement(subseq(x = reference[as.character(dels.1bp[,'CHROM'])], 
                                                                 start = as.numeric(dels.1bp[,'POS']) - 9, 
                                                                 end = as.numeric(dels.1bp[,'POS']) + 11)))
    dels.1bp.context.rc.middle <- paste0('[', str_split_fixed(dels.1bp.context.rc, '', 21)[,11,drop=F], ']')
    dels.1bp.context.rc.start <- str_split_fixed(dels.1bp.context.rc, '', 11)[,1:10,drop=F]
    dels.1bp.context.rc.start <- paste(dels.1bp.context.rc.start[,1], dels.1bp.context.rc.start[,2], dels.1bp.context.rc.start[,3],
                                       dels.1bp.context.rc.start[,4], dels.1bp.context.rc.start[,5], dels.1bp.context.rc.start[,6],
                                       dels.1bp.context.rc.start[,7], dels.1bp.context.rc.start[,8], dels.1bp.context.rc.start[,9],
                                       dels.1bp.context.rc.start[,10], sep = '')
    dels.1bp.context.rc.end <- str_split_fixed(dels.1bp.context.rc, '', 12)[,12,drop=F]
    dels.1bp <- cbind(dels.1bp[,1:5,drop=F], paste(dels.1bp.context.rc.start, dels.1bp.context.rc.middle, dels.1bp.context.rc.end, sep = ''))
    colnames(dels.1bp)[6] <- 'CONTEXT RC'
    
    ## summarise 1 bp deletions in matrix format
    dels.1bp.summary <- matrix(0, ncol = 2, nrow = 6)
    colnames(dels.1bp.summary) <- c('C', 'T') ## pyrimidine-centred deleted base
    rownames(dels.1bp.summary) <- c('0 bp', '1 bp', '2 bp', '3 bp', '4 bp', '5+ bp') ## contextual homopolymer-length
    fw <- str_split_fixed(str_split_fixed(dels.1bp[,'CONTEXT FW'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
    rc <- str_split_fixed(str_split_fixed(dels.1bp[,'CONTEXT RC'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
    fw.if <- fw == 'C' | fw == 'T'
    for (i in 1:nrow(dels.1bp)){
      
      if(fw.if[i] == T){
        
        ## look at forward context
        upstream.tmp <- str_split_fixed(str_split_fixed(as.character(dels.1bp[i,'CONTEXT FW']), '\\[', 2)[,1], '', 10)
        upstream.tmp <- which(upstream.tmp == fw[i])
        
        ### check all 6 categories for upstream bases
        if(any(upstream.tmp %in% 10)){
          
          if(any(upstream.tmp %in% 9)){
            
            if(any(upstream.tmp %in% 8)){
              
              if(any(upstream.tmp %in% 7)){
                
                if(any(upstream.tmp %in% 6)){
                  
                  upstream.tmp <- 6
                  
                }else{
                  upstream.tmp <- 5
                }
                
              }else{
                upstream.tmp <- 4
              }
              
            }else{
              upstream.tmp <- 3
            }
            
          }else{
            upstream.tmp <- 2
          }
          
        }else{
          upstream.tmp <- 1
        }
        
        downstream.tmp <- str_split_fixed(str_split_fixed(as.character(dels.1bp[i,'CONTEXT FW']), '\\]', 2)[,2], '', 10)
        downstream.tmp <- which(downstream.tmp == fw[i])
        
        ### check all 6 categories for downstream bases
        if(any(downstream.tmp %in% 1)){
          
          if(any(downstream.tmp %in% 2)){
            
            if(any(downstream.tmp %in% 3)){
              
              if(any(downstream.tmp %in% 4)){
                
                if(any(downstream.tmp %in% 5)){
                  
                  downstream.tmp <- 6
                  
                }else{
                  downstream.tmp <- 5
                }
                
              }else{
                downstream.tmp <- 4
              }
              
            }else{
              downstream.tmp <- 3
            }
            
          }else{
            downstream.tmp <- 2
          }
          
        }else{
          downstream.tmp <- 1
        }
        
        ## summarise, which homopolymer (upstream vs. downstream context) is longer
        ### MODIFIED (ao7) - handle 'N' bases
        if (fw[i] %in% c("A","C","G","T")) {
          dels.1bp.summary[max(c(upstream.tmp, downstream.tmp)),fw[i]] <- dels.1bp.summary[max(c(upstream.tmp, downstream.tmp)),fw[i]] + 1
        }
        ###
        
      } else {
        
        ## look at reverse complement context
        upstream.tmp <- str_split_fixed(str_split_fixed(as.character(dels.1bp[i,'CONTEXT RC']), '\\[', 2)[,1], '', 10)
        upstream.tmp <- which(upstream.tmp == rc[i])
        
        ### check all 6 categories for upstream bases
        if(any(upstream.tmp %in% 10)){
          
          if(any(upstream.tmp %in% 9)){
            
            if(any(upstream.tmp %in% 8)){
              
              if(any(upstream.tmp %in% 7)){
                
                if(any(upstream.tmp %in% 6)){
                  
                  upstream.tmp <- 6
                  
                }else{
                  upstream.tmp <- 5
                }
                
              }else{
                upstream.tmp <- 4
              }
              
            }else{
              upstream.tmp <- 3
            }
            
          }else{
            upstream.tmp <- 2
          }
          
        }else{
          upstream.tmp <- 1
        }
        
        downstream.tmp <- str_split_fixed(str_split_fixed(as.character(dels.1bp[i,'CONTEXT RC']), '\\]', 2)[,2], '', 10)
        downstream.tmp <- which(downstream.tmp == rc[i])
        
        ### check all 6 categories for downstream bases
        if(any(downstream.tmp %in% 1)){
          
          if(any(downstream.tmp %in% 2)){
            
            if(any(downstream.tmp %in% 3)){
              
              if(any(downstream.tmp %in% 4)){
                
                if(any(downstream.tmp %in% 5)){
                  
                  downstream.tmp <- 6
                  
                }else{
                  downstream.tmp <- 5
                }
                
              }else{
                downstream.tmp <- 4
              }
              
            }else{
              downstream.tmp <- 3
            }
            
          }else{
            downstream.tmp <- 2
          }
          
        }else{
          downstream.tmp <- 1
        }
        
        ## summarise, which homopolymer (upstream vs. downstream context) is longer
        ### MODIFIED (ao7) - handle 'N' bases
        if (rc[i] %in% c("A","C","G","T")) {
          dels.1bp.summary[max(c(upstream.tmp, downstream.tmp)),rc[i]] <- dels.1bp.summary[max(c(upstream.tmp, downstream.tmp)),rc[i]] + 1
        }
        ###
      }
      
    }
    
  }else{

    ## summarise 1 bp deletions in matrix format
    dels.1bp.summary <- matrix(0, ncol = 2, nrow = 6)
    colnames(dels.1bp.summary) <- c('C', 'T') ## pyrimidine-centred deleted base
    rownames(dels.1bp.summary) <- c('0 bp', '1 bp', '2 bp', '3 bp', '4 bp', '5+ bp') ## contextual homopolymer-length
    
  }
  
  # (ii) 1 bp insertions at homopolymers (length 0 == "no neighbouring homopolymers")
  if(nrow(ins.1bp) > 0){
    
    ## extract 10 bp upstream/downstream sequence context from reference
    ins.1bp.context <- as.character(subseq(x = reference[as.character(ins.1bp[,'CHROM'])], 
                                           start = as.numeric(ins.1bp[,'POS']) - 9, 
                                           end = as.numeric(ins.1bp[,'POS']) + 10))
    ## insertion after base pos. 10
    ins.1bp.context.middle <- paste0('[', str_split_fixed(as.character(ins.1bp[,'ALT']), '', 2)[,2,drop=F], ']')
    ins.1bp.context.start <- str_split_fixed(ins.1bp.context, '', 11)[,1:10,drop=F]
    ins.1bp.context.start <- paste(ins.1bp.context.start[,1], ins.1bp.context.start[,2], ins.1bp.context.start[,3],
                                   ins.1bp.context.start[,4], ins.1bp.context.start[,5], ins.1bp.context.start[,6],
                                   ins.1bp.context.start[,7], ins.1bp.context.start[,8], ins.1bp.context.start[,9],
                                   ins.1bp.context.start[,10], sep = '')
    ins.1bp.context.end <- str_split_fixed(ins.1bp.context, '', 11)[,11,drop=F]
    ins.1bp[,'TRIPLET'] <- paste(ins.1bp.context.start, ins.1bp.context.middle, ins.1bp.context.end, sep = '')
    colnames(ins.1bp)[5] <- 'CONTEXT FW'
    
    ## need the central base to be pyrimidine-centred, i.e. C or T, hence also run reverseComplements on full contexts
    ins.1bp.context.rc <- as.character(reverseComplement(subseq(x = reference[as.character(ins.1bp[,'CHROM'])], 
                                                                start = as.numeric(ins.1bp[,'POS']) - 9, 
                                                                end = as.numeric(ins.1bp[,'POS']) + 10)))
    ins.1bp.context.rc.middle <- paste0('[', as.character(reverseComplement(DNAStringSet(str_split_fixed(as.character(ins.1bp[,'ALT']), '', 2)[,2,drop=F]))), ']')
    ins.1bp.context.rc.start <- str_split_fixed(ins.1bp.context.rc, '', 11)[,1:10,drop=F]
    ins.1bp.context.rc.start <- paste(ins.1bp.context.rc.start[,1], ins.1bp.context.rc.start[,2], ins.1bp.context.rc.start[,3],
                                      ins.1bp.context.rc.start[,4], ins.1bp.context.rc.start[,5], ins.1bp.context.rc.start[,6],
                                      ins.1bp.context.rc.start[,7], ins.1bp.context.rc.start[,8], ins.1bp.context.rc.start[,9],
                                      ins.1bp.context.rc.start[,10], sep = '')
    ins.1bp.context.rc.end <- str_split_fixed(ins.1bp.context.rc, '', 11)[,11,drop=F]
    ins.1bp <- cbind(ins.1bp[,1:5,drop=F], paste(ins.1bp.context.rc.start, ins.1bp.context.rc.middle, ins.1bp.context.rc.end, sep = ''))
    colnames(ins.1bp)[6] <- 'CONTEXT RC'
    
    ## summarise 1 bp insertions in matrix format
    ins.1bp.summary <- matrix(0, ncol = 2, nrow = 6)
    colnames(ins.1bp.summary) <- c('C', 'T') ## pyrimidine-centred inserted base
    rownames(ins.1bp.summary) <- c('0 bp', '1 bp', '2 bp', '3 bp', '4 bp', '5+ bp') ## contextual homopolymer-length
    fw <- str_split_fixed(str_split_fixed(ins.1bp[,'CONTEXT FW'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
    rc <- str_split_fixed(str_split_fixed(ins.1bp[,'CONTEXT RC'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
    fw.if <- fw == 'C' | fw == 'T'
    for (i in 1:nrow(ins.1bp)){
      
      if(fw.if[i] == T){
        
        ## look at forward context
        upstream.tmp <- str_split_fixed(str_split_fixed(as.character(ins.1bp[i,'CONTEXT FW']), '\\[', 2)[,1], '', 10)
        upstream.tmp <- which(upstream.tmp == fw[i])
        
        ### check all 6 categories for upstream bases
        if(any(upstream.tmp %in% 10)){
          
          if(any(upstream.tmp %in% 9)){
            
            if(any(upstream.tmp %in% 8)){
              
              if(any(upstream.tmp %in% 7)){
                
                if(any(upstream.tmp %in% 6)){
                  
                  upstream.tmp <- 6
                  
                }else{
                  upstream.tmp <- 5
                }
                
              }else{
                upstream.tmp <- 4
              }
              
            }else{
              upstream.tmp <- 3
            }
            
          }else{
            upstream.tmp <- 2
          }
          
        }else{
          upstream.tmp <- 1
        }
        
        downstream.tmp <- str_split_fixed(str_split_fixed(as.character(ins.1bp[i,'CONTEXT FW']), '\\]', 2)[,2], '', 10)
        downstream.tmp <- which(downstream.tmp == fw[i])
        
        ### check all 6 categories for downstream bases
        if(any(downstream.tmp %in% 1)){
          
          if(any(downstream.tmp %in% 2)){
            
            if(any(downstream.tmp %in% 3)){
              
              if(any(downstream.tmp %in% 4)){
                
                if(any(downstream.tmp %in% 5)){
                  
                  downstream.tmp <- 6
                  
                }else{
                  downstream.tmp <- 5
                }
                
              }else{
                downstream.tmp <- 4
              }
              
            }else{
              downstream.tmp <- 3
            }
            
          }else{
            downstream.tmp <- 2
          }
          
        }else{
          downstream.tmp <- 1
        }
        
        ## summarise, which homopolymer is longer
        ### MODIFIED (ao7) - handle 'N' bases
        if (fw[i] %in% c("A","C","G","T")) {
          ins.1bp.summary[max(c(upstream.tmp, downstream.tmp)),fw[i]] <- ins.1bp.summary[max(c(upstream.tmp, downstream.tmp)),fw[i]] + 1
        }
        ###
        
      } else {
        
        ## look at reverse complement context
        upstream.tmp <- str_split_fixed(str_split_fixed(as.character(ins.1bp[i,'CONTEXT RC']), '\\[', 2)[,1], '', 10)
        upstream.tmp <- which(upstream.tmp == rc[i])
        
        ### check all 6 categories for upstream bases
        if(any(upstream.tmp %in% 10)){
          
          if(any(upstream.tmp %in% 9)){
            
            if(any(upstream.tmp %in% 8)){
              
              if(any(upstream.tmp %in% 7)){
                
                if(any(upstream.tmp %in% 6)){
                  
                  upstream.tmp <- 6
                  
                }else{
                  upstream.tmp <- 5
                }
                
              }else{
                upstream.tmp <- 4
              }
              
            }else{
              upstream.tmp <- 3
            }
            
          }else{
            upstream.tmp <- 2
          }
          
        }else{
          upstream.tmp <- 1
        }
        
        downstream.tmp <- str_split_fixed(str_split_fixed(as.character(ins.1bp[i,'CONTEXT RC']), '\\]', 2)[,2], '', 10)
        downstream.tmp <- which(downstream.tmp == rc[i])
        
        ### check all 6 categories for downstream bases
        if(any(downstream.tmp %in% 1)){
          
          if(any(downstream.tmp %in% 2)){
            
            if(any(downstream.tmp %in% 3)){
              
              if(any(downstream.tmp %in% 4)){
                
                if(any(downstream.tmp %in% 5)){
                  
                  downstream.tmp <- 6
                  
                }else{
                  downstream.tmp <- 5
                }
                
              }else{
                downstream.tmp <- 4
              }
              
            }else{
              downstream.tmp <- 3
            }
            
          }else{
            downstream.tmp <- 2
          }
          
        }else{
          downstream.tmp <- 1
        }
        
        ## summarise, which homopolymer is longer
        ### MODIFIED (ao7) - handle 'N' bases
        if (rc[i] %in% c("A","C","G","T")) {
          ins.1bp.summary[max(c(upstream.tmp, downstream.tmp)),rc[i]] <- ins.1bp.summary[max(c(upstream.tmp, downstream.tmp)),rc[i]] + 1
        }
        ###
        
      }
      
    } 
    
  }else{
    
    ## summarise 1 bp insertions in matrix format
    ins.1bp.summary <- matrix(0, ncol = 2, nrow = 6)
    colnames(ins.1bp.summary) <- c('C', 'T') ## pyrimidine-centred inserted base
    rownames(ins.1bp.summary) <- c('0 bp', '1 bp', '2 bp', '3 bp', '4 bp', '5+ bp') ## contextual homopolymer-length
    
  }
  
  ## 3. classify >=2 bp deletions into:
  
  # (i) 2 bp deletions at simple repeats (length 1 == "no neighbouring simple repeat")
  if(nrow(dels.2bp) > 0){
    
    ## extract 5 x 2 bp upstream/downstream sequence context from reference
    dels.2bp.context <- as.character(subseq(x = reference[as.character(dels.2bp[,'CHROM'])], 
                                            start = as.numeric(dels.2bp[,'POS']) - 9, 
                                            end = as.numeric(dels.2bp[,'POS']) + 2 + 10))
    dels.2bp.context.middle <- str_split_fixed(dels.2bp.context, '', 22)[,11:12,drop=F]
    dels.2bp.context.middle <- paste0('[', dels.2bp.context.middle[,1], dels.2bp.context.middle[,2], ']')
    dels.2bp.context.start <- str_split_fixed(dels.2bp.context, '', 11)[,1:10,drop=F]
    dels.2bp.context.start <- paste(dels.2bp.context.start[,1], dels.2bp.context.start[,2], dels.2bp.context.start[,3],
                                    dels.2bp.context.start[,4], dels.2bp.context.start[,5], dels.2bp.context.start[,6],
                                    dels.2bp.context.start[,7], dels.2bp.context.start[,8], dels.2bp.context.start[,9],
                                    dels.2bp.context.start[,10], sep = '')
    dels.2bp.context.end <- str_split_fixed(dels.2bp.context, '', 13)[,13,drop=F]
    dels.2bp[,'TRIPLET'] <- paste(dels.2bp.context.start, dels.2bp.context.middle, dels.2bp.context.end, sep = '')
    colnames(dels.2bp)[5] <- 'CONTEXT' 
    
  }
  
  # (ii) 3 bp deletions at simple repeats (length 1 == "no neighbouring simple repeat")
  if(nrow(dels.3bp) > 0){
   
    ## extract 5 x 3 bp upstream/downstream sequence context from reference
    dels.3bp.context <- as.character(subseq(x = reference[as.character(dels.3bp[,'CHROM'])], 
                                            start = as.numeric(dels.3bp[,'POS']) - 14, 
                                            end = as.numeric(dels.3bp[,'POS']) + 3 + 15))
    dels.3bp.context.middle <- str_split_fixed(dels.3bp.context, '', 33)[,16:18,drop=F]
    dels.3bp.context.middle <- paste0('[', dels.3bp.context.middle[,1], dels.3bp.context.middle[,2], dels.3bp.context.middle[,3], ']')
    dels.3bp.context.start <- str_split_fixed(dels.3bp.context, '', 33)[,1:15,drop=F]
    dels.3bp.context.start <- paste(dels.3bp.context.start[,1], dels.3bp.context.start[,2], dels.3bp.context.start[,3],
                                    dels.3bp.context.start[,4], dels.3bp.context.start[,5], dels.3bp.context.start[,6],
                                    dels.3bp.context.start[,7], dels.3bp.context.start[,8], dels.3bp.context.start[,9],
                                    dels.3bp.context.start[,10], dels.3bp.context.start[,11], dels.3bp.context.start[,12], 
                                    dels.3bp.context.start[,13], dels.3bp.context.start[,14], dels.3bp.context.start[,15], sep = '')
    dels.3bp.context.end <- str_split_fixed(dels.3bp.context, '', 19)[,19,drop=F]
    dels.3bp[,'TRIPLET'] <- paste(dels.3bp.context.start, dels.3bp.context.middle, dels.3bp.context.end, sep = '')
    colnames(dels.3bp)[5] <- 'CONTEXT' 
    
  }
  
  # (iii) 4 bp deletions at simple repeats (length 1 == "no neighbouring simple repeat")
  if(nrow(dels.4bp) > 0){
    
    ## extract 5 x 4 bp upstream/downstream sequence context from reference
    dels.4bp.context <- as.character(subseq(x = reference[as.character(dels.4bp[,'CHROM'])], 
                                            start = as.numeric(dels.4bp[,'POS']) - 19, 
                                            end = as.numeric(dels.4bp[,'POS']) + 4 + 20))
    dels.4bp.context.middle <- str_split_fixed(dels.4bp.context, '', 41)[,21:24,drop=F]
    dels.4bp.context.middle <- paste0('[', dels.4bp.context.middle[,1], 
                                      dels.4bp.context.middle[,2], 
                                      dels.4bp.context.middle[,3], 
                                      dels.4bp.context.middle[,4], ']')
    dels.4bp.context.start <- str_split_fixed(dels.4bp.context, '', 21)[,1:20,drop=F]
    dels.4bp.context.start <- paste(dels.4bp.context.start[,1], dels.4bp.context.start[,2], dels.4bp.context.start[,3],
                                    dels.4bp.context.start[,4], dels.4bp.context.start[,5], dels.4bp.context.start[,6],
                                    dels.4bp.context.start[,7], dels.4bp.context.start[,8], dels.4bp.context.start[,9],
                                    dels.4bp.context.start[,10], dels.4bp.context.start[,11], dels.4bp.context.start[,12], 
                                    dels.4bp.context.start[,13], dels.4bp.context.start[,14], dels.4bp.context.start[,15], 
                                    dels.4bp.context.start[,16], dels.4bp.context.start[,17], dels.4bp.context.start[,18], 
                                    dels.4bp.context.start[,19], dels.4bp.context.start[,20], sep = '')
    dels.4bp.context.end <- str_split_fixed(dels.4bp.context, '', 25)[,25,drop=F]
    dels.4bp[,'TRIPLET'] <- paste(dels.4bp.context.start, dels.4bp.context.middle, dels.4bp.context.end, sep = '')
    colnames(dels.4bp)[5] <- 'CONTEXT' 
    
  }
  
  # (iv) 5+ bp deletions at simple repeats (length 1 == "no neighbouring simple repeat")
  if(nrow(dels.5bp) > 0){
    
    ## extract 5 x 100 bp upstream/downstream sequence context from reference
    ## i.e. max. 100 bp repeat motif
    dels.5bp.context <- as.character(subseq(x = reference[as.character(dels.5bp[,'CHROM'])], 
                                            start = as.numeric(dels.5bp[,'POS']) - 499, 
                                            end = as.numeric(dels.5bp[,'POS']) + 5 + 500))
    
    ### account for the different repeat lengths: iterate
    dels.5bp.context.middle <- as.character(dels.5bp[,'REF'])
    dels.5bp.context.middle <- str_split_fixed(dels.5bp.context.middle, '', 2)[,2,drop=F]
    dels.5bp.context.middle <- paste0('[', dels.5bp.context.middle, ']')
    dels.5bp.context.middle.lengths <- nchar(as.character(dels.5bp[,'REF'])) - 1
    dels.5bp.context.start <- rep(NA, nrow(dels.5bp))
    for (i in 1:length(dels.5bp.context.start)){
      tmp.dels.5bp.context.start <- str_split_fixed(dels.5bp.context[i], '', 501)[,c(1 + 500-c(5*dels.5bp.context.middle.lengths[i])):500,drop=F]
      dels.5bp.context.start[i] <- paste(tmp.dels.5bp.context.start, collapse = '')
    }
    dels.5bp.context.end <- rep(NA, nrow(dels.5bp))
    for (i in 1:length(dels.5bp.context.end)){
      tmp.dels.5bp.context.end <- str_split_fixed(dels.5bp.context[i], '', 1006)[,c(501+dels.5bp.context.middle.lengths[i]):c(501 + dels.5bp.context.middle.lengths[i]*6 - 1),drop=F]
      dels.5bp.context.end[i] <- paste(tmp.dels.5bp.context.end, collapse = '')
    }
    dels.5bp[,'TRIPLET'] <- paste(dels.5bp.context.start, dels.5bp.context.middle, dels.5bp.context.end, sep = '')
    colnames(dels.5bp)[5] <- 'CONTEXT'
    
  }
  
  ## summarise >=2 bp deletions at simple repeats in matrix format
  ## also create sub-tables for repeat length = 1 hits, to later assess these for microhomologies
  dels.greater.2bp.summary <- matrix(0, ncol = 4, nrow = 6)
  colnames(dels.greater.2bp.summary) <- c('2 bp', '3 bp', '4 bp', '5+ bp') ## deletion size
  rownames(dels.greater.2bp.summary) <- c('1 or MH', '2', '3', '4', '5', '6+') ## number of repeats
  
  ## 2 bp
  if(nrow(dels.2bp) > 0){
    
    repeat.nts <- str_split_fixed(str_split_fixed(dels.2bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
    downstream.context <- str_split_fixed(str_split_fixed(dels.2bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
    dels.2bp.pot.MH.ind <- c()
    for (i in 1:nrow(dels.2bp)){
      
      ## look at repeat
      tmp.repeat.nts <- repeat.nts[i]
      tmp.repeat.length <- nchar(tmp.repeat.nts)
      
      ## how often does it match consecutively in the immediate downstream context?
      tmp.downstream.context <- downstream.context[i]
      tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]]
      
      ## group
      tmp.downstream.context <- c(paste(tmp.downstream.context[c(tmp.repeat.length-c(tmp.repeat.length - 1)):tmp.repeat.length], collapse = ''),
                                  paste(tmp.downstream.context[c(2*tmp.repeat.length-c(tmp.repeat.length - 1)):c(2*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(3*tmp.repeat.length-c(tmp.repeat.length - 1)):c(3*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(4*tmp.repeat.length-c(tmp.repeat.length - 1)):c(4*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(5*tmp.repeat.length-c(tmp.repeat.length - 1)):c(5*tmp.repeat.length)], collapse = ''))
      
      ### check all 6 categories for upstream bases
      if(tmp.downstream.context[1] %in% tmp.repeat.nts){
        
        if(tmp.downstream.context[2] %in% tmp.repeat.nts){
          
          if(tmp.downstream.context[3] %in% tmp.repeat.nts){
            
            if(tmp.downstream.context[4] %in% tmp.repeat.nts){
              
              if(tmp.downstream.context[5] %in% tmp.repeat.nts){
                
                dels.greater.2bp.summary['6+', '2 bp'] <- dels.greater.2bp.summary['6+', '2 bp'] + 1
                
              } else{
                
                dels.greater.2bp.summary['5', '2 bp'] <- dels.greater.2bp.summary['5', '2 bp'] + 1
                
              }
              
            }else{
              
              dels.greater.2bp.summary['4', '2 bp'] <- dels.greater.2bp.summary['4', '2 bp'] + 1
              
            }
            
          }else{
            
            dels.greater.2bp.summary['3', '2 bp'] <- dels.greater.2bp.summary['3', '2 bp'] + 1
            
          }
          
        }else{
          
          dels.greater.2bp.summary['2', '2 bp'] <- dels.greater.2bp.summary['2', '2 bp'] + 1
          
        }
        
      }else{
        
        dels.greater.2bp.summary['1 or MH', '2 bp'] <- dels.greater.2bp.summary['1 or MH', '2 bp'] + 1
        dels.2bp.pot.MH.ind <- c(dels.2bp.pot.MH.ind, i)
        
      }
      
    }
    dels.2bp.pot.MH <- dels.2bp[dels.2bp.pot.MH.ind,,drop=F] 
    
  }
  
  ## 3 bp
  if(nrow(dels.3bp) > 0){
    
    repeat.nts <- str_split_fixed(str_split_fixed(dels.3bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
    downstream.context <- str_split_fixed(str_split_fixed(dels.3bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
    dels.3bp.pot.MH.ind <- c()
    for (i in 1:nrow(dels.3bp)){
      
      ## look at repeat
      tmp.repeat.nts <- repeat.nts[i]
      tmp.repeat.length <- nchar(tmp.repeat.nts)
      
      ## how often does it match consecutively in the immediate downstream context?
      tmp.downstream.context <- downstream.context[i]
      tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]]
      
      ## group
      tmp.downstream.context <- c(paste(tmp.downstream.context[c(tmp.repeat.length-c(tmp.repeat.length - 1)):tmp.repeat.length], collapse = ''),
                                  paste(tmp.downstream.context[c(2*tmp.repeat.length-c(tmp.repeat.length - 1)):c(2*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(3*tmp.repeat.length-c(tmp.repeat.length - 1)):c(3*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(4*tmp.repeat.length-c(tmp.repeat.length - 1)):c(4*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(5*tmp.repeat.length-c(tmp.repeat.length - 1)):c(5*tmp.repeat.length)], collapse = ''))
      
      ### check all 6 categories for upstream bases
      if(tmp.downstream.context[1] %in% tmp.repeat.nts){
        
        if(tmp.downstream.context[2] %in% tmp.repeat.nts){
          
          if(tmp.downstream.context[3] %in% tmp.repeat.nts){
            
            if(tmp.downstream.context[4] %in% tmp.repeat.nts){
              
              if(tmp.downstream.context[5] %in% tmp.repeat.nts){
                
                dels.greater.2bp.summary['6+', '3 bp'] <- dels.greater.2bp.summary['6+', '3 bp'] + 1
                
              } else{
                
                dels.greater.2bp.summary['5', '3 bp'] <- dels.greater.2bp.summary['5', '3 bp'] + 1
                
              }
              
            }else{
              
              dels.greater.2bp.summary['4', '3 bp'] <- dels.greater.2bp.summary['4', '3 bp'] + 1
              
            }
            
          }else{
            
            dels.greater.2bp.summary['3', '3 bp'] <- dels.greater.2bp.summary['3', '3 bp'] + 1
            
          }
          
        }else{
          
          dels.greater.2bp.summary['2', '3 bp'] <- dels.greater.2bp.summary['2', '3 bp'] + 1
          
        }
        
      }else{
        
        dels.greater.2bp.summary['1 or MH', '3 bp'] <- dels.greater.2bp.summary['1 or MH', '3 bp'] + 1
        dels.3bp.pot.MH.ind <- c(dels.3bp.pot.MH.ind, i)
        
      }
      
    }
    dels.3bp.pot.MH <- dels.3bp[dels.3bp.pot.MH.ind,,drop=F] 
    
  }
  
  ## 4 bp
  if(nrow(dels.4bp) > 0){
    
    repeat.nts <- str_split_fixed(str_split_fixed(dels.4bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
    downstream.context <- str_split_fixed(str_split_fixed(dels.4bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
    dels.4bp.pot.MH.ind <- c()
    for (i in 1:nrow(dels.4bp)){
      
      ## look at repeat
      tmp.repeat.nts <- repeat.nts[i]
      tmp.repeat.length <- nchar(tmp.repeat.nts)
      
      ## how often does it match consecutively in the immediate downstream context?
      tmp.downstream.context <- downstream.context[i]
      tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]]
      
      ## group
      tmp.downstream.context <- c(paste(tmp.downstream.context[c(tmp.repeat.length-c(tmp.repeat.length - 1)):tmp.repeat.length], collapse = ''),
                                  paste(tmp.downstream.context[c(2*tmp.repeat.length-c(tmp.repeat.length - 1)):c(2*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(3*tmp.repeat.length-c(tmp.repeat.length - 1)):c(3*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(4*tmp.repeat.length-c(tmp.repeat.length - 1)):c(4*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(5*tmp.repeat.length-c(tmp.repeat.length - 1)):c(5*tmp.repeat.length)], collapse = ''))
      
      ### check all 6 categories for upstream bases
      if(tmp.downstream.context[1] %in% tmp.repeat.nts){
        
        if(tmp.downstream.context[2] %in% tmp.repeat.nts){
          
          if(tmp.downstream.context[3] %in% tmp.repeat.nts){
            
            if(tmp.downstream.context[4] %in% tmp.repeat.nts){
              
              if(tmp.downstream.context[5] %in% tmp.repeat.nts){
                
                dels.greater.2bp.summary['6+', '4 bp'] <- dels.greater.2bp.summary['6+', '4 bp'] + 1
                
              } else{
                
                dels.greater.2bp.summary['5', '4 bp'] <- dels.greater.2bp.summary['5', '4 bp'] + 1
                
              }
              
            }else{
              
              dels.greater.2bp.summary['4', '4 bp'] <- dels.greater.2bp.summary['4', '4 bp'] + 1
              
            }
            
          }else{
            
            dels.greater.2bp.summary['3', '4 bp'] <- dels.greater.2bp.summary['3', '4 bp'] + 1
            
          }
          
        }else{
          
          dels.greater.2bp.summary['2', '4 bp'] <- dels.greater.2bp.summary['2', '4 bp'] + 1
          
        }
        
      }else{
        
        dels.greater.2bp.summary['1 or MH', '4 bp'] <- dels.greater.2bp.summary['1 or MH', '4 bp'] + 1
        dels.4bp.pot.MH.ind <- c(dels.4bp.pot.MH.ind, i)
        
      }
      
    }
    dels.4bp.pot.MH <- dels.4bp[dels.4bp.pot.MH.ind,,drop=F]
    
  }
  
  ## 5+ bp
  if(nrow(dels.5bp) > 0){
    
    repeat.nts <- str_split_fixed(str_split_fixed(dels.5bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
    downstream.context <- str_split_fixed(str_split_fixed(dels.5bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
    dels.5bp.pot.MH.ind <- c()
    for (i in 1:nrow(dels.5bp)){
      
      ## look at repeat
      tmp.repeat.nts <- repeat.nts[i]
      tmp.repeat.length <- nchar(tmp.repeat.nts)
      
      ## how often does it match consecutively in the immediate downstream context?
      tmp.downstream.context <- downstream.context[i]
      tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]]
      
      ## group
      tmp.downstream.context <- c(paste(tmp.downstream.context[c(tmp.repeat.length-c(tmp.repeat.length - 1)):tmp.repeat.length], collapse = ''),
                                  paste(tmp.downstream.context[c(2*tmp.repeat.length-c(tmp.repeat.length - 1)):c(2*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(3*tmp.repeat.length-c(tmp.repeat.length - 1)):c(3*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(4*tmp.repeat.length-c(tmp.repeat.length - 1)):c(4*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(5*tmp.repeat.length-c(tmp.repeat.length - 1)):c(5*tmp.repeat.length)], collapse = ''))
      
      ### check all 6 categories for upstream bases
      if(tmp.downstream.context[1] %in% tmp.repeat.nts){
        
        if(tmp.downstream.context[2] %in% tmp.repeat.nts){
          
          if(tmp.downstream.context[3] %in% tmp.repeat.nts){
            
            if(tmp.downstream.context[4] %in% tmp.repeat.nts){
              
              if(tmp.downstream.context[5] %in% tmp.repeat.nts){
                
                dels.greater.2bp.summary['6+', '5+ bp'] <- dels.greater.2bp.summary['6+', '5+ bp'] + 1
                
              } else{
                
                dels.greater.2bp.summary['5', '5+ bp'] <- dels.greater.2bp.summary['5', '5+ bp'] + 1
                
              }
              
            }else{
              
              dels.greater.2bp.summary['4', '5+ bp'] <- dels.greater.2bp.summary['4', '5+ bp'] + 1
              
            }
            
          }else{
            
            dels.greater.2bp.summary['3', '5+ bp'] <- dels.greater.2bp.summary['3', '5+ bp'] + 1
            
          }
          
        }else{
          
          dels.greater.2bp.summary['2', '5+ bp'] <- dels.greater.2bp.summary['2', '5+ bp'] + 1
          
        }
        
      }else{
        
        dels.greater.2bp.summary['1 or MH', '5+ bp'] <- dels.greater.2bp.summary['1 or MH', '5+ bp'] + 1
        dels.5bp.pot.MH.ind <- c(dels.5bp.pot.MH.ind, i)
        
      }
      
    }
    dels.5bp.pot.MH <- dels.5bp[dels.5bp.pot.MH.ind,,drop=F] 
    
  }
  
  # (v) >=2 bp deletions at microhomologies
  ## go directly into in matrix format, thereby also re-classify repeat length = 1 deletions
  dels.greater.2bp.MH.summary <- matrix(0, ncol = 4, nrow = 5)
  colnames(dels.greater.2bp.MH.summary) <- c('2 bp', '3 bp', '4 bp', '5+ bp') ## deletion size
  rownames(dels.greater.2bp.MH.summary) <- c('1 bp MH', '2 bp MH', '3 bp MH', '4 bp MH', '5+ bp MH') ## microhomology length
  dels.greater.2bp.MH.summary[2:5, '2 bp'] <- NA
  dels.greater.2bp.MH.summary[3:5, '3 bp'] <- NA
  dels.greater.2bp.MH.summary[4:5, '4 bp'] <- NA
  
  ## 2 bp
  if(nrow(dels.2bp) > 0){
    
    if(nrow(dels.2bp.pot.MH) > 0){
      
      repeat.nts <- str_split_fixed(str_split_fixed(dels.2bp.pot.MH[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
      upstream.context <- str_split_fixed(dels.2bp.pot.MH[,'CONTEXT'], '\\[', 2)[,1,drop=F]
      downstream.context <- str_split_fixed(str_split_fixed(dels.2bp.pot.MH[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
      for (i in 1:nrow(dels.2bp.pot.MH)){
        
        ## look at repeat
        tmp.repeat.nts <- repeat.nts[i]
        
        ## isolate first nucleotide in the immediate upstream context
        tmp.upstream.context <- upstream.context[i]
        tmp.upstream.context <- strsplit(tmp.upstream.context, '')[[1]][length(strsplit(tmp.upstream.context, '')[[1]])]
        
        ## isolate first nucleotide in the immediate downstream context
        tmp.downstream.context <- downstream.context[i]
        tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]][1]
        
        if(strsplit(tmp.repeat.nts, '')[[1]][2] == tmp.upstream.context | strsplit(tmp.repeat.nts, '')[[1]][1] == tmp.downstream.context){
          
          dels.greater.2bp.MH.summary['1 bp MH', '2 bp'] <- dels.greater.2bp.MH.summary['1 bp MH', '2 bp'] + 1
          dels.greater.2bp.summary['1 or MH', '2 bp'] <- dels.greater.2bp.summary['1 or MH', '2 bp'] - 1
          
        }
      } 
      
    }
    
  }
  
  ## 3 bp
  if(nrow(dels.3bp) > 0){
    
    if(nrow(dels.3bp.pot.MH) > 0){
      
      repeat.nts <- str_split_fixed(str_split_fixed(dels.3bp.pot.MH[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
      upstream.context <- str_split_fixed(dels.3bp.pot.MH[,'CONTEXT'], '\\[', 2)[,1,drop=F]
      downstream.context <- str_split_fixed(str_split_fixed(dels.3bp.pot.MH[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
      for (i in 1:nrow(dels.3bp.pot.MH)){
        
        ## look at repeat
        tmp.repeat.nts <- repeat.nts[i]
        
        ## isolate first two nucleotides in the immediate upstream context
        tmp.upstream.context <- upstream.context[i]
        tmp.upstream.context <- strsplit(tmp.upstream.context, '')[[1]][c(length(strsplit(tmp.upstream.context, '')[[1]]) - 1):length(strsplit(tmp.upstream.context, '')[[1]])]
        
        ## isolate first two nucleotides in the immediate downstream context
        tmp.downstream.context <- downstream.context[i]
        tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]][1:2]
        
        ### MH length = 2 (important to start with the highest possible MH length)
        
        if(all(c(strsplit(tmp.repeat.nts, '')[[1]][2:3] == tmp.upstream.context) == T) | all(c(strsplit(tmp.repeat.nts, '')[[1]][1:2] == tmp.downstream.context) == T)){
          
          dels.greater.2bp.MH.summary['2 bp MH', '3 bp'] <- dels.greater.2bp.MH.summary['2 bp MH', '3 bp'] + 1
          dels.greater.2bp.summary['1 or MH', '3 bp'] <- dels.greater.2bp.summary['1 or MH', '3 bp'] - 1
          
        }else{
          
          ### MH length = 1
          
          if(strsplit(tmp.repeat.nts, '')[[1]][3] == tmp.upstream.context[2] | strsplit(tmp.repeat.nts, '')[[1]][1] == tmp.downstream.context[1]){
            
            dels.greater.2bp.MH.summary['1 bp MH', '3 bp'] <- dels.greater.2bp.MH.summary['1 bp MH', '3 bp'] + 1
            dels.greater.2bp.summary['1 or MH', '3 bp'] <- dels.greater.2bp.summary['1 or MH', '3 bp'] - 1
            
          }
          
        }
        
      } 
      
    }
      
  }
  
  ## 4 bp
  if(nrow(dels.4bp) > 0){
    
    if(nrow(dels.4bp.pot.MH) > 0){
      
      repeat.nts <- str_split_fixed(str_split_fixed(dels.4bp.pot.MH[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
      upstream.context <- str_split_fixed(dels.4bp.pot.MH[,'CONTEXT'], '\\[', 2)[,1,drop=F]
      downstream.context <- str_split_fixed(str_split_fixed(dels.4bp.pot.MH[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
      for (i in 1:nrow(dels.4bp.pot.MH)){
        
        ## look at repeat
        tmp.repeat.nts <- repeat.nts[i]
        
        ## isolate first three nucleotides in the immediate upstream context
        tmp.upstream.context <- upstream.context[i]
        tmp.upstream.context <- strsplit(tmp.upstream.context, '')[[1]][c(length(strsplit(tmp.upstream.context, '')[[1]]) - 2):length(strsplit(tmp.upstream.context, '')[[1]])]
        
        ## isolate first three nucleotides in the immediate downstream context
        tmp.downstream.context <- downstream.context[i]
        tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]][1:3]
        
        ### MH length = 3 (important to start with the highest possible MH length)
        
        if(all(c(strsplit(tmp.repeat.nts, '')[[1]][2:4] == tmp.upstream.context) == T) | all(c(strsplit(tmp.repeat.nts, '')[[1]][1:3] == tmp.downstream.context) == T)){
          
          dels.greater.2bp.MH.summary['3 bp MH', '4 bp'] <- dels.greater.2bp.MH.summary['3 bp MH', '4 bp'] + 1
          dels.greater.2bp.summary['1 or MH', '4 bp'] <- dels.greater.2bp.summary['1 or MH', '4 bp'] - 1
          
        }else{
          
          ### MH length = 2
          
          if(all(c(strsplit(tmp.repeat.nts, '')[[1]][3:4] == tmp.upstream.context[2:3]) ==T) | all(c(strsplit(tmp.repeat.nts, '')[[1]][1:2] == tmp.downstream.context[1:2]) == T)){
            
            dels.greater.2bp.MH.summary['2 bp MH', '4 bp'] <- dels.greater.2bp.MH.summary['2 bp MH', '4 bp'] + 1
            dels.greater.2bp.summary['1 or MH', '4 bp'] <- dels.greater.2bp.summary['1 or MH', '4 bp'] - 1
            
          }else{
            
            ### MH length = 1
            
            if(strsplit(tmp.repeat.nts, '')[[1]][4] == tmp.upstream.context[3] | strsplit(tmp.repeat.nts, '')[[1]][1] == tmp.downstream.context[1]){
              
              dels.greater.2bp.MH.summary['1 bp MH', '4 bp'] <- dels.greater.2bp.MH.summary['1 bp MH', '4 bp'] + 1
              dels.greater.2bp.summary['1 or MH', '4 bp'] <- dels.greater.2bp.summary['1 or MH', '4 bp'] - 1
              
            }
            
          }
          
        }
        
      } 
      
    }
    
  }
  
  ## 5+ bp
  if(nrow(dels.5bp) > 0){
    
    if(nrow(dels.5bp.pot.MH) > 0){
      
      repeat.nts <- str_split_fixed(str_split_fixed(dels.5bp.pot.MH[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
      upstream.context <- str_split_fixed(dels.5bp.pot.MH[,'CONTEXT'], '\\[', 2)[,1,drop=F]
      downstream.context <- str_split_fixed(str_split_fixed(dels.5bp.pot.MH[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
      for (i in 1:nrow(dels.5bp.pot.MH)){
        
        ## look at repeat
        tmp.repeat.nts <- repeat.nts[i]
        tmp.repeat.length <- nchar(tmp.repeat.nts)
        
        ## isolate first X (= repeat length - 1) nucleotides in the immediate upstream context
        tmp.upstream.context <- upstream.context[i]
        tmp.upstream.context <- strsplit(tmp.upstream.context, '')[[1]][c(length(strsplit(tmp.upstream.context, '')[[1]]) - 4):length(strsplit(tmp.upstream.context, '')[[1]])]
        
        ## isolate first X (= repeat length - 1) nucleotides in the immediate downstream context
        tmp.downstream.context <- downstream.context[i]
        tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]][1:5]
        
        ### MH length = 5+ (important to start with the highest possible MH length); only need to look at first 5 up/downstream NTs
        
        if(all(c(strsplit(tmp.repeat.nts, '')[[1]][c(tmp.repeat.length-4):tmp.repeat.length] == tmp.upstream.context[c(length(tmp.upstream.context) - 4):length(tmp.upstream.context)]) == T) | all(c(strsplit(tmp.repeat.nts, '')[[1]][1:5] == tmp.downstream.context[1:5]) == T)){
          
          dels.greater.2bp.MH.summary['5+ bp MH', '5+ bp'] <- dels.greater.2bp.MH.summary['5+ bp MH', '5+ bp'] + 1
          dels.greater.2bp.summary['1 or MH', '5+ bp'] <- dels.greater.2bp.summary['1 or MH', '5+ bp'] - 1
          
        }else{
          
          ### MH length = 4
          
          if(all(c(strsplit(tmp.repeat.nts, '')[[1]][c(tmp.repeat.length-3):tmp.repeat.length] == tmp.upstream.context[c(length(tmp.upstream.context) - 3):length(tmp.upstream.context)]) == T) | all(c(strsplit(tmp.repeat.nts, '')[[1]][1:4] == tmp.downstream.context[1:4]) == T)){
            
            dels.greater.2bp.MH.summary['4 bp MH', '5+ bp'] <- dels.greater.2bp.MH.summary['4 bp MH', '5+ bp'] + 1
            dels.greater.2bp.summary['1 or MH', '5+ bp'] <- dels.greater.2bp.summary['1 or MH', '5+ bp'] - 1
            
          }else{
            
            ### MH length = 3
            
            if(all(c(strsplit(tmp.repeat.nts, '')[[1]][c(tmp.repeat.length-2):tmp.repeat.length] == tmp.upstream.context[c(length(tmp.upstream.context) - 2):length(tmp.upstream.context)]) == T) | all(c(strsplit(tmp.repeat.nts, '')[[1]][1:3] == tmp.downstream.context[1:3]) == T)){
              
              dels.greater.2bp.MH.summary['3 bp MH', '5+ bp'] <- dels.greater.2bp.MH.summary['3 bp MH', '5+ bp'] + 1
              dels.greater.2bp.summary['1 or MH', '5+ bp'] <- dels.greater.2bp.summary['1 or MH', '5+ bp'] - 1
              
            }else{
              
              ### MH length = 2
              
              if(all(c(strsplit(tmp.repeat.nts, '')[[1]][c(tmp.repeat.length-1):tmp.repeat.length] == tmp.upstream.context[c(length(tmp.upstream.context) - 1):length(tmp.upstream.context)]) == T) | all(c(strsplit(tmp.repeat.nts, '')[[1]][1:2] == tmp.downstream.context[1:2]) == T)){
                
                dels.greater.2bp.MH.summary['2 bp MH', '5+ bp'] <- dels.greater.2bp.MH.summary['2 bp MH', '5+ bp'] + 1
                dels.greater.2bp.summary['1 or MH', '5+ bp'] <- dels.greater.2bp.summary['1 or MH', '5+ bp'] - 1
                
              }else{
                
                if(strsplit(tmp.repeat.nts, '')[[1]][tmp.repeat.length] == tmp.upstream.context[length(tmp.upstream.context)] | strsplit(tmp.repeat.nts, '')[[1]][1] == tmp.downstream.context[1]){
                  
                  dels.greater.2bp.MH.summary['1 bp MH', '5+ bp'] <- dels.greater.2bp.MH.summary['1 bp MH', '5+ bp'] + 1
                  dels.greater.2bp.summary['1 or MH', '5+ bp'] <- dels.greater.2bp.summary['1 or MH', '5+ bp'] - 1
                  
                }
                
              }
              
            }
            
          }
          
        }
        
      } 
      
    }
    
  }
  
  ## after MH cases have been taken out, rename row of deletions at 1-unit repeats in original table
  rownames(dels.greater.2bp.summary)[1] <- '1'
  
  ## 4. classify >=2 bp insertions into:
  
  # (i) 2 bp insertions at simple repeats (length 0 == "no neighbouring simple repeat")
  if(nrow(ins.2bp) > 0){
    
    ins.2bp.context <- as.character(subseq(x = reference[as.character(ins.2bp[,'CHROM'])], 
                                           start = as.numeric(ins.2bp[,'POS']) - 9, 
                                           end = as.numeric(ins.2bp[,'POS']) + 10))
    ins.2bp.context.middle <- as.character(ins.2bp[,'ALT'])
    ins.2bp.context.middle <- paste0('[', str_split_fixed(ins.2bp.context.middle, '', 2)[,2,drop=F], ']')
    ins.2bp.context.start <- str_split_fixed(ins.2bp.context, '', 11)[,1:10,drop=F]
    ins.2bp.context.start <- paste(ins.2bp.context.start[,1], ins.2bp.context.start[,2], ins.2bp.context.start[,3],
                                   ins.2bp.context.start[,4], ins.2bp.context.start[,5], ins.2bp.context.start[,6],
                                   ins.2bp.context.start[,7], ins.2bp.context.start[,8], ins.2bp.context.start[,9],
                                   ins.2bp.context.start[,10], sep = '')
    ins.2bp.context.end <- str_split_fixed(ins.2bp.context, '', 11)[,11,drop=F]
    ins.2bp[,'TRIPLET'] <- paste(ins.2bp.context.start, ins.2bp.context.middle, ins.2bp.context.end, sep = '')
    colnames(ins.2bp)[5] <- 'CONTEXT' 
    
  }
  
  # (ii) 3 bp insertions at simple repeats (length 0 == "no neighbouring simple repeat")
  if(nrow(ins.3bp) > 0){
    
    ins.3bp.context <- as.character(subseq(x = reference[as.character(ins.3bp[,'CHROM'])], 
                                           start = as.numeric(ins.3bp[,'POS']) - 14, 
                                           end = as.numeric(ins.3bp[,'POS']) + 15))
    ins.3bp.context.middle <- as.character(ins.3bp[,'ALT'])
    ins.3bp.context.middle <- paste0('[', str_split_fixed(ins.3bp.context.middle, '', 2)[,2,drop=F], ']')
    ins.3bp.context.start <- str_split_fixed(ins.3bp.context, '', 16)[,1:15,drop=F]
    ins.3bp.context.start <- paste(ins.3bp.context.start[,1], ins.3bp.context.start[,2], ins.3bp.context.start[,3],
                                   ins.3bp.context.start[,4], ins.3bp.context.start[,5], ins.3bp.context.start[,6],
                                   ins.3bp.context.start[,7], ins.3bp.context.start[,8], ins.3bp.context.start[,9],
                                   ins.3bp.context.start[,10], ins.3bp.context.start[,11], ins.3bp.context.start[,12], 
                                   ins.3bp.context.start[,13], ins.3bp.context.start[,14], ins.3bp.context.start[,15], sep = '')
    ins.3bp.context.end <- str_split_fixed(ins.3bp.context, '', 16)[,16,drop=F]
    ins.3bp[,'TRIPLET'] <- paste(ins.3bp.context.start, ins.3bp.context.middle, ins.3bp.context.end, sep = '')
    colnames(ins.3bp)[5] <- 'CONTEXT' 
    
  }
  
  # (iii) 4 bp insertions at simple repeats (length 0 == "no neighbouring simple repeat")
  if(nrow(ins.4bp) > 0){
    
    ins.4bp.context <- as.character(subseq(x = reference[as.character(ins.4bp[,'CHROM'])], 
                                           start = as.numeric(ins.4bp[,'POS']) - 19, 
                                           end = as.numeric(ins.4bp[,'POS']) + 20))
    ins.4bp.context.middle <- as.character(ins.4bp[,'ALT'])
    ins.4bp.context.middle <- paste0('[', str_split_fixed(ins.4bp.context.middle, '', 2)[,2,drop=F], ']')
    ins.4bp.context.start <- str_split_fixed(ins.4bp.context, '', 21)[,1:20,drop=F]
    ins.4bp.context.start <- paste(ins.4bp.context.start[,1], ins.4bp.context.start[,2], ins.4bp.context.start[,3],
                                   ins.4bp.context.start[,4], ins.4bp.context.start[,5], ins.4bp.context.start[,6],
                                   ins.4bp.context.start[,7], ins.4bp.context.start[,8], ins.4bp.context.start[,9],
                                   ins.4bp.context.start[,10], ins.4bp.context.start[,11], ins.4bp.context.start[,12], 
                                   ins.4bp.context.start[,13], ins.4bp.context.start[,14], ins.4bp.context.start[,15], 
                                   ins.4bp.context.start[,16], ins.4bp.context.start[,17], ins.4bp.context.start[,18], 
                                   ins.4bp.context.start[,19], ins.4bp.context.start[,20], sep = '')
    ins.4bp.context.end <- str_split_fixed(ins.4bp.context, '', 21)[,21,drop=F]
    ins.4bp[,'TRIPLET'] <- paste(ins.4bp.context.start, ins.4bp.context.middle, ins.4bp.context.end, sep = '')
    colnames(ins.4bp)[5] <- 'CONTEXT' 
    
  }
  
  # (iv) 5+ bp insertions at simple repeats (length 0 == "no neighbouring simple repeat")
  if(nrow(ins.5bp) > 0){
    
    ins.5bp.context <- as.character(subseq(x = reference[as.character(ins.5bp[,'CHROM'])], 
                                           start = as.numeric(ins.5bp[,'POS']) - 499, 
                                           end = as.numeric(ins.5bp[,'POS']) + 500))
    
    ### account for the different repeat lengths: iterate
    ins.5bp.context.middle <- as.character(ins.5bp[,'ALT'])
    ins.5bp.context.middle <- str_split_fixed(ins.5bp.context.middle, '', 2)[,2,drop=F]
    ins.5bp.context.middle <- paste0('[', ins.5bp.context.middle, ']')
    ins.5bp.context.middle.lengths <- nchar(as.character(ins.5bp[,'ALT'])) - 1
    ins.5bp.context.start <- rep(NA, nrow(ins.5bp))
    for (i in 1:length(ins.5bp.context.start)){
      tmp.ins.5bp.context.start <- str_split_fixed(ins.5bp.context[i], '', 501)[,c(1 + 500-c(5*ins.5bp.context.middle.lengths[i])):500,drop=F]
      ins.5bp.context.start[i] <- paste(tmp.ins.5bp.context.start, collapse = '')
    }
    ins.5bp.context.end <- rep(NA, nrow(ins.5bp))
    for (i in 1:length(ins.5bp.context.end)){
      tmp.ins.5bp.context.end <- str_split_fixed(ins.5bp.context[i], '', 1000)[,501:c(501 + ins.5bp.context.middle.lengths[i]*6 - 1),drop=F]
      ins.5bp.context.end[i] <- paste(tmp.ins.5bp.context.end, collapse = '')
    }
    ins.5bp[,'TRIPLET'] <- paste(ins.5bp.context.start, ins.5bp.context.middle, ins.5bp.context.end, sep = '')
    colnames(ins.5bp)[5] <- 'CONTEXT'
    
  }
  
  ## summarise >=2 bp insertions at simple repeats in matrix format
  ins.greater.2bp.summary <- matrix(0, ncol = 4, nrow = 6)
  colnames(ins.greater.2bp.summary) <- c('2 bp', '3 bp', '4 bp', '5+ bp') ## insertion size
  rownames(ins.greater.2bp.summary) <- c('0', '1', '2', '3', '4', '5+') ## number of repeats
  
  ## 2 bp
  if(nrow(ins.2bp) > 0){
    
    repeat.nts <- str_split_fixed(str_split_fixed(ins.2bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
    downstream.context <- str_split_fixed(str_split_fixed(ins.2bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
    for (i in 1:nrow(ins.2bp)){
      
      ## look at repeat
      tmp.repeat.nts <- repeat.nts[i]
      tmp.repeat.length <- nchar(tmp.repeat.nts)
      
      ## how often does it match consecutively in the immediate downstream context?
      tmp.downstream.context <- downstream.context[i]
      tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]]
      
      ## group
      tmp.downstream.context <- c(paste(tmp.downstream.context[c(tmp.repeat.length-c(tmp.repeat.length - 1)):tmp.repeat.length], collapse = ''),
                                  paste(tmp.downstream.context[c(2*tmp.repeat.length-c(tmp.repeat.length - 1)):c(2*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(3*tmp.repeat.length-c(tmp.repeat.length - 1)):c(3*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(4*tmp.repeat.length-c(tmp.repeat.length - 1)):c(4*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(5*tmp.repeat.length-c(tmp.repeat.length - 1)):c(5*tmp.repeat.length)], collapse = ''))
      
      ### check all 6 categories for upstream bases
      if(tmp.downstream.context[1] %in% tmp.repeat.nts){
        
        if(tmp.downstream.context[2] %in% tmp.repeat.nts){
          
          if(tmp.downstream.context[3] %in% tmp.repeat.nts){
            
            if(tmp.downstream.context[4] %in% tmp.repeat.nts){
              
              if(tmp.downstream.context[5] %in% tmp.repeat.nts){
                
                ins.greater.2bp.summary['5+', '2 bp'] <- ins.greater.2bp.summary['5+', '2 bp'] + 1
                
              } else{
                
                ins.greater.2bp.summary['4', '2 bp'] <- ins.greater.2bp.summary['4', '2 bp'] + 1
                
              }
              
            }else{
              
              ins.greater.2bp.summary['3', '2 bp'] <- ins.greater.2bp.summary['3', '2 bp'] + 1
              
            }
            
          }else{
            
            ins.greater.2bp.summary['2', '2 bp'] <- ins.greater.2bp.summary['2', '2 bp'] + 1
            
          }
          
        }else{
          
          ins.greater.2bp.summary['1', '2 bp'] <- ins.greater.2bp.summary['1', '2 bp'] + 1
          
        }
        
      }else{
        
        ins.greater.2bp.summary['0', '2 bp'] <- ins.greater.2bp.summary['0', '2 bp'] + 1
        
      }
      
    } 
    
  }
  
  ## 3 bp
  if(nrow(ins.3bp) > 0){
    
    repeat.nts <- str_split_fixed(str_split_fixed(ins.3bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
    downstream.context <- str_split_fixed(str_split_fixed(ins.3bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
    for (i in 1:nrow(ins.3bp)){
      
      ## look at repeat
      tmp.repeat.nts <- repeat.nts[i]
      tmp.repeat.length <- nchar(tmp.repeat.nts)
      
      ## how often does it match consecutively in the immediate downstream context?
      tmp.downstream.context <- downstream.context[i]
      tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]]
      
      ## group
      tmp.downstream.context <- c(paste(tmp.downstream.context[c(tmp.repeat.length-c(tmp.repeat.length - 1)):tmp.repeat.length], collapse = ''),
                                  paste(tmp.downstream.context[c(2*tmp.repeat.length-c(tmp.repeat.length - 1)):c(2*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(3*tmp.repeat.length-c(tmp.repeat.length - 1)):c(3*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(4*tmp.repeat.length-c(tmp.repeat.length - 1)):c(4*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(5*tmp.repeat.length-c(tmp.repeat.length - 1)):c(5*tmp.repeat.length)], collapse = ''))
      
      ### check all 6 categories for upstream bases
      if(tmp.downstream.context[1] %in% tmp.repeat.nts){
        
        if(tmp.downstream.context[2] %in% tmp.repeat.nts){
          
          if(tmp.downstream.context[3] %in% tmp.repeat.nts){
            
            if(tmp.downstream.context[4] %in% tmp.repeat.nts){
              
              if(tmp.downstream.context[5] %in% tmp.repeat.nts){
                
                ins.greater.2bp.summary['5+', '3 bp'] <- ins.greater.2bp.summary['5+', '3 bp'] + 1
                
              } else{
                
                ins.greater.2bp.summary['4', '3 bp'] <- ins.greater.2bp.summary['4', '3 bp'] + 1
                
              }
              
            }else{
              
              ins.greater.2bp.summary['3', '3 bp'] <- ins.greater.2bp.summary['3', '3 bp'] + 1
              
            }
            
          }else{
            
            ins.greater.2bp.summary['2', '3 bp'] <- ins.greater.2bp.summary['2', '3 bp'] + 1
            
          }
          
        }else{
          
          ins.greater.2bp.summary['1', '3 bp'] <- ins.greater.2bp.summary['1', '3 bp'] + 1
          
        }
        
      }else{
        
        ins.greater.2bp.summary['0', '3 bp'] <- ins.greater.2bp.summary['0', '3 bp'] + 1
        
      }
      
    } 
    
  }
  
  ## 4 bp
  if(nrow(ins.4bp) > 0){
    
    repeat.nts <- str_split_fixed(str_split_fixed(ins.4bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
    downstream.context <- str_split_fixed(str_split_fixed(ins.4bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
    for (i in 1:nrow(ins.4bp)){
      
      ## look at repeat
      tmp.repeat.nts <- repeat.nts[i]
      tmp.repeat.length <- nchar(tmp.repeat.nts)
      
      ## how often does it match consecutively in the immediate downstream context?
      tmp.downstream.context <- downstream.context[i]
      tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]]
      
      ## group
      tmp.downstream.context <- c(paste(tmp.downstream.context[c(tmp.repeat.length-c(tmp.repeat.length - 1)):tmp.repeat.length], collapse = ''),
                                  paste(tmp.downstream.context[c(2*tmp.repeat.length-c(tmp.repeat.length - 1)):c(2*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(3*tmp.repeat.length-c(tmp.repeat.length - 1)):c(3*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(4*tmp.repeat.length-c(tmp.repeat.length - 1)):c(4*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(5*tmp.repeat.length-c(tmp.repeat.length - 1)):c(5*tmp.repeat.length)], collapse = ''))
      
      ### check all 6 categories for upstream bases
      if(tmp.downstream.context[1] %in% tmp.repeat.nts){
        
        if(tmp.downstream.context[2] %in% tmp.repeat.nts){
          
          if(tmp.downstream.context[3] %in% tmp.repeat.nts){
            
            if(tmp.downstream.context[4] %in% tmp.repeat.nts){
              
              if(tmp.downstream.context[5] %in% tmp.repeat.nts){
                
                ins.greater.2bp.summary['5+', '4 bp'] <- ins.greater.2bp.summary['5+', '4 bp'] + 1
                
              } else{
                
                ins.greater.2bp.summary['4', '4 bp'] <- ins.greater.2bp.summary['4', '4 bp'] + 1
                
              }
              
            }else{
              
              ins.greater.2bp.summary['3', '4 bp'] <- ins.greater.2bp.summary['3', '4 bp'] + 1
              
            }
            
          }else{
            
            ins.greater.2bp.summary['2', '4 bp'] <- ins.greater.2bp.summary['2', '4 bp'] + 1
            
          }
          
        }else{
          
          ins.greater.2bp.summary['1', '4 bp'] <- ins.greater.2bp.summary['1', '4 bp'] + 1
          
        }
        
      }else{
        
        ins.greater.2bp.summary['0', '4 bp'] <- ins.greater.2bp.summary['0', '4 bp'] + 1
        
      }
      
    } 
    
  }
  
  ## 5+ bp
  if(nrow(ins.5bp) > 0){
    
    repeat.nts <- str_split_fixed(str_split_fixed(ins.5bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
    downstream.context <- str_split_fixed(str_split_fixed(ins.5bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
    for (i in 1:nrow(ins.5bp)){
      
      ## look at repeat
      tmp.repeat.nts <- repeat.nts[i]
      tmp.repeat.length <- nchar(tmp.repeat.nts)
      
      ## how often does it match consecutively in the immediate downstream context?
      tmp.downstream.context <- downstream.context[i]
      tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]]
      
      ## group
      tmp.downstream.context <- c(paste(tmp.downstream.context[c(tmp.repeat.length-c(tmp.repeat.length - 1)):tmp.repeat.length], collapse = ''),
                                  paste(tmp.downstream.context[c(2*tmp.repeat.length-c(tmp.repeat.length - 1)):c(2*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(3*tmp.repeat.length-c(tmp.repeat.length - 1)):c(3*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(4*tmp.repeat.length-c(tmp.repeat.length - 1)):c(4*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(5*tmp.repeat.length-c(tmp.repeat.length - 1)):c(5*tmp.repeat.length)], collapse = ''))
      
      ### check all 6 categories for upstream bases
      if(tmp.downstream.context[1] %in% tmp.repeat.nts){
        
        if(tmp.downstream.context[2] %in% tmp.repeat.nts){
          
          if(tmp.downstream.context[3] %in% tmp.repeat.nts){
            
            if(tmp.downstream.context[4] %in% tmp.repeat.nts){
              
              if(tmp.downstream.context[5] %in% tmp.repeat.nts){
                
                ins.greater.2bp.summary['5+', '5+ bp'] <- ins.greater.2bp.summary['5+', '5+ bp'] + 1
                
              } else{
                
                ins.greater.2bp.summary['4', '5+ bp'] <- ins.greater.2bp.summary['4', '5+ bp'] + 1
                
              }
              
            }else{
              
              ins.greater.2bp.summary['3', '5+ bp'] <- ins.greater.2bp.summary['3', '5+ bp'] + 1
              
            }
            
          }else{
            
            ins.greater.2bp.summary['2', '5+ bp'] <- ins.greater.2bp.summary['2', '5+ bp'] + 1
            
          }
          
        }else{
          
          ins.greater.2bp.summary['1', '5+ bp'] <- ins.greater.2bp.summary['1', '5+ bp'] + 1
          
        }
        
      }else{
        
        ins.greater.2bp.summary['0', '5+ bp'] <- ins.greater.2bp.summary['0', '5+ bp'] + 1
        
      }
      
    } 
    
  }
  
  ## 5. summarise outputs as five lists, each featuring one table per major ID category
  out <- list('1bp del' = dels.1bp.summary,
              '1bp ins' = ins.1bp.summary,
              '>2bp del' = dels.greater.2bp.summary,
              '>2bp ins' = ins.greater.2bp.summary,
              'MH dels' = dels.greater.2bp.MH.summary)
  
  ### MODIFIED (ao7) - return vector instead of list
  out <- c(out$`1bp del`[,1], out$`1bp del`[,2],
           out$`1bp ins`[,1], out$`1bp ins`[,2],
           out$`>2bp del`[,1], out$`>2bp del`[,2], out$`>2bp del`[,3], out$`>2bp del`[,4],
           out$`>2bp ins`[,1], out$`>2bp ins`[,2], out$`>2bp ins`[,3], out$`>2bp ins`[,4],
           out$`MH dels`[1,1], out$`MH dels`[1:2,2], out$`MH dels`[1:3,3], out$`MH dels`[1:5,4])
  ###
  
  return(out)
  
}

### MODIFIED (ao7)
### - deleted code below this point, including function 'plot.indel.spectrum'
