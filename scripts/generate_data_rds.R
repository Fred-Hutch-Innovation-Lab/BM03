library(tidyverse)  ## General R logic
library(here)       ## Easier specification of file locations
library(immunarch)  ## MixCR output parsing and stats

main <- function() {
  parse_name <- function(x) {
    x <- as.data.frame(str_match(x, '(^.+-F\\d[^_]*)_+([^_]*)?'))
    colnames(x) <- c(
      'File', 'Sample', 'Chain'
    )
    x
  }
  
  load_data <- function(data_dir, metadata=metadata, ...) {
    # rep_data <- repLoad(data_dir)
    rep_data <- repLoad(data_dir, ...)
    
    rep_data$meta <- parse_name(rep_data$meta$Sample)
    rep_data
  }
  
  metadata <- read.csv(here('config/metadata.csv'), header = TRUE)
  rep_data1 <- load_data(here('../processing/mixcr/RNA/downsampled_data/by_reads/'), .count='uniqueMoleculeCount')
  rep_data1$meta <- rbindlist(lapply(names(rep_data1$data), parse_name))
  
  rep_data2 <- load_data('/fh/fast/_IRC/FHIL/grp/BM03/processing/adaptive')
  
  rep_data <- list()
  rep_data$data <- c(rep_data1$data, rep_data2$data)
  rep_data$meta <- rbind(rep_data1$meta, rep_data2$meta)
  
  rep_data$meta <- merge(
    rep_data$meta, 
    metadata,
    by='Sample'
  ) %>%
    mutate(Chain = ifelse(is.na(Chain), 'other', Chain)) %>%
    mutate(Chain = factor(Chain, levels=c(
      'TRA', 'TRB',
      'TRD', 'TRG',
      'ILH', 'IGK', 'IGH', 'IGL',
      'other'
    )
   )
  )
  
  saveRDS(rep_data, here('rds/rep_data_downsampled_umi.rds'), compress=FALSE)  
  return(rep_data)
}

rep_data <- main()
rm(main)