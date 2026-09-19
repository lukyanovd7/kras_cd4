ccc <- c('#FFFF99', '#17becf', '#E6AB02', '#1B9E77', '#f7b6d2', '#8DD3C7',
         '#E31A1C', '#BC80A0', '#FF7F00', '#33A02C', '#666666','#7570B3', 
         '#1F78B4', '#A6D854', '#FFFF33', '#A65628', '#FDB462', '#d6616b',
         '#F0027F','#a0a0a0', '#9edae5', '#b8a35a','#5b9928', '#7d6231', 
         '#c3e089', '#984ea3', '#BEAED4')
#-----------------------------------------------------------------------COLORS FOR BATCHES----------------------------------------------------------------------

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors,
                           rownames(qual_col_pals)))

col_vector[4] <- '#FDDA0D'

#-----------------------------------------------------------------------COLORS FOR CELLTYPES----------------------------------------------------------------------

subsets_cols_level1 <- c(
  'Tcm/Naive'           = '#f7b6d2',  # Cool sky blue
  'Tem Th17 like'       = '#33A02C',  # Medium green
  'Trm Th1 like'        = '#d6616b',  # Reddish-orange
  'Tem Th1 like'        = '#E66101',  # Pumpkin orange
  'Cytotoxic TEMRA'     = '#E603DC',  # Brownish-orange
  'Cytotoxic PD1 high'  = '#E31A1C',  # Bright red
  'Trm Tfh like'        = '#b8a35a',  # Deep blue
  'Tem Tfh like'        = '#33AACC',  # Cyan-teal
  'PD1 high'            = '#FDDA0D',  # Yellow
  'T reg'               = '#984EA3',  # Dark purple
  'IFN-induced'         = '#999999',  # Neutral gray
  'Cycling'             = '#E6AB02',  # Mustard yellow
  'Unknown'             = '#CCCCCC'   # Light gray
)

subsets_cols_level2 <- c(
  'Tcm/Naive'                                     = '#f7b6d2',  # Sky blue
  
  # Th17-like → Greens
  'Tem Th17 like CCR6 high IL17 low'                   = '#6BBE44',  # Light green
  'Tem Th17 like CCR6 high IL17 high'                  = '#228B22',  # Forest green
  
  # Th1-like → Oranges/Reds
  'Trm Th1 like TNF high'                         = '#8DD3C7',
  'Tem Th1 like GZMA low'                         = '#F07F44',
  'Tem Th1 like GZMA high'                        = '#d6616b',
  'Tem Th1 like GZMK high'                        = '#C94D00',
  
  # Cytotoxic
  'Cytotoxic TEMRA'                               = '#E603DC',
  'Cytotoxic PD1 high'                            = '#E31A1C',
  
  # Tfh-like → Blues/Teals
  'Trm Tfh like CXCR5 high CXCL13 low'            = '#b8a35a',
  'Tem Tfh like CXCR5 low CXCL13 high'            = '#33AACC',
  'Tem Tfh like CXCR5 high CXCL13 high'           = '#005F87',
  
  # PD1 high
  'PD1 high'                                      = '#FDDA0D',
  
  # T regs
  'T reg CD25 low'                                = '#984EA3',
  'T reg CD25 high'                               = '#BEAED4',
  
  # IFN-induced
  'IFN-induced'                                   = '#999999',
  
  # Cycling
  'Cycling'                                       = '#E6AB02',
  
  # Unknown
  'Unknown'                                       = '#CCCCCC'
)


