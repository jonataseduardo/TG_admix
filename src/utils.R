#' Estimate the admix thereshold 
#'
#'\code{admix_th()} estimates the alpha for every variant 
#' 
#' @import data.table
#' @export
#' 
#' @param  pop_list: list of where each element is a data.table with the
#' following columns: CHR, SNP, CM, POS, NCHROBS, POP, VAR, AF
#' 
#' @return data.table with Fst estimated for every SNP 
admix_th <-
    function(pop1_dt, pop2_dt){

      on_cols = c("CHR", "CM", "POS", "VAR", "SNP")
      scols = c(on_cols, "AF", "NCHROBS", "POP")
      pop_dt <- merge(pop1_dt[, ..scols], 
                      pop2_dt[, ..scols], 
                      by = on_cols,
                      all = TRUE)

      pop_dt[, N := .N, by = setdiff(on_cols, c('SNP'))]
      pop_dt <- pop_dt[N == 1]
      pop_dt[, N := NULL]

      pop_dt[, `:=`(HTZ.x = 2 * AF.x * (1 - AF.x), 
                    HTZ.y = 2 * AF.y * (1 - AF.y),
                    HTZ.xy = AF.x * (1 - AF.y) + AF.y * (1 - AF.x))]
      pop_dt[HTZ.x > HTZ.y, ALPHA := (HTZ.x - HTZ.y) / (HTZ.xy - HTZ.x - HTZ.y)]
      pop_dt[HTZ.y > HTZ.x, ALPHA := (HTZ.y - HTZ.x) / (HTZ.xy - HTZ.x - HTZ.y)]
      pop_dt[ALPHA < 0 | ALPHA > 1, ALPHA := NA]
      fst_dt <- hudson_fst(list(pop1_dt, pop2_dt))
      out_dt <- merge(pop_dt, fst_dt, by = setdiff(on_cols, c('VAR')), all = TRUE)
      return(out_dt[])
    }


#' Read CADD annotations files 
#'
#'\code{read_cadd_annotation()} read TG cadd annot by CHR
#' 
#' @import data.table
#' @export
#' 
#' @param fsuffix: path with the file suffix  
#' @param chr: chromossome number or name
#' @param is_compressed: boolean defaul TRUE
#' @param list of columns to be extrated
#' 
#' @return data.table with cadd annotations SNP 
read_cadd_annotation <- 
  function(fsuffix = '../data/TG_phase3/cadd_annotation/1000G_phase3_inclAnno_chr',
           chr, 
           is_compressed = TRUE, 
           select_first = TRUE, 
           annot_cols = 
             c('isTv', 'isDerived', 
                'GeneName', 'Consequence', 
                'GerpRS', 'SIFTcat', 
                'PolyPhenCat', 'RawScore','PHRED')
           ){

   fn <- 
     paste0( fsuffix,
             chr, 
             '.tsv')
   if(is_compressed)
     fn <- paste0('zcat ', fn)

    annot_cols <- 
      c('Chrom', 'Pos', 'Ref', 'Anc', 'Alt', annot_cols)

    annot_dt <- 
      fread(fn, select = annot_cols)

    if(select_first){
      uannot_dt <- unique(annot_dt[, .(Chrom, Pos)])
      annot_dt <- annot_dt[uannot_dt, 
                           on = .(Chrom, Pos), 
                           mult = 'first']
    }

    return(annot_dt)
  }



#' Read CADD annotations files 
#'
#'\code{merge_freq_and_annot()} read TG cadd annot by CHR
#' 
#' @import data.table
#' @export
#' 
#' @param fsuffix: path with the file suffix  
#' @param chr: chromossome number or name
#' 
#' @return data.table with 
merge_freq_and_annot <-
  function(freq_dt, annot_dt){

  danf <- 
    annot_dt[freq_dt[], on = c(Chrom = 'CHR', Pos = 'POS')]

  setnames(danf, c('Chrom', 'Pos'), c('CHR', 'POS'))

  return(danf[VAR == Alt])
  }


#' Calculate Derived allele frequency 
#'
#'\code{polarize_derived_allele()} read TG cadd annot by CHR
#' 
#' @import data.table
#' @export
#' 
#' @param data_in: data.table with alleles frequencies and annotations
#' @param inplace: default TRUE
#' 
#' @return data.table with 
polarize_derived_allele <-
  function(data_in, inplace = TRUE){

    if(!inplace) data_in <- copy(data_in)
    data_in[isDerived == FALSE, AF := 1 - AF]

  return(data_in[!is.na(isDerived)])
  }
