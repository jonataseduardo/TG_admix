library(data.table)
library(diverdt)
library(ggplot2)
library(ggridges)
library(devtools)
library(pkgload)
library(hexbin)
library(ggExtra)

devtools::install_github('jonataseduardo/diverdt')
reload(pkgload::inst('diverdt'))
help(package = 'diverdt')

source('utils.R')

#pops <- c('AFR', 'EUR', 'ACB', 'ASW')
pops <- c('YRI', 'CEU', 'ASW')
l_chr <- c(9)

pop_freqs <- 
    lapply(
      pops, 
      function(pop_name){
        rbindlist(
          lapply( l_chr, 
            function(chr){
              fname <- paste0('../data/TG_phase3/',
                         pop_name, '/', pop_name, '_chr', chr)
              load_bim_frq(fname, pop_name = pop_name )
            })
        )
      }
    )

system.time(
annot_dt <- 
  rbindlist(lapply(l_chr, function(chr) read_cadd_annotation(chr = chr)))
)

annot_dt
pop_freqs[[1]]

yri_dt <- 
  polarize_derived_allele(
    merge_freq_and_annot(pop_freqs[[1]], annot_dt)
    )

ceu_dt <- 
  polarize_derived_allele(
    merge_freq_and_annot(pop_freqs[[2]], annot_dt)
    )

asw_dt <- 
  polarize_derived_allele(
    merge_freq_and_annot(pop_freqs[[3]], annot_dt)
    )
 
source('utils.R')
pop1_dt <- yri_dt[1:100]
pop2_dt <- ceu_dt[1:100]

on_cols = c("CHR", "CM", "POS", "VAR", "SNP")
scols = c(on_cols, "AF", "NCHROBS", "POP")
system.time(
alpha_dt <- merge(asw_dt[, ..scols], 
                  admix_th(yri_dt, ceu_dt), 
                  by = on_cols, 
                  all = TRUE)
)
alpha_dt[, HTZ := 2 * AF * (1 - AF)]
alpha_dt[, D_HTZ := HTZ - max(HTZ.x, HTZ.y), by = on_cols]

alpha_dt[ D_HTZ > 0.1 & FST < 0.25, .(mean(AF), mean(AF.x), mean(AF.y))]

p <- 
  ggplot(alpha_dt[!is.na(D_HTZ)], aes(x = FST, y = D_HTZ)) + 
    theme_bw() + 
    theme(legend.position = "bottom") +
    geom_point(alpha = 0)  + 
    geom_hex(aes(fill = log10(..count..)))
p2 <- ggMarginal(p, type = 'histogram')
  ggsave('teste.png', plot = p2)
