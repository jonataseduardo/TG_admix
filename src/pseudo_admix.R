library(data.table)
library(diverdt)
library(ggplot2)
library(ggridges)
library(devtools)
library(pkgload)
library(hexbin)
library(ggExtra)
library(RcppRoll)

devtools::install_github('jonataseduardo/diverdt')
reload(pkgload::inst('diverdt'))
help(package = 'diverdt')

source('utils.R')

#pops <- c('AFR', 'EUR', 'ACB', 'ASW')
pops <- c('YRI', 'CEU', 'ASW')
l_chr <- c(18)

pop_freqs <- 
    lapply(
      pops, 
      function(pop_name){
        rbindlist(
          lapply( l_chr, 
            function(chr){
              fname <- paste0('../data/TG_phase3/',
                         pop_name, '61/', pop_name, '_chr', chr)
              load_bim_frq(fname, pop_name = pop_name )
            })
        )
      }
    )

system.time(
annot_dt <- 
  rbindlist(lapply(l_chr, function(chr) read_cadd_annotation(chr = chr)))
)

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
alpha_dt[, `:=`(T1_m = frollmean(T1, 150, na.rm = TRUE),
                T2_m = frollmean(T2, 150, na.rm = TRUE),
                D_HTZ_m = roll_sumr(D_HTZ, 150, na.rm = TRUE) / 150.),
         by = CHR]

S_wd <- 
  alpha_dt[, .(S_asw= length(.I[which(!is.na(AF))]), 
               S_yri= length(.I[which(!is.na(AF.x))]), 
               S_ceu= length(.I[which(!is.na(AF.y))]), 
               D_HTZ = mean(D_HTZ, na.rm = TRUE), 
               FST = sum(T1, na.rm = TRUE)/sum(T2, na.rm = TRUE) 
              )
           ,by = floor(POS / 2e4)][1:(.N - 1)]


S_wd[,cor(S_asw / S_yri , FST)]
S_wd[, summary(S_asw)]

S_wd_bulk <- 
  alpha_dt[!is.na(AF.x) & !is.na(AF.y), 
           .(S_asw= length(.I[which(!is.na(AF))]), 
             S_yri= length(.I[which(!is.na(AF.x))]), 
             S_ceu= length(.I[which(!is.na(AF.y))]), 
             D_HTZ = mean(D_HTZ, na.rm = TRUE), 
             FST = sum(T1, na.rm = TRUE)/sum(T2, na.rm = TRUE) 
              )
           ,by = floor(POS / 2e4)][1:(.N - 1)]


alpha_dt[D_HTZ_m > 0, .N] / alpha_dt[, .N]

alpha_dt[, .(sum(HTZ, na.rm = TRUE) / .N, sum(HTZ.x, na.rm = TRUE) /.N , sum(HTZ.y, na.rm = TRUE) /.N)]
alpha_dt[, .(sum(HTZ, na.rm = TRUE) / (max(POS) -  min(POS)) , 
             sum(HTZ.x, na.rm = TRUE) /(max(POS) -  min(POS)) , 
             sum(HTZ.y, na.rm = TRUE) /(max(POS) -  min(POS)))]
alpha_dt[is.na(HTZ), .N] / alpha_dt[, .N]

S_wd_bulk

{
p0 <- 
  ggplot(alpha_dt[!is.na(D_HTZ) & !is.na(FST)], aes(x = FST, y = D_HTZ)) + 
    theme_bw() + 
    theme(legend.position = "bottom") +
    geom_hex(aes(fill = log10(..count..))) + 
    geom_smooth(method = 'lm', se = FALSE, color = 'red') + 
    geom_hline(yintercept = 0, color = 'black') 
ggsave('~/teste0.png', plot = p0)

p1 <- 
  ggplot(alpha_dt[!is.na(D_HTZ_m) & !is.na(FST_m)], aes(x = FST_m, y = D_HTZ_m)) + 
    theme_bw() + 
    theme(legend.position = "bottom") +
    geom_hex(aes(fill = log10(..count..))) + 
    geom_smooth(method = 'lm', se = FALSE, color = 'red') + 
    geom_hline(yintercept = 0, color = 'black') 
ggsave('~/teste1.png', plot = p1)

p2 <- 
  ggplot(S_wd, aes(x = S_asw / S_yri, y = S_asw / S_ceu)) + 
    theme_bw() + 
    theme(legend.position = "bottom") +
    geom_hex(aes(fill = log10(..count..))) +
    geom_smooth(method = 'lm', se = FALSE, color = 'red') + 
    geom_hline(yintercept = 1, color = 'black') +
    geom_vline(xintercept = 1, color = 'black') 
ggsave('~/teste2.png', plot = p2)

p3 <- 
  ggplot(S_wd, aes(x = D_HTZ, y = S_asw / S_yri)) + 
    theme_bw() + 
    theme(legend.position = "bottom") +
    geom_hex(aes(fill = log10(..count..))) +
    geom_smooth(method = 'lm', se = FALSE, color = 'red') + 
    geom_hline(yintercept = 1, color = 'black') 
ggsave('~/teste3.png', plot = p3)

p4 <- 
  ggplot(S_wd, aes(x = FST, y = S_asw / S_yri)) + 
    theme_bw() + 
    theme(legend.position = "bottom") +
    geom_hex(aes(fill = log10(..count..))) +
    geom_smooth(method = 'lm', se = FALSE, color = 'red') + 
    geom_hline(yintercept = 1, color = 'black') 
ggsave('~/teste4.png', plot = p4)
}
