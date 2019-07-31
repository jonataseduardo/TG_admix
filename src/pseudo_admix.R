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
l_chr <- c(10:13)

{
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
                D_HTZ_m = roll_sumr(D_HTZ, 150, na.rm = TRUE) / 150.,
                HTZ_asw = roll_sumr(HTZ, 150, na.rm = TRUE) / 150.,
                HTZ_yri = roll_sumr(HTZ.x, 150, na.rm = TRUE) / 150.,
                HTZ_ceu = roll_sumr(HTZ.y, 150, na.rm = TRUE) / 150.),
         by = CHR]
alpha_dt[, FST_m := T1_m / T2_m]

alpha_dt[, .(sum(HTZ, na.rm = TRUE )/.N, 
             sum(HTZ.x, na.rm = TRUE )/.N, 
             sum(HTZ.y, na.rm = TRUE )/.N,
             mean(HTZ, na.rm = TRUE ), 
             mean(HTZ.x, na.rm = TRUE ), 
             mean(HTZ.y, na.rm = TRUE ))]


S_wd <- 
  alpha_dt[, .(S_asw= length(.I[which(!is.na(AF))]), 
               S_yri= length(.I[which(!is.na(AF.x))]), 
               S_ceu= length(.I[which(!is.na(AF.y))]), 
               HTZ_asw = mean(HTZ, na.rm = TRUE), 
               HTZ_yri = mean(HTZ.x, na.rm = TRUE), 
               HTZ_ceu = mean(HTZ.y, na.rm = TRUE), 
               D_HTZ = mean(D_HTZ, na.rm = TRUE), 
               HTZn_asw = sum(HTZ, na.rm = TRUE) / 2e4, 
               HTZn_yri = sum(HTZ.x, na.rm = TRUE) / 2e4, 
               HTZn_ceu = sum(HTZ.y, na.rm = TRUE) / 2e4, 
               D_HTZn =sum(D_HTZ, na.rm = TRUE) / 2e4, 
               FST = sum(T1, na.rm = TRUE)/sum(T2, na.rm = TRUE) 
              )
           ,by = floor(POS / 2e4)][1:(.N - 1)]

S_wd

S_wd[,cor(S_asw / S_yri , FST)]
S_wd[, summary(S_asw)]
S_wd[ S_asw > S_yri, .N] / S_wd[, .N]
S_wd[ S_asw > S_ceu, .N] / S_wd[, .N]

S_wd_bulk <- 
  alpha_dt[!is.na(AF.x) & !is.na(AF.y), 
           .(S_asw= length(.I[which(!is.na(AF))]), 
             S_yri= length(.I[which(!is.na(AF.x))]), 
             S_ceu= length(.I[which(!is.na(AF.y))]), 
             D_HTZ = mean(D_HTZ, na.rm = TRUE), 
             HTZ_asw = mean(HTZ, na.rm = TRUE), 
             HTZ_yri = mean(HTZ.x, na.rm = TRUE), 
             HTZ_ceu = mean(HTZ.y, na.rm = TRUE), 
             FST = sum(T1, na.rm = TRUE)/sum(T2, na.rm = TRUE) 
              )
           ,by = floor(POS / 2e4)][1:(.N - 1)]


alpha_dt[D_HTZ > 0, .N] / alpha_dt[, .N]

alpha_dt[, .(sum(HTZ, na.rm = TRUE) / .N, 
             sum(HTZ.x, na.rm = TRUE) / .N, 
             sum(HTZ.y, na.rm = TRUE) / .N, 
             sum(D_HTZ_m, na.rm = TRUE) / .N)]

alpha_dt[, .(mean(HTZ, na.rm = TRUE) , 
             mean(HTZ.x, na.rm = TRUE) , 
             mean(HTZ.y, na.rm = TRUE) , 
             mean(D_HTZ_m, na.rm = TRUE) )]

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

p5 <- 
  ggplot(alpha_dt[!is.na(D_HTZ_m) & !is.na(FST_m)], aes(x = FST_m, y = HTZ_asw - HTZ_yri)) + 
    theme_bw() + 
    theme(legend.position = "bottom") +
    geom_hex(aes(fill = log10(..count..))) + 
    geom_smooth(method = 'lm', se = FALSE, color = 'red') + 
    geom_hline(yintercept = 0, color = 'black') 
ggsave('~/teste5.png', plot = p5)

p6 <- 
  ggplot(alpha_dt[!is.na(D_HTZ_m) & !is.na(FST_m)], aes(x = FST_m, y = HTZ_asw - HTZ_ceu)) + 
    theme_bw() + 
    theme(legend.position = "bottom") +
    geom_hex(aes(fill = log10(..count..))) + 
    geom_smooth(method = 'lm', se = FALSE, color = 'red') + 
    geom_hline(yintercept = 0, color = 'black') 
ggsave('~/teste6.png', plot = p6)

p7 <- 
  ggplot(S_wd, aes(x = HTZ_asw, y = S_asw / S_yri)) + 
    theme_bw() + 
    theme(legend.position = "bottom") +
    geom_hex(aes(fill = log10(..count..))) +
    geom_smooth(method = 'lm', se = FALSE, color = 'red') + 
    geom_hline(yintercept = 1, color = 'black') 
ggsave('~/teste7.png', plot = p7)

p8 <- 
  ggplot(S_wd, aes(x = HTZ_yri, y = S_asw / S_yri)) + 
    theme_bw() + 
    theme(legend.position = "bottom") +
    geom_hex(aes(fill = log10(..count..))) +
    geom_smooth(method = 'lm', se = FALSE, color = 'red') + 
    geom_hline(yintercept = 1, color = 'black') 
ggsave('~/teste8.png', plot = p8)

p9 <- 
  ggplot(S_wd, aes(x = HTZ_ceu, y = S_asw / S_ceu)) + 
    theme_bw() + 
    theme(legend.position = "bottom") +
    geom_hex(aes(fill = log10(..count..))) +
    geom_smooth(method = 'lm', se = FALSE, color = 'red') + 
    geom_hline(yintercept = 1, color = 'black') 
ggsave('~/teste9.png', plot = p9)

p10 <- 
  ggplot(alpha_dt[!is.na(D_HTZ_m) & !is.na(FST_m)], aes(x = FST_m, y = HTZ_yri - HTZ_ceu)) + 
    theme_bw() + 
    theme(legend.position = "bottom") +
    geom_hex(aes(fill = log10(..count..))) + 
    geom_smooth(method = 'lm', se = FALSE, color = 'red') + 
    geom_hline(yintercept = 0, color = 'black') 
ggsave('~/teste10.png', plot = p10)

p11 <- 
  ggplot(S_wd, aes(x = HTZ_asw - HTZ_yri, y = S_asw / S_yri)) + 
    theme_bw() + 
    theme(legend.position = "bottom") +
    geom_hex(aes(fill = log10(..count..))) +
    geom_smooth(method = 'lm', se = FALSE, color = 'red') + 
    geom_hline(yintercept = 1, color = 'black') 
ggsave('~/teste11.png', plot = p11)


p12 <- 
  ggplot(S_wd, aes(x = HTZ_asw - HTZ_ceu, y = S_asw / S_ceu)) + 
    theme_bw() + 
    theme(legend.position = "bottom") +
    geom_hex(aes(fill = log10(..count..))) +
    geom_smooth(method = 'lm', se = FALSE, color = 'red') + 
    geom_hline(yintercept = 1, color = 'black') 
ggsave('~/teste12.png', plot = p12)

p13 <- 
  ggplot(S_wd, aes(y = HTZ_asw - HTZ_yri, x = FST)) + 
    theme_bw() + 
    theme(legend.position = "bottom") +
    geom_hex(aes(fill = log10(..count..))) +
    geom_smooth(method = 'lm', se = FALSE, color = 'red') + 
    geom_hline(yintercept = 0, color = 'black') 
ggsave('~/teste13.png', plot = p13)


p14 <- 
  ggplot(S_wd, aes(y = HTZ_asw - HTZ_ceu, x = FST)) + 
    theme_bw() + 
    theme(legend.position = "bottom") +
    geom_hex(aes(fill = log10(..count..))) +
    geom_smooth(method = 'lm', se = FALSE, color = 'red') + 
    geom_hline(yintercept = 0, color = 'black') 
ggsave('~/teste14.png', plot = p14)

p15 <- 
  ggplot(S_wd, aes(y = HTZn_asw - HTZn_ceu, x = FST)) + 
    theme_bw() + 
    theme(legend.position = "bottom") +
    geom_hex(aes(fill = log10(..count..))) +
    geom_smooth(method = 'lm', se = FALSE, color = 'red') + 
    geom_hline(yintercept = 0, color = 'black') 
ggsave('~/teste15.png', plot = p15)

p16 <- 
  ggplot(S_wd, aes(y = HTZn_asw - HTZn_yri, x = FST)) + 
    theme_bw() + 
    theme(legend.position = "bottom") +
    geom_hex(aes(fill = log10(..count..))) +
    geom_smooth(method = 'lm', se = FALSE, color = 'red') + 
    geom_hline(yintercept = 0, color = 'black') 
ggsave('~/teste16.png', plot = p16)

p17 <- 
  ggplot(S_wd, aes(y = D_HTZn, x = FST)) + 
    theme_bw() + 
    theme(legend.position = "bottom") +
    geom_hex(aes(fill = log10(..count..))) +
    geom_smooth(method = 'lm', se = FALSE, color = 'red') + 
    geom_hline(yintercept = 0, color = 'black') 
ggsave('~/teste17.png', plot = p17)

p18 <- 
  ggplot(S_wd, aes(y = D_HTZ, x = FST)) + 
    theme_bw() + 
    theme(legend.position = "bottom") +
    geom_hex(aes(fill = log10(..count..))) +
    geom_smooth(method = 'lm', se = FALSE, color = 'red') + 
    geom_hline(yintercept = 0, color = 'black') 
ggsave('~/teste18.png', plot = p18)

p19 <- 
  ggplot(S_wd, aes(x = HTZn_asw - HTZn_yri, y = S_asw / S_yri)) + 
    theme_bw() + 
    theme(legend.position = "bottom") +
    geom_hex(aes(fill = log10(..count..))) +
    geom_smooth(method = 'lm', se = FALSE, color = 'red') + 
    geom_hline(yintercept = 1, color = 'black') 
ggsave('~/teste19.png', plot = p19)


p20 <- 
  ggplot(S_wd, aes(x = HTZn_asw - HTZn_ceu, y = S_asw / S_ceu)) + 
    theme_bw() + 
    theme(legend.position = "bottom") +
    geom_hex(aes(fill = log10(..count..))) +
    geom_smooth(method = 'lm', se = FALSE, color = 'red') + 
    geom_hline(yintercept = 1, color = 'black') 
ggsave('~/teste20.png', plot = p20)
}

S_wd
