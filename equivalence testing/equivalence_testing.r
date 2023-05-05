library(TOSTER)
rm(list = ls())


df <- read.csv(
  "C:\\My Drive\\Projects\\cgrc md\\codebase\\data\\03_cgrc_data\\sbmd_cgrA__cgrc_data.csv")

scales = c('PANAS', 'QIDS', 'mood', 'creativity')
eqbs <- list(PANAS=8.2, mood=17.8, creativity=15.4, QIDS=2.9)

for (scale in scales){
  
  eval(parse(text = paste('eq_bound <- eqbs$', scale, sep="")))
  is_sig <- c()
  
  for (cgr_trial_id in unique(df$cgr_trial_id)) {
    
    # do not understand why, but in subset()  scale==scale does not work
    # a, b as intermediary variables are a temporary workaround
    a <- cgr_trial_id 
    b <- scale
    
    pls <- subset(df, condition=='PL' & cgr_trial_id==a & scale==b)
    acs <- subset(df, condition=='AC' & cgr_trial_id==a & scale==b)
    
    res <- tsum_TOST(
      hypothesis = 'EQU',
      paired = FALSE,
      var.equal = FALSE,
      eqb = eq_bound,
      eqbound_type = 'raw',
      m1 = mean(pls$delta_score),
      m2 = mean(acs$delta_score),
      sd1 = sd(pls$delta_score),
      sd2 = sd(acs$delta_score),
      n1 = nrow(pls),
      n2 = nrow(acs),)
    
    print(res$decision$TOST)
    is_sig <- c(is_sig, res$TOST$p.value[3] < 0.05)
  }
  print(paste(scale,':',round(100*sum(is_sig)/length(is_sig),2),'% of CGR-adjs indicate equivalence with eq. bound of', eq_bound))
}