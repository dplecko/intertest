
ricu:::init_proj()

# test if it runs through
# one_step_debias(mtcars[rep(1:nrow(mtcars), 10), ], 
#                 X = "vs", Z = NULL, W = "disp", Y = "mpg")
# 
# 
# one_step_debias(mtcars[rep(1:nrow(mtcars), 10), ], 
#                 X = "vs", Z = NULL, W = NULL, Y = "am", log_risk = TRUE)

n <- 10^3
log_risk <- TRUE
ret <- NULL
for (sclass in LETTERS[c(6, 7, 10, 11)]) {
  
  for (rep in seq_len(10)) {
    
    if (is.element(sclass, c("F", "G", "J", "K"))) {
      
      c(data, SFM, gt, prob_y) %<-% gen_from_scm(sclass, n = n)
    } else c(data, SFM, gt) %<-% gen_from_scm(sclass, n = n)
    
    c(X, Z, W, Y) %<-% SFM
    gt <- ia_gt(sclass, log_risk = log_risk)
    for (opt in c(1, 2)) {
      
      if (opt == 1 & is.element(sclass, c("J", "K"))) next
      osd_est <- if (opt == 1) {
        
        one_step_debias(data, X, Z, W, Y, log_risk = log_risk)
      } else one_step_debias2(data, X, Z, W, Y, log_risk = log_risk)
      
      iter <- merge(osd_est, gt, by = c("measure", "scale"))
      iter$rep <- rep
      iter$option <- opt
      ret <- rbind(ret, iter)
      
      cat(sclass, "-", rep, "-", opt, "\n", sep = "")
    }
  }
}

ret <- as.data.table(ret)
ret[, pval := 2 * pnorm(-abs((gt_value - value) / sd))]
ret[gt_value == 0 & value == 0, pval := runif(.N)]

ggplot(ret, aes(x = pval, color = factor(option))) +
  stat_ecdf() + theme_bw() +
  facet_wrap(~measure) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "orange")
