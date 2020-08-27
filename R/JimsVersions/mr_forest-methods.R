#Jim's version of mr_forest
#Initially dependent on GGPlot2
#Note package currently depends on GGPlot2 so I will continue as such

setMethod("mr_forest",
          signature = "MRInput",
          function(object,
                   alpha = 0.05,
                   snp_estimates = TRUE,
                   method = "default",
                   ordered = FALSE){
            
            
            Bx = object@betaX
            By = object@betaY
            Bxse = object@betaXse
            Byse = object@betaYse
            
            snps = object@snps
            #In case of named variants, below code
            if("snp" %in% snps) {
              snps <- paste("snp", 1:length(Bx), sep = "_")
            } else {
              snps <- snps
            }
            
            Interval_type <- paste(100*(1-alpha), "% CI)", sep = "")
            Y_label <-  paste("Variant-specific causal estimate (", Interval_type, sep = "")

            if (snp_estimates) {
              #Calculate the estimates and CIs
              estimates = By / Bx
              CI_range = qnorm(1-alpha/2)
              #Note Steve asked for / Bx but am using abs cause want CI_lower < CI_upper
              CI_lower = estimates - (CI_range * Byse)/abs(Bx)
              CI_upper = estimates + (CI_range * Byse) / abs(Bx)
              #Create the dataframe with these values
              df = data.frame(snps, estimates, CI_lower, CI_upper)
            } else {
              #Don't calculate the SNPs
              #Still need to create the data frame so we can append
              df = data.frame(matrix(ncol = 4, nrow = 0))
              names(df) = c("snps", "estimates", "CI_lower", "CI_upper")
            }
            
            

            
            
            ivw_output = mr_ivw(mr_input(Bx,Bxse,By,Byse), alpha = alpha)
            ivw_estimate = ivw_output$Estimate
            ivw_CI_lower = ivw_output$CILower
            ivw_CI_upper = ivw_output$CIUpper
            
            ivw_row = data.frame("IVW Estimate", ivw_estimate, ivw_CI_lower, ivw_CI_upper)
            names(ivw_row) = names(df)
            
            
            if (ordered == TRUE){
              factor_order = rev(df$snps[order(df$estimates)])
              df$snps = factor(df$snps, levels = factor_order)
              df = rbind(df, ivw_row)
              df$snps = factor(df$snps, levels= c("IVW Estimate", factor_order))
            } else{
              df = rbind(df, ivw_row)
              df$snps = factor(df$snps, rev(df$snps))
            }
            
            forest <- ggplot(data = df, aes(y=snps, x=estimates, xmin=CI_lower, xmax = CI_upper)) +
              geom_point(shape = 15) +
              geom_point(data = ivw_row, aes(y=snps, x = estimates), shape = 18, size = 3) +
              geom_linerange() +
              geom_vline(xintercept = 0, lty= 2) +
              ylab("Variants") + xlab(Y_label) +
              theme_classic()
            
            print(forest)
            
          })
