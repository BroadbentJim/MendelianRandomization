#Jim's version of mr_mvmedian-methods.R

#' @docType methods
#' @rdname mr_mvmedian

setMethod(mr_mvmedian,
          "MRMVInput",
          function(object,
                   boot = FALSE,
                   boot_it = 1000,
                   model = "default",
                   correl = FALSE,
                   distribution = "normal",
                   alpha = 0.05, ...){
            
            Bx = object@betaX
            By = object@betaY
            Bxse = object@betaXse
            Byse = object@betaYse
            rho = object@correlation
            
            nsnps <- dim(Bx)[1]
            
            qr_mod = rq(By ~ Bx - 1, weights = Byse^-2)
            if (boot == TRUE){
              #Boot
              est = sapply(1:boot_it, function(i){
                p = length(By)
                k = dim(Bx)[2]
                Sx = lapply(1:p, function(j){diag(Bxse[j, ]^2)})
                bxboot = sapply(1:p, function(j){mvrnorm(1, Bx[j, ], Sx[[j]])})
                bxboot = t(bxboot)
                byboot = rnorm(p, By, Byse)
                rq(byboot ~ bxboot - 1, weights = Byse^-2)$coefficients
              })
              se = apply(est, 1, sd)
              thetaMed = qr_mod$coefficients
              thetaMedse = se
              
              if(distribution == "normal"){
                ciLower <- ci_normal("l", thetaMed, thetaMedse, alpha)
                ciUpper <- ci_normal("u", thetaMed, thetaMedse, alpha)
              } else if (distribution == "t-dist"){
                ciLower <- ci_t("l", thetaMed, thetaMedse, nsnps - dim(Bx)[2], alpha)
                ciUpper <- ci_t("u", thetaMed, thetaMedse, nsnps - dim(Bx)[2], alpha)
              }
              
              if (distribution == "normal") { pvalue <- 2*pnorm(-abs(thetaMed/thetaMedse)) }
              if (distribution == "t-dist") { pvalue <- 2*pt(-abs(thetaMed/thetaMedse), df=nsnps-dim(Bx)[2]) }
              return(new("MVMedian",
                         Model = model,
                         Exposure = object@exposure,
                         Outcome = object@outcome,
                         
                         Correlation = object@correlation,
                         
                         Estimate = as.numeric(thetaMed),
                         StdError = as.numeric(thetaMedse),
                         CILower =  as.numeric(ciLower),
                         CIUpper = as.numeric(ciUpper),
                         
                         SNPs = nsnps,
                         Pvalue = as.numeric(pvalue),
                         
                         Alpha = alpha,
                         
                         RSE = 0,
                         Heter.Stat = c(0,
                                        0)))

              }
              
              else {
              #Not boot
              #SE not calculated
              
              return(list("coefficients" = qr_mod$coefficients))
            }
            
            
            return(new("MVMedian",
                       Model = model,
                       ))
            
            
            
            
          }
)
