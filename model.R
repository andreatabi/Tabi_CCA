rm(list=ls())

library(foreach)
library(doParallel)

numCores <- detectCores()-1
cl <- makeCluster(numCores, outfile="")
registerDoParallel(cl)  
clusterCall(cl, function() {
  library(dplyr)
  library(entropy)
  library(scales) 
  library(varrank)
})

source("~/toolbox.R")     

## STOCHASTIC SIMULATIONS ---------------------------------------------------------------------

set.seed(19850827)
time <-  31*1e04          
S <- 55         
mu <- sort(rep(c(0, 0.1, 0.5, 1, 1.5), 100))  
frac <- 0.4  
sims <- length(mu)

results <- foreach(u = 1:sims, .combine = rbind) %dopar% {
      
      print(u)
      size <- 10^rnorm(S,1,0.3)  
      size <- sort(size)
      int <- interaction_matrix(size, frac=frac, mu=mu[u])
      A <- int$A
      pred_all <- int$pred
      fr <- length(pred_all)/S
      q <- size^(-0.25 + rnorm(S, 0, 0.12) )        
      K <- 1000 * size^(-0.75 + rnorm(S, 0, 0.1))    
      d <- rep(0.5,S)                                 
      x <- rep(10,S)                                
      out <- array(0, dim = c(time, S+1) )
      
      for(i in 1:time){
        out[i,-1] <- x
        r <- size^(-0.25 + rnorm(S, 0, 0.12))      
        B <- q + r * x * ( 1 - x / K )
        D <- d * x + x * ( c(A %*% x) / K )
        totalrate <- B + D
        sumtotalrate <- sum(totalrate)
        if(i == max(time)) break
        z1 <- runif(1)
        tau <- 1/sumtotalrate * log(1/z1)
        out[i+1,1] <- out[i,1] + tau
        z2 <- runif(1)
        for(j in 1:S){
          ss <- sum(totalrate[1:j-1]) <= z2*sumtotalrate && z2*sumtotalrate < sum(totalrate[1:j])
          if( ss==TRUE ){
            ifelse(sum(totalrate[1:j-1]) + B[j] > z2*sumtotalrate, 
                   x[j] <-  x[j] + 1, 
                   x[j] <-  x[j] - 1)
          }
        }
      
      }

      n_sample <- 30                            
      inter <- max(out[,1])/(n_sample)              
      res <- matrix(NA, nrow=c(n_sample), ncol=c(S+1))
      res[,1] <- seq(inter, max(out[,1]), by=inter)
      for(i in 1:n_sample){
        ss <- min(which(res[i,1] <= out[,1]))
        res[i,-1] <- out[ss,-1]
      }

      rm(out)
      ti <- seq(n_sample)
      gl <- GL_func(res[,-1],ti)
      co <- tryCatch(cor.test(gl$gains, gl$losses, method = "spearman")$est , error = function(e) {NA})
      mi <- tryCatch(mut.info(gl$gains, gl$losses, discretization.method = "sturges"), error = function(e) {NA})
      TO <- tryCatch(mean(gl$turnover, na.rm=T), error = function(e) {NA})
      SR <- apply(res[,-1], 1, function(x)  length(which(x>0)))[-1]
      even <- apply(res[,-1], 1, function(x) evenness(x )$J )[-1]
      corR <- cor.test(size, 1/r, method="spearman")$est 
      corDisp <- cor.test(size, q, method="spearman")$est 
      corCC <- cor.test(size, K, method="spearman")$est 
      preds <- apply(res[,-1], 1, function(x) length(which( which(x>0) %in% pred_all)) )[-1]
      pred_prop <- preds / SR
      
      diag(A) <- NA
      
      df <- data.frame(sim = u,
                       SP_pool=S, 
                       mu=mu[u], 
                       pred_frac=fr,
                       sdA = sd(A, na.rm=T),
                       meanA = mean(A, na.rm=T),
                       cor = co, 
                       mi = mi,
                       turnover = TO, 
                       rich = mean(SR, na.rm=T),
                       even= mean(even, na.rm=T),
                       pred_prop = mean(pred_prop, na.rm=T),
                       corVolt.size = mean(corR, na.rm=T),  
                       corDisp.size = mean(corDisp, na.rm=T),
                       corCC.size = mean(corCC, na.rm=T),
                       meanM = mean(size)
      )
      return(df)
}

stopCluster(cl)

results

save(results, file="~/sim.RData") 

## FIGURE 4 -----------------------------------------------------------------------------------
load("~/sim.RData") 
load("~/emp.RData")  
library(ggplot2)

df1 <- lapply(seq(length(db)), function(x) db[[x]][[1]] )
res_df <- do.call("rbind", df1)

df2 <- lapply(seq(length(db)), function(x) db[[x]][[2]] )
res_df2 <- do.call("rbind", df2)

dd <- results
head(dd)
panelA <- ggplot(dd, aes(x=factor(mu), y=cor, fill=factor(mu))) + 
  geom_violin( color="white"  ) +
  geom_boxplot(outlier.shape = NA, width=0.2, color="grey80" )+
  theme_classic()+
  ylim(-1,1)+
  scale_fill_manual(values=scales::alpha(c("#E40303", "salmon", "#FF8C00", "darkmagenta", "#24408E"),0.6))+
  geom_hline(aes(yintercept=median(res_df$cor, na.rm = T)), linewidth=0.8, color=scales::alpha("grey20",1), linetype="dashed")+
  ylab(expression(paste( "Spearman's coefficient ",  sep="")))+
  xlab(expression(paste( "Interaction strength ",  sep="")))+
  theme(legend.position="none", axis.title = element_text(size=17), 
        axis.text = element_text(size=12)  )
panelA

panelB <- ggplot(dd, aes(x=factor(mu), y=mi, fill=factor(mu) )) + 
  geom_violin( color="white") +
  geom_boxplot( outlier.shape = NA, width=0.2, color="grey80" )+
  theme_classic()+
  scale_fill_manual(values=scales::alpha(c("#E40303", "salmon", "#FF8C00", "darkmagenta", "#24408E"),0.6))+
  ylim(0,1)+
  geom_hline(aes(yintercept=median(res_df$mi, na.rm = T)), linewidth=0.8, color=scales::alpha("grey20",1), linetype="dashed")+
  ylab(expression(paste( "Mutual information",  sep="")))+
  xlab(expression(paste( "Interaction strength ",  sep="")))+
  theme(legend.position="none", axis.title = element_text(size=17), 
        axis.text = element_text(size=12)  )
panelB

panelC <- ggplot(dd, aes(x=factor(mu), y=rich/SP_pool, fill=factor(mu))) + 
  geom_violin( color="white") +
  geom_boxplot( outlier.shape = NA, width=0.2, color="grey80" )+
  theme_classic()+
  ylim(0,1)+
  scale_fill_manual(values=scales::alpha(c("#E40303", "salmon", "#FF8C00", "darkmagenta", "#24408E"),0.6))+
  geom_hline(aes(yintercept=median(res_df$sr/res_df$sp.pool, na.rm = T)), linewidth=0.8, color=scales::alpha("grey20",1), linetype="dashed")+
  ylab(expression(paste( "Relative richness",  sep="")))+
  xlab(expression(paste( "Interaction strength " ,  sep="")))+
  theme(legend.position="none", axis.title = element_text(size=17), 
        axis.text = element_text(size=12)  )
panelC

panelD <- ggplot(dd, aes(x=factor(mu), y=even, fill=factor(mu))) + 
  geom_violin(color="white") +
  geom_boxplot( outlier.shape = NA, width=0.2, color="grey70" )+
  theme_classic()+
  scale_fill_manual(values=scales::alpha(c("#E40303", "salmon", "#FF8C00", "darkmagenta", "#24408E"),0.6))+
  ylim(0,1)+
  geom_hline(aes(yintercept=median(res_df$even, na.rm = T)), linewidth=0.8, color=scales::alpha("grey20",1), linetype="dashed")+
  ylab(expression(paste( "Species evenness",  sep="")))+
  xlab(expression(paste( "Interaction strength ",  sep="")))+
  theme(legend.position="none", axis.title = element_text(size=17), 
        axis.text = element_text(size=12)  )
panelD

panelE <- ggplot(dd, aes(x=factor(mu),y=turnover, fill=factor(mu))) + 
  geom_violin(color="white") +
  geom_boxplot(outlier.shape = NA, width=0.2, color="grey80" )+
  theme_classic()+
  scale_fill_manual(values=scales::alpha(c("#E40303", "salmon", "#FF8C00", "darkmagenta", "#24408E"),0.6))+
  ylim(0,1)+
  geom_hline(aes(yintercept=quantile(res_df$turnover, na.rm = T)[3]), linewidth=0.8, color=scales::alpha("grey20",1), linetype="dashed")+
  ylab(expression(paste( "Species turnover",  sep="")))+
  xlab(expression(paste( "Interaction strength ",  sep="")))+
  theme(legend.position="none", axis.title = element_text(size=17), 
        axis.text = element_text(size=12)  )
panelE

library(cowplot)
library(grid)
library(gridExtra)

mygrid <- plot_grid( panelA, panelB, panelC, panelD, panelE,
                     labels = c("a", "b", "c", "d", "e"), 
                     nrow=1, align="hv", label_size = 25, 
                     label_fontface = "plain" )
mygrid

save_plot("~/Fig4.pdf", mygrid,
          ncol = 5, 
          nrow = 1, 
          base_aspect_ratio = 0.9)
