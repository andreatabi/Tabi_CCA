rm(list=ls())
dev.off()

library(dplyr)
library(ggplot2)
library(ppcor)
library(pcalg)
library(ggdag)      
library(dagitty)    
library(discharge)
library(tseriesChaos)

## LOAD DATA  
load("~/emp.RData")  
source("~/toolbox.R")     

df <- lapply(seq(length(db)), function(x) db[[x]][[2]] )
res_df <- do.call("rbind", df)

df2 <- lapply(seq(length(db)), function(x) db[[x]][[1]] )
res_df2 <- do.call("rbind", df2)

## FIGURE 2

TO <- ggplot(res_df) + 
  geom_line(aes(y=turnover, x=time2, by=Site_code), color=alpha("#24408E", 0.5) ) + 
  geom_smooth(aes(y=turnover, x=time2), method="lm", color="#FF8C00", size=1 )+
  theme_classic()+
  theme(legend.position="none", axis.title = element_text(size=18), 
        axis.text = element_text(size=12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1) )+
  ylab("Species turnover")+
  xlab("Time")
TO

SR <- ggplot(res_df) + 
  geom_line(aes(y=sr/sp_pool, x=time2, by=Site_code), color=alpha("#24408E", 0.5) ) + 
  geom_smooth(aes(y=sr/sp_pool, x=time2), method="lm", color="#FF8C00", size=1 )+
  theme_classic()+
  theme(legend.position="none", axis.title = element_text(size=18), 
        axis.text= element_text(size=12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1) )+
  ylab("Relative richness")+
  xlab("Time")
SR

EV <- ggplot(res_df) + 
  geom_line(aes(y=even, x=time2, by=Site_code), color=alpha("#24408E", 0.5) ) + 
  geom_smooth(aes(y=even, x=time2), method="lm", color="#FF8C00", size=1 )+
  theme_classic()+
  theme(legend.position="none", axis.title = element_text(size=18), 
        axis.text = element_text(size=12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1) )+
  ylab("Species evenness")+
  xlab("Time")
EV

m4 <- data.frame(m=median(res_df2$cor , na.rm=T))
hist_GL <- ggplot(res_df2) + 
  geom_histogram(aes(x=cor), fill=alpha("#24408E",0.6),  color="white", bins = 20) + 
  theme_classic()+
  geom_vline(aes(xintercept=median(res_df2$cor, na.rm=T) ), size=1, color="#FF8C00", linetype="dashed")+
  xlab(expression(paste( "Spearman's coefficient ")))+
  ylab("Sampling sites")+  
  theme(axis.title = element_text(size=18), axis.text = element_text(size=12)  )
hist_GL

m3 <- data.frame(m=median(res_df2$mi , na.rm=T))
hist_mi <- ggplot(res_df2) + 
  geom_histogram(aes(x=mi), fill=alpha("#24408E",0.6),  color="white", bins = 20) + 
  theme_classic()+
  geom_vline(aes(xintercept=median(res_df2$mi, na.rm=T) ), size=1, color="#FF8C00", linetype="dashed")+
  xlab(expression(paste( "Mutual information")))+
  ylab("Sampling sites")+
  theme(axis.title = element_text(size=18), axis.text = element_text(size=12)  )
hist_mi

mygrid <- plot_grid( SR,  EV,  TO, hist_GL, hist_mi,
                     labels = c("a", "b", "c", "d", "e"), 
                     nrow=1, align="hv", label_size = 25, # hjust=-4, 
                     label_fontface = "plain" )

save_plot("~/Fig2.pdf", mygrid,
          ncol = 5, 
          nrow = 1, 
          base_aspect_ratio = 0.9)

## CAUSAL DISCOVERY ------------------------------------------------------------------

sites <- sort(unique(res_df$Site_code))
res0 <- lapply(sites, function(x){ 
  temp <- res_df %>% filter(Site_code==x)    
  temp_t1 <- temp[-1,-2] 
  colnames(temp_t1)[-c(1)] <- c("time2",  paste0(colnames(temp_t1)[-c(1:2)], "_t1"))
  colnames(temp)[-c(1)] <- c("time2",  "time1", paste0(colnames(temp)[-c(1:3)], "_t"))
  res <- merge(temp, temp_t1, by=c("Site_code", "time2"), all = T)
  return(res)
})

res <- do.call("rbind", res0)

## Figure S1
pdf("~/FigS1.pdf")
par(mfrow=c(3,3))
acf(na.omit(res$losses_t), type = "correlation", main="Losses")
acf(na.omit(res$gains_t), type = "correlation", main = "Gains")
acf(na.omit(res$noise_t), type = "correlation", main="Noise")
acf(na.omit(res$pred_t), type = "correlation", main="Predation")
acf(na.omit(res$disp_t), type = "correlation", main="Dispersal")
acf(na.omit(res$volt_t), type = "correlation", main="Voltinism")
acf(na.omit(res$size_t), type = "correlation", main="Size")
dev.off()

## SHD analysis for causal discovery

a_vector <- c(1e-05, 5e-05, 1e-04, 5e-04, 0.001, 0.005, 0.01)
n <- 5000
shd_res <- lapply(a_vector, function(y){   
    shd.av <- lapply(seq(n), function(x){
    var <- c(  "gains_t", "losses_t", "pred_t", "noise_t", "disp_t", "volt_t", "size_t",
                 "gains_t1", "losses_t1" , "pred_t1",  "disp_t1" , "volt_t1", "size_t1")
    dat <- res[,sample(var)] %>% na.omit() 
    lab <- colnames(dat)
    suffStat <- list(C = cor(dat, method = "spearman"), n = nrow(dat))
    pc.g <- pc(suffStat,
               indepTest = gaussCItest, labels = lab, 
               alpha = y, skel.method = "stable",
               u2pd = "retry")
    dat2 <- res[,sample(var)] %>% na.omit() 
    lab2 <- colnames(dat2)
    suffStat2 <- list(C = cor(dat2, method = "spearman"), n = nrow(dat2))
    pc.g2 <- pc(suffStat2,
                indepTest = gaussCItest, labels = lab2, 
                alpha = y, skel.method = "stable",
                u2pd = "retry")
    
    shdist(pc.g, pc.g2)
  })
  out <- unlist(shd.av)
})

df.shd <- do.call("cbind", shd_res)
colnames(df.shd) <- a_vector
df.shd2 <- reshape2::melt(df.shd)
df.shd3 <- df.shd2 %>% group_by(Var2) %>% summarise(m=mean(value),s=sd(value))

## Figure S2
ggplot(df.shd3, aes(x=factor(Var2), y=m))+ 
  geom_point(color="#24408E")+
  geom_errorbar(aes(ymin=m-1.96*s, ymax=m+1.96*s), width=.2, color="#24408E") +
  theme_classic()+
  ylab("SHD")+xlab(expression(paste(alpha)))+
  theme(legend.position="none", axis.title = element_text(size=18), 
        axis.text = element_text(size=12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1) )

ggsave("~/FigS2.pdf", width=5, height=4)

## causal graph 
alpha <-df.shd3$Var2[which(df.shd3$m == min(df.shd3$m))]
var <- c(  "gains_t", "losses_t", "pred_t", "noise_t", "disp_t", "volt_t", "size_t",
            "gains_t1", "losses_t1" , "pred_t1",  "disp_t1" , "volt_t1", "size_t1")

dat <- res[,var] %>% na.omit() 
lab <- colnames(dat)
suffStat <- list(C = cor(dat, method = "spearman"), n = nrow(dat))
pc.g <- pc(suffStat,
           indepTest = gaussCItest, labels = lab, #verbose = T,
           alpha = alpha, skel.method = "stable",
           u2pd = "retry")
# require(Rgraphviz)
# plot(pc.g, main = "CPDAG")
# dev.off()

n <- nrow(dat)
V <- lab 
amat <- as(pc.g, "amat")    
amat[amat != 0] <- 1
if(requireNamespace("dagitty")) {
  gpc <- pcalg2dagitty(amat,V,type="dag")
  dagitty::is.dagitty(gpc)
}
coordinates(gpc) <- list( x=c(gains_t=1, gains_t1=2.5, losses_t=1.5, losses_t1=3, pred_t=1, pred_t1=2.5, 
                              disp_t=1.5, disp_t1=3, volt_t=1.5, volt_t1=2.5, noise_t=2, size_t=1, size_t1=3),
                          y=c(gains_t=1, gains_t1=1, losses_t=2, losses_t1=2, pred_t=3, pred_t1=3, 
                              disp_t=4, disp_t1=4, volt_t=5, volt_t1=5, noise_t=6, size_t=7, size_t1=7) )

## Figure S3
plot(gpc)



## Calculate structural causal effects ----------------------------------------------------
### Gains and Losses
S1 <- dagitty::adjustmentSets(gpc,"gains_t","losses_t",type = "minimal", effect = c("direct"))
S1
GL <- cor.test(dat$gains_t, dat$losses_t, method="spearman" )
GL
round(GL$est,3)

S1.2 <- dagitty::adjustmentSets(gpc,"gains_t1","losses_t1",type = "minimal", effect = c("direct"))
S1.2
GL2 <- pcor.test(dat$gains_t1, dat$losses_t1, dat$gains_t, method="spearman" )
GL2.2 <- pcor.test(dat$gains_t1, dat$losses_t1, dat$losses_t, method="spearman" )
GL2

### Predation and Losses
S2 <- dagitty::adjustmentSets(gpc,"pred_t1","losses_t1",type = "minimal", effect = c("direct"))
S2
PL <- cor.test(dat$pred_t1, dat$losses_t1, method="spearman" )
PL
round(PL$est,2)

### Predation and size
S4 <- dagitty::adjustmentSets(gpc,"pred_t","size_t",type = "minimal", effect = c("direct"))
S4
PS <- cor.test(dat$pred_t, dat$size_t, method="spearman" )
PS$est
round(PS$est,2)
S4.2 <- dagitty::adjustmentSets(gpc,"pred_t1","size_t1",type = "minimal", effect = c("direct"))
S4.2
PS2 <- pcor.test(dat$pred_t1, dat$size_t1, dat[,c("disp_t1", "size_t")], method="spearman" )
round(PS2$est,2)

### Noise and voltinism
S5 <- dagitty::adjustmentSets(gpc,"noise_t","volt_t1",type = "minimal", effect = c("direct"))
S5
NV <- cor.test(dat$noise_t, dat$volt_t1, method="spearman" )
NV
round(NV$est,2)

### Dispersal and size
S6 <- dagitty::adjustmentSets(gpc,"disp_t","size_t",type = "minimal", effect = c("direct"))
S6
DS <- cor.test(dat$disp_t, dat$size_t,  method="spearman" )
DS
round(DS$est,2)

S6.2 <- dagitty::adjustmentSets(gpc,"disp_t1","size_t1",type = "minimal", effect = c("direct"))
S6.2
DS2 <- pcor.test(dat$disp_t1, dat$size_t1, dat[,c("pred_t1", "size_t")],  method="spearman" )
DS2 <- pcor.test(dat$disp_t1, dat$size_t1, dat[,c("pred_t", "size_t")],  method="spearman" )
DS2 <- pcor.test(dat$disp_t1, dat$size_t1, dat[,c("disp_t")],  method="spearman" )
DS2

### Size and voltanism
S7 <- dagitty::adjustmentSets(gpc,"size_t","volt_t",type = "minimal", effect = c("direct"))
S7
SV <- cor.test(dat$size_t, dat$volt_t,  method="spearman" )
SV
round(SV$est,2)

S7.2 <- dagitty::adjustmentSets(gpc,"size_t1","volt_t1",type = "minimal", effect = c("direct"))
S7.2
SV2 <- pcor.test(dat$size_t1, dat$volt_t1, dat$volt_t,  method="spearman" )
SV2 <- pcor.test(dat$size_t1, dat$volt_t1, dat$size_t,  method="spearman" )
SV2

