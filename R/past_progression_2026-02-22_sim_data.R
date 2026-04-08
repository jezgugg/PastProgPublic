# Power calculation after sample selection
# ----------------------------------------

rm(list=ls())

save_path         <- "C:/past_grogression_"
save_version      <- 3

############################################

library(dplyr)
library(plyr)
library(rstatix)
library(ggplot2)
library(ggtext)

############################################

rowVars <- function(x, na.rm=F) {
    rowSums((x - rowMeans(x, na.rm=na.rm))^2, na.rm=na.rm) / (ncol(x) - 1)
}

theme_grs <- function (base_size = 12, base_family = "") {
    theme_gray(base_size = base_size, base_family = base_family) %+replace% 
        theme(
            axis.text = element_text(colour = "black"),
            axis.title.x = element_text(colour = "black", size=rel(1), margin = margin(t = 5)),
            axis.title.y = element_text(colour = "black", size=rel(1), angle=90, margin = margin(r = 5)),
            axis.text.x = element_text(colour = "black", size=rel(0.8)),
            axis.text.y = element_text(colour = "black", size=rel(0.8), margin = margin(r = 5)),
            plot.title = element_text(hjust = 0.5, size=rel(1.1), margin=margin(t = 5, b = 5)),
            panel.grid.minor = element_blank(), 
            panel.grid.major = element_blank(),
            plot.background = element_blank(),
            panel.background = element_rect(fill="white"),
            legend.text = element_text(size=rel(0.9)),
            legend.title = element_text(size=rel(0.9), hjust=0),
            legend.key.height = unit(0.35, 'cm'),
            plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"),
    )   
}

############################################

# Simulations (empirical)
# -----------------------

num_iterations         <- 1000
num_subs               <- 1000
cor_ser_types          <- c(0.000, 0.125, 0.250, 0.375, 0.500)
prog1_ser_mean         <- (-0.45)
prog1_ser_sd           <- (-0.45)
prog2_ser_mean         <- (-0.40)
prog2_ser_sd           <- (-0.40)
min_subs_per_age       <- 160
num_in_rct             <- 80
select_criteria        <- c("unselected","select_on_ser")
effect_types           <- c("treatment_perc","treatment_absolute")
responder_types        <- c("all","half","mixed")
effect_perc            <- 0.3
effect_abs_ser         <- (-0.13)
num_cors               <- length(cor_ser_types)
num_types              <- length(effect_types)
num_crieria            <- length(select_criteria)
num_responders         <- length(responder_types)

results                <- as.data.frame(matrix(nrow=1, ncol=14))
names(results)         <- c("Selection_trait","Corr_pre_vs_post","Threshold","Responder","Iteration","Obs_corPrePost",
                            "sd_deltaSER_PRE","sd_deltaSER_POST","sd_new_deltaSER_POST",
                            "deltaSER_PRE","deltaSER_POST","new_deltaSER_POST",
                            "Delta","Pval")
myrow                  <- 1

pb = txtProgressBar(min = 0, max = num_iterations, initial = 0) 

for(myit  in 1:num_iterations){
  setTxtProgressBar(pb,myit)
for(mycor in 1:num_cors){
for(mytype in 1:num_types){
for(myresp in 1:num_responders){
for(mycrit in 1:num_crieria){
  currCor              <- cor_ser_types[mycor]
  currCriterion        <- select_criteria[mycrit]
  currType             <- effect_types[mytype]
  currResp             <- responder_types[myresp]
  mymat                <- matrix(cbind(1.000, currCor, currCor, 1.000),nrow=2)
  rownames(mymat)      <- c("SER1","SER2")
  colnames(mymat)      <- c("SER1","SER2")
  U                    <- t(chol(mymat))
  nvars                <- 2
  r.norm               <- matrix(rnorm(nvars*num_subs, mean=0,sd=1), nrow=nvars, ncol=num_subs);
  UU                   <- U %*% r.norm
  newU                 <- t(UU)
  data5                <- as.data.frame(newU)
  data5$deltaSER_PRE   <- (data5$SER1*prog1_ser_sd) + prog1_ser_mean
  data5$deltaSER_POST  <- (data5$SER2*prog2_ser_sd) + prog2_ser_mean
  data5$DoseConstant   <- 1
  PRE_STATS            <- as.data.frame(data5 %>% 
                          dplyr::summarise(meanSER   = mean(deltaSER_PRE), 
                                           medianSER = median(deltaSER_PRE)))
  data7                <- data5
  data7$SERFastProg    <- ifelse(data7$deltaSER_PRE < PRE_STATS$medianSER,1,0)
  mydataA              <- data7
  if(currCriterion=="unselected"){ mydataF <- mydataA }
  if(currCriterion=="select_on_ser"){ mydataF <- mydataA[which(mydataA$SERFastProg==1),] }
  mysampleF            <- mydataF[sample(nrow(mydataF), num_in_rct, replace = TRUE), ]
  mysampleF$treated    <- rbinom(n=num_in_rct, size=1, prob=0.5)
  mysampleF$new_deltaSER_POST <- mysampleF$deltaSER_POST
  if(currResp=="all"){
    if(currType=="treatment_absolute"){ mysampleF[which(mysampleF$treated==1),]$new_deltaSER_POST <- (mysampleF[which(mysampleF$treated==1),]$deltaSER_POST - effect_abs_ser) }
    if(currType=="treatment_perc")    { mysampleF[which(mysampleF$treated==1),]$new_deltaSER_POST <- (mysampleF[which(mysampleF$treated==1),]$deltaSER_POST * (1 - effect_perc)) }
  }
  if(currResp=="half"){
    mysampleF$responder <- rbinom(n=num_in_rct, size=1, prob=0.5)
    if(currType=="treatment_absolute"){ mysampleF[which(mysampleF$treated==1 & mysampleF$responder),]$new_deltaSER_POST <- (mysampleF[which(mysampleF$treated==1 & mysampleF$responder),]$deltaSER_POST - (2*effect_abs_ser)) }
    if(currType=="treatment_perc")    { mysampleF[which(mysampleF$treated==1 & mysampleF$responder),]$new_deltaSER_POST <- (mysampleF[which(mysampleF$treated==1 & mysampleF$responder),]$deltaSER_POST * (1 - (2*effect_perc))) }
  }
  if(currResp=="mixed"){
    mysampleF$responder <- runif(n=num_in_rct, min=0, max=2)
    mysampleF$eff1      <- mysampleF$responder*effect_abs_ser
    mysampleF$eff2      <- mysampleF$deltaSER_POST - mysampleF$eff1
    if(currType=="treatment_absolute"){ mysampleF[which(mysampleF$treated==1),]$new_deltaSER_POST <- (mysampleF[which(mysampleF$treated==1),]$eff2) }
    mysampleF$eff1      <- mysampleF$responder*effect_perc
    mysampleF$eff2      <- (1 -  mysampleF$eff1)
    mysampleF$eff3      <- mysampleF$deltaSER_POST*mysampleF$eff2
    if(currType=="treatment_perc")    { mysampleF[which(mysampleF$treated==1),]$new_deltaSER_POST <- (mysampleF[which(mysampleF$treated==1),]$eff3) }
  }
  ser_ttF              <- t.test(x=mysampleF[which(mysampleF$treated==1),]$new_deltaSER_POST,
                                 y=mysampleF[which(mysampleF$treated==0),]$new_deltaSER_POST,
                                 alternative = "greater", paired = FALSE, var.equal = FALSE, conf.level = 0.95)
  
  results[myrow,]$Selection_trait      <- currCriterion
  results[myrow,]$Corr_pre_vs_post     <- currCor
  results[myrow,]$Threshold            <- currType
  results[myrow,]$Responder            <- currResp
  results[myrow,]$Iteration            <- myit
  results[myrow,]$Obs_corPrePost       <- cor(mysampleF$deltaSER_PRE,mysampleF$deltaSER_POST)
  results[myrow,]$deltaSER_PRE         <- mean(mysampleF$deltaSER_PRE)
  results[myrow,]$deltaSER_POST        <- mean(mysampleF$deltaSER_POST)
  results[myrow,]$sd_new_deltaSER_POST <- sd(mysampleF$new_deltaSER_POST)
  results[myrow,]$sd_deltaSER_PRE      <- sd(mysampleF$deltaSER_PRE)
  results[myrow,]$sd_deltaSER_POST     <- sd(mysampleF$deltaSER_POST)
  results[myrow,]$Delta                <- as.numeric(ser_ttF$estimate[1] - ser_ttF$estimate[2])
  results[myrow,]$Pval                 <- ser_ttF$p.value
  myrow                                <- myrow + 1
  }
  }
  }
  }
  }
close(pb)
results$signif        <- ifelse(results$Pval < 0.05, 1, 0)
res_sum               <- as.data.frame(results %>% group_by(Selection_trait,Corr_pre_vs_post,Threshold,Responder) %>% 
                          get_summary_stats(signif, type="mean_ci"))
res_sum
mylab                     <- paste0((effect_perc*100),"% slowing")
res_sum$Threshold         <- ordered(res_sum$Threshold, levels=c("treatment_absolute","treatment_perc"))
res_sum$Threshold         <- revalue(res_sum$Threshold, c("treatment_absolute"="Absolute threshold", "treatment_perc"=mylab)) 
res_sum$Selection_trait   <- ordered(res_sum$Selection_trait, levels=c("unselected","select_on_ser"))
res_sum$Selection_trait   <- revalue(res_sum$Selection_trait, c("unselected"="Unselected", "select_on_ser"="Fast SER prog."))
res_sum$Responder         <- ordered(res_sum$Responder, levels=c("all","half","mixed"))
res_sum$Responder         <- revalue(res_sum$Responder, c("all"="All respond equally", "half"="Half respond; Half do not",
                                                                                           "mixed"="All respond, but variably")) 

theme_set(theme_grs())
plot1 <- ggplot(res_sum, aes(x=Corr_pre_vs_post,y=mean,colour=Selection_trait))+
  theme(panel.grid.major = element_line(linetype="dotted", colour = "grey", linewidth = 0.5)) +
  theme(axis.title.x = element_markdown())+
  coord_cartesian(ylim=c(0,1))+
  labs(x="Correlation between &Delta;SER<sub>PRE</sub> and &Delta;SER<sub>POST</sub> (*&rho;*<sub>SER</sub>)",y="Statistical power")+
  scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1),labels = scales::percent)+
  scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5))+
  scale_colour_manual(values=c("dark red","blue"),name="Selection scheme")+
  geom_point(size=2, position=position_dodge(width=0.05))+
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), position=position_dodge(width=0.05), linewidth=0.6, width=0)+
  theme(strip.background = element_rect(fill="#333333", colour="#333333"))+
  theme(strip.text = element_text(colour="white", face="bold"))+
  theme(legend.position = c(0.18,0.9))+
  facet_grid(Threshold ~ Responder)
plot1
 
out_file=paste0(save_path,"fig5_v",save_version,".tiff")
ggsave(out_file, units="cm", width=18, height=11, dpi=300, compression = 'lzw')
dev.off()

out_file=paste0(save_path,"simulation_results_v",save_version,".csv")
write.csv(results, file=out_file, row.names=FALSE)














