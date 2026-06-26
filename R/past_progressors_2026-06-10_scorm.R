# Power of a myopia control RCT after selecting for fast progressors
# ------------------------------------------------------------------

# Data from SCORM
# ===============

rm(list=ls())

save_path         <- "C:/"
save_version      <- 3
data_file         <- "C:/changshiVisit_1-4.dta"

############################################

library(dplyr)
library(plyr)
library(rstatix)
library(ggplot2)
library(cowplot)
library(ggtext)
library(ggh4x)
library(haven)
library(data.table)

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
theme_set(theme_grs())

scaleFUN2        <- function(x) sprintf("%.2f", x)

############################################

# Data preparation (SCORM)
# ------------------------

data1 <- as.data.frame(read_dta(data_file))

names(data1)[which(names(data1) == "sno")] <- "SubjectID"

# Convert Stata Gender variable to plain R character values: "F" / "M"
# Assuming original coding is 1 = male, 2 = female
data1$Gender <- ifelse(as.numeric(data1$Gender) == 2, "F", "M")

table(data1$Gender, useNA = "ifany")

# Map RIGHT-eye SCORM variables to the names the script expects
data1$lsemean1 <- data1$rsemean1
data1$lsemean2 <- data1$rsemean2
data1$lsemean3 <- data1$rsemean3
data1$lsemean4 <- data1$rsemean4
data1$lal1     <- data1$ral1
data1$lal2     <- data1$ral2
data1$lal3     <- data1$ral3
data1$lal4     <- data1$ral4

# Years 1-3 (PRE = Years 1-2; POST = Years 2-3)
data1$BaselineAge      <- data1$age + 1
data1$BaselineSER      <- data1$lsemean2
data1$BaselineAXL      <- data1$lal2
data1$deltaSER_PRE     <- data1$lsemean2 - data1$lsemean1
data1$deltaSER_POST    <- data1$lsemean3 - data1$lsemean2
data1$deltaAXL_PRE     <- data1$lal2 - data1$lal1
data1$deltaAXL_POST    <- data1$lal3 - data1$lal2
data1$SER1             <- data1$lsemean1
data1$DoseConstant     <- 1 # Not used for SCORM
data2a                 <- data1[which(data1$SER1 < -0.50),
                          c("SubjectID","Gender","BaselineAge","BaselineSER","BaselineAXL","DoseConstant","deltaSER_PRE","deltaSER_POST","deltaAXL_PRE","deltaAXL_POST")]

# Years 2-4 (PRE = Years 2-3; POST = Years 3-4)
data1$BaselineAge      <- data1$age + 2
data1$BaselineSER      <- data1$lsemean3
data1$BaselineAXL      <- data1$lal3
data1$deltaSER_PRE     <- data1$lsemean3 - data1$lsemean2
data1$deltaSER_POST    <- data1$lsemean4 - data1$lsemean3
data1$deltaAXL_PRE     <- data1$lal3 - data1$lal2
data1$deltaAXL_POST    <- data1$lal4 - data1$lal3
data1$SER1             <- data1$lsemean2
data1$DoseConstant     <- 1 # Not used for SCORM
data2b                 <- data1[which(data1$SER1 < -0.50),
                          c("SubjectID","Gender","BaselineAge","BaselineSER","BaselineAXL","DoseConstant","deltaSER_PRE","deltaSER_POST","deltaAXL_PRE","deltaAXL_POST")]
data3                  <- rbind(data2a,data2b)
data4                  <- data3[complete.cases(data3),]
modQC_PRE              <- lm(deltaSER_PRE  ~ deltaAXL_PRE  + BaselineSER + BaselineAXL + BaselineAge + Gender, data=data4) 
modQC_POST             <- lm(deltaSER_POST ~ deltaAXL_POST + BaselineSER + BaselineAXL + BaselineAge + Gender, data=data4)
data4$resQC_PRE        <- modQC_PRE$residuals
data4$resQC_POST       <- modQC_POST$residuals
data4$outlierPRE       <- ifelse((data4$resQC_PRE  > -1 & data4$resQC_PRE  < 1),0,1)
data4$outlierPOST      <- ifelse((data4$resQC_POST > -1 & data4$resQC_POST < 1),0,1)
data5                  <- data4[which(data4$outlierPRE==0 & data4$outlierPOST==0),]
plot1                  <- ggplot(data4, aes(deltaSER_PRE, deltaAXL_PRE, fill=factor(outlierPRE)))+
                          theme(axis.title.x = element_markdown())+
                          theme(axis.title.y = element_markdown())+
                          theme(legend.position = c(0.2,0.2))+
                          labs(x="&Delta;SER<sub>PRE</sub> (D/year)",y="&Delta;AL<sub>PRE</sub> (mm/year)")+
                          geom_hline(aes(yintercept=0),linetype="dashed", colour="grey")+
                          geom_vline(aes(xintercept=0),linetype="dashed", colour="grey")+
                          scale_y_continuous(breaks=c(-3,-2,-1,0,1,2,3), labels=scaleFUN2)+
                          scale_x_continuous(breaks=c(-4,-3,-2,-1,0,1,2,3), labels=scaleFUN2)+
                          scale_fill_manual(name="Outlier point", values=c("white","red"),labels=c("No","Yes"))+
                          geom_point(size=2, shape=21, colour="dark blue")
plot1
plot2                  <- ggplot(data4, aes(deltaSER_POST, deltaAXL_POST, fill=factor(outlierPOST)))+
                          theme(axis.title.x = element_markdown())+
                          theme(axis.title.y = element_markdown())+
                          theme(legend.position = c(0.2,0.2))+
                          labs(x="&Delta;SER<sub>POST</sub> (D/year)",y="&Delta;AL<sub>POST</sub> (mm/year)")+
                          geom_hline(aes(yintercept=0),linetype="dashed", colour="grey")+
                          geom_vline(aes(xintercept=0),linetype="dashed", colour="grey")+
                          scale_y_continuous(breaks=c(-3,-2,-1,0,1,2,3), labels=scaleFUN2)+
                          scale_x_continuous(breaks=c(-4,-3,-2,-1,0,1,2,3), labels=scaleFUN2)+
                          #scale_colour_manual(name="Outlier point", values=c("dark blue","red"),labels=c("No","Yes") )+
                          scale_fill_manual(name="Outlier point", values=c("white","red"),labels=c("No","Yes") )+
                          geom_point(size=2, shape=21, colour="dark blue")
plot2

out_file=paste0(save_path, "outliers_v", save_version, ".tiff")
tiff(out_file, width = 18, height = 12, units = "cm", compression = "lzw", res=300)
print(plot_grid(plot1, plot2, ncol=2, rel_widths=c(1,1), labels=c("A","B")))
dev.off()


# Simulations (empirical)
# -----------------------

num_iterations         <- 1000 
min_subs_per_age       <- 160
num_in_rct             <- round(min_subs_per_age/2)
select_criteria        <- c("unselected","select_on_ser","select_on_axl")
effect_types           <- c("treatment_perc","treatment_absolute")
responder_types        <- c("all","half","mixed")
effect_perc            <- 0.3
effect_abs_ser         <- (-0.13)
effect_abs_axl         <- (0.06)
num_types              <- length(effect_types)
num_crieria            <- length(select_criteria)
num_responders         <- length(responder_types)

mytab                  <- table(data4$BaselineAge)
mytab2                 <- mytab[which(mytab > min_subs_per_age)]
data6                  <- data5[which(data5$BaselineAge %in% unlist(dimnames(mytab2))),]
PRE_STATS              <- as.data.frame(data6 %>% group_by(BaselineAge) %>% 
                          dplyr::summarise(meanSERpre   = mean(deltaSER_PRE), 
                                           sdSERpre     = sd(deltaSER_PRE),
                                           meanSERpost  = mean(deltaSER_POST),
                                           sdSERpost    = sd(deltaSER_POST),
                                           medianSERpre = median(deltaSER_PRE),
                                           meanAXLpre   = mean(deltaAXL_PRE),
                                           sdAXLpre     = sd(deltaAXL_PRE),
                                           meanAXLpost  = mean(deltaAXL_POST),
                                           sdAXLpost    = sd(deltaAXL_POST),
                                           medianAXLpre = median(deltaAXL_PRE)))

num_ages               <- nrow(PRE_STATS)
data7                  <- merge(data6,PRE_STATS[,c("BaselineAge","medianSERpre","medianAXLpre")], by="BaselineAge")
data7$SERFastProg      <- ifelse(data7$deltaSER_PRE < data7$medianSERpre,1,0)
data7$AXLFastProg      <- ifelse(data7$deltaAXL_PRE > data7$medianAXLpre,1,0)

demog1                 <- as.data.frame(data7 %>% group_by(BaselineAge) %>% 
                          dplyr::summarise(N               = n(),
                                           femaleN         = n()*mean(Gender=="F"),
                                           femalePerc      = 100*mean(Gender=="F"),
                                           meanSER_PRE     = mean(deltaSER_PRE), 
                                           sdSER_PRE       = sd(deltaSER_PRE),
                                           meanSER_POST    = mean(deltaSER_POST),
                                           sdSER_POST      = sd(deltaSER_POST),
                                           meanAXL_PRE     = mean(deltaAXL_PRE),
                                           sdAXL_PRE       = sd(deltaAXL_PRE),
                                           meanAXL_POST    = mean(deltaAXL_POST),
                                           sdAXL_POST      = sd(deltaAXL_POST),
                                           corSER_PRE_POST = cor(deltaSER_PRE,deltaSER_POST),
                                           corSER_pval     = cor.test(deltaSER_PRE,deltaSER_POST)$p.value,
                                           corAXL_PRE_POST = cor(deltaAXL_PRE,deltaAXL_POST),
                                           corAXL_pval     = cor.test(deltaAXL_PRE,deltaAXL_POST)$p.value))

demog2                 <- as.data.frame(data7 %>%  
                          dplyr::summarise(N               = n(),
                                           femaleN         = n()*mean(Gender=="F"),
                                           femalePerc      = 100*mean(Gender=="F"),
                                           meanSER_PRE     = mean(deltaSER_PRE), 
                                           sdSER_PRE       = sd(deltaSER_PRE),
                                           meanSER_POST    = mean(deltaSER_POST),
                                           sdSER_POST      = sd(deltaSER_POST),
                                           meanAXL_PRE     = mean(deltaAXL_PRE),
                                           sdAXL_PRE       = sd(deltaAXL_PRE),
                                           meanAXL_POST    = mean(deltaAXL_POST),
                                           sdAXL_POST      = sd(deltaAXL_POST),
                                           corSER_PRE_POST = cor(deltaSER_PRE,deltaSER_POST),
                                           corSER_pval     = cor.test(deltaSER_PRE,deltaSER_POST)$p.value,
                                           corAXL_PRE_POST = cor(deltaAXL_PRE,deltaAXL_POST),
                                           corAXL_pval     = cor.test(deltaAXL_PRE,deltaAXL_POST)$p.value))

demog2$BaselineAge       <- NA
demog                    <- rbind(demog1,demog2)
demog[,c(5:13,15)]       <- lapply(demog[,c(5:13,15)], sprintf, fmt = "%.3f")
demog[,c(14,16)]         <- lapply(demog[,c(14,16)], sprintf, fmt = "%.2e")
demog[,4]                <- round(as.numeric(demog[,4]),1)

out_file=paste0(save_path,"table2_v",save_version,".csv")
write.csv(demog, file=out_file, row.names=FALSE)

results                <- as.data.frame(matrix(nrow=1, ncol=21))
names(results)         <- c("Selection_trait","RCT_trait","Threshold","Responder","BaselineAge","Iteration","corPrePost",
                            "sd_deltaSER_PRE","sd_deltaSER_POST","sd_new_deltaSER_POST",
                            "sd_deltaAXL_PRE","sd_deltaAXL_POST","sd_new_deltaAXL_POST",
                            "deltaSER_PRE","deltaSER_POST","new_deltaSER_POST",
                            "deltaAXL_PRE","deltaAXL_POST","new_deltaAXL_POST",
                            "Delta","Pval")
myrow                  <- 1

pb = txtProgressBar(min = 0, max = num_iterations, initial = 0) 

for(myit  in 1:num_iterations){
  setTxtProgressBar(pb,myit)
for(myage in 1:num_ages){
for(mytype in 1:num_types){
for(myresp in 1:num_responders){
for(mycrit in 1:num_crieria){
  currCriterion        <- select_criteria[mycrit]
  currType             <- effect_types[mytype]
  currAge              <- PRE_STATS$BaselineAge[myage]
  currResp             <- responder_types[myresp]
  mydataA              <- data7[which(data7$BaselineAge==currAge),]
  if(currCriterion=="unselected"){ mydataF <- mydataA }
  if(currCriterion=="select_on_ser"){ mydataF <- mydataA[which(mydataA$SERFastProg==1),] }
  if(currCriterion=="select_on_axl"){ mydataF <- mydataA[which(mydataA$AXLFastProg==1),] }
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
  results[myrow,]$RCT_trait            <- "RCT_for_SER"
  results[myrow,]$Threshold            <- currType
  results[myrow,]$BaselineAge          <- currAge
  results[myrow,]$Responder            <- currResp
  results[myrow,]$Iteration            <- myit
  results[myrow,]$corPrePost           <- cor(mysampleF$deltaSER_PRE,mysampleF$deltaSER_POST)
  results[myrow,]$deltaSER_PRE         <- mean(mysampleF$deltaSER_PRE)
  results[myrow,]$deltaSER_POST        <- mean(mysampleF$deltaSER_POST)
  results[myrow,]$deltaAXL_PRE         <- mean(mysampleF$deltaAXL_PRE)
  results[myrow,]$deltaAXL_POST        <- mean(mysampleF$deltaAXL_POST)
  results[myrow,]$sd_deltaSER_PRE      <- sd(mysampleF$deltaSER_PRE)
  results[myrow,]$sd_deltaSER_POST     <- sd(mysampleF$deltaSER_POST)
  results[myrow,]$sd_new_deltaSER_POST <- sd(mysampleF$new_deltaSER_POST)
  results[myrow,]$sd_deltaAXL_PRE      <- sd(mysampleF$deltaAXL_PRE)
  results[myrow,]$sd_deltaAXL_POST     <- sd(mysampleF$deltaAXL_POST)
  results[myrow,]$sd_new_deltaAXL_POST <- sd(mysampleF$new_deltaAXL_POST)
  results[myrow,]$Delta                <- as.numeric(ser_ttF$estimate[1] - ser_ttF$estimate[2])
  results[myrow,]$Pval                 <- ser_ttF$p.value
  myrow                                <- myrow + 1
  
  mydataA              <- data7[which(data7$BaselineAge==currAge),]
  if(currCriterion=="unselected"){ mydataF <- mydataA }
  if(currCriterion=="select_on_ser"){ mydataF <- mydataA[which(mydataA$SERFastProg==1),] }
  if(currCriterion=="select_on_axl"){ mydataF <- mydataA[which(mydataA$AXLFastProg==1),] }
  mysampleF            <- mydataF[sample(nrow(mydataF), num_in_rct, replace = TRUE), ]
  mysampleF$treated    <- rbinom(n=num_in_rct, size=1, prob=0.5)
  mysampleF$new_deltaAXL_POST <- mysampleF$deltaAXL_POST
  if(currResp=="all"){
    if(currType=="treatment_absolute"){ mysampleF[which(mysampleF$treated==1),]$new_deltaAXL_POST <- (mysampleF[which(mysampleF$treated==1),]$deltaAXL_POST - effect_abs_axl) }
    if(currType=="treatment_perc")    { mysampleF[which(mysampleF$treated==1),]$new_deltaAXL_POST <- (mysampleF[which(mysampleF$treated==1),]$deltaAXL_POST * (1 - effect_perc)) }
  }
  if(currResp=="half"){
    mysampleF$responder <- rbinom(n=num_in_rct, size=1, prob=0.5)
    if(currType=="treatment_absolute"){ mysampleF[which(mysampleF$treated==1 & mysampleF$responder),]$new_deltaAXL_POST <- (mysampleF[which(mysampleF$treated==1 & mysampleF$responder),]$deltaAXL_POST - (2*effect_abs_axl)) }
    if(currType=="treatment_perc")    { mysampleF[which(mysampleF$treated==1 & mysampleF$responder),]$new_deltaAXL_POST <- (mysampleF[which(mysampleF$treated==1 & mysampleF$responder),]$deltaAXL_POST * (1 - (2*effect_perc))) }
  }
  if(currResp=="mixed"){
    mysampleF$responder <- runif(n=num_in_rct, min=0, max=2)
    mysampleF$eff1      <- mysampleF$responder*effect_abs_axl
    mysampleF$eff2      <- mysampleF$deltaAXL_POST - mysampleF$eff1
    if(currType=="treatment_absolute"){ mysampleF[which(mysampleF$treated==1),]$new_deltaAXL_POST <- (mysampleF[which(mysampleF$treated==1),]$eff2) }
    mysampleF$eff1      <- mysampleF$responder*effect_perc
    mysampleF$eff2      <- (1 -  mysampleF$eff1)
    mysampleF$eff3      <- mysampleF$deltaAXL_POST*mysampleF$eff2
    if(currType=="treatment_perc")    { mysampleF[which(mysampleF$treated==1),]$new_deltaAXL_POST <- (mysampleF[which(mysampleF$treated==1),]$eff3) }
  }
  axl_ttF              <- t.test(x=mysampleF[which(mysampleF$treated==1),]$new_deltaAXL_POST,
                                 y=mysampleF[which(mysampleF$treated==0),]$new_deltaAXL_POST,
                                 alternative = "less", paired = FALSE, var.equal = FALSE, conf.level = 0.95)

  results[myrow,]$Selection_trait  <- currCriterion
  results[myrow,]$RCT_trait        <- "RCT_for_AXL"
  results[myrow,]$Threshold        <- currType
  results[myrow,]$BaselineAge      <- currAge
  results[myrow,]$Responder           <- currResp
  results[myrow,]$Iteration        <- myit
  results[myrow,]$corPrePost       <- cor(mysampleF$deltaAXL_PRE,mysampleF$deltaAXL_POST)
  results[myrow,]$deltaSER_PRE         <- mean(mysampleF$deltaSER_PRE)
  results[myrow,]$deltaSER_POST        <- mean(mysampleF$deltaSER_POST)
  results[myrow,]$deltaAXL_PRE         <- mean(mysampleF$deltaAXL_PRE)
  results[myrow,]$deltaAXL_POST        <- mean(mysampleF$deltaAXL_POST)
  results[myrow,]$sd_deltaSER_PRE      <- sd(mysampleF$deltaSER_PRE)
  results[myrow,]$sd_deltaSER_POST     <- sd(mysampleF$deltaSER_POST)
  results[myrow,]$sd_new_deltaSER_POST <- sd(mysampleF$new_deltaSER_POST)
  results[myrow,]$sd_deltaAXL_PRE      <- sd(mysampleF$deltaAXL_PRE)
  results[myrow,]$sd_deltaAXL_POST     <- sd(mysampleF$deltaAXL_POST)
  results[myrow,]$sd_new_deltaAXL_POST <- sd(mysampleF$new_deltaAXL_POST)
  results[myrow,]$Delta            <- as.numeric(axl_ttF$estimate[1] - axl_ttF$estimate[2])
  results[myrow,] $Pval            <- axl_ttF$p.value
  myrow                            <- myrow + 1
  }
  }
  }
  }
  }
close(pb)
results$signif        <- ifelse(results$Pval < 0.05, 1, 0)

out_file=paste0(save_path,"scorm_results_v",save_version,".csv")
write.csv(results, file=out_file, row.names=FALSE)

#-----------------------------------------------------

save_version           <- 3
in_file                <- paste0(save_path,"scorm_results_v",save_version,".csv")
results                <- as.data.frame(fread(in_file, header=TRUE))

save_version           <- 4
num_iterations         <- 1000
min_subs_per_age       <- 160
num_in_rct             <- round(min_subs_per_age/2)
select_criteria        <- c("unselected","select_on_ser","select_on_axl")
effect_types           <- c("treatment_perc","treatment_absolute")
responder_types        <- c("all","half","mixed")
effect_perc            <- 0.3
effect_abs_ser         <- (-0.13)
effect_abs_axl         <- (0.06)
num_types              <- length(effect_types)
num_crieria            <- length(select_criteria)
num_responders         <- length(responder_types)

res_sum                   <- as.data.frame(results %>% group_by(Selection_trait,RCT_trait,Threshold,Responder,BaselineAge) %>% 
                          get_summary_stats(signif, type="mean_ci"))
mylab                     <- paste0((effect_perc*100),"% slowing")
res_sum$RCT_trait         <- ordered(res_sum$RCT_trait, levels=c("RCT_for_AXL","RCT_for_SER"))
res_sum$RCT_trait         <- revalue(res_sum$RCT_trait, c("RCT_for_AXL"="RCT outcome trait: AL", "RCT_for_SER"="RCT outcome trait: SER"))
res_sum$Threshold         <- ordered(res_sum$Threshold, levels=c("treatment_absolute","treatment_perc"))
res_sum$Threshold         <- revalue(res_sum$Threshold, c("treatment_absolute"="Absolute threshold", "treatment_perc"=mylab)) 
res_sum$Selection_trait   <- ordered(res_sum$Selection_trait, levels=c("unselected","select_on_axl","select_on_ser"))
res_sum$Selection_trait   <- revalue(res_sum$Selection_trait, c("unselected"="Unselected", "select_on_axl"="Fast AL prog.",
                                                                                           "select_on_ser"="Fast SER prog."))
res_sum$Responder         <- ordered(res_sum$Responder, levels=c("all","half","mixed"))
res_sum$Responder         <- revalue(res_sum$Responder, c("all"="All respond equally", "half"="Half respond; Half do not",
                                                                                           "mixed"="All respond, but variably")) 

plot1 <- ggplot(res_sum, aes(x=BaselineAge,y=mean,colour=Selection_trait))+
  theme(panel.grid.major = element_line(linetype="dotted", colour = "grey", linewidth = 0.5)) +
  coord_cartesian(ylim=c(0,1))+
  labs(x="Age at baseline (years)",y="Statistical power")+
  scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1),labels = scales::percent)+
  scale_colour_manual(values=c("dark red","blue","orange"),name="Selection scheme")+
  geom_point(size=2, position=position_dodge(width=0.2))+
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), position=position_dodge(width=0.2), linewidth=0.6, width=0)+
  theme(strip.background = element_rect(fill="#333333", colour="#333333"))+
  theme(strip.text = element_text(colour="white", face="bold"))+
  theme(legend.position = c(0.83,0.93))+
  facet_grid(Threshold ~ RCT_trait ~ Responder)
plot1

out_file=paste0(save_path,"scorm_fig1_v",save_version,".tiff")
ggsave(out_file, units="cm", width=18, height=20, dpi=300, compression = 'lzw')
dev.off()

res_sum2               <- as.data.frame(results %>% group_by(Selection_trait,RCT_trait,Threshold,Responder,BaselineAge) %>% 
                          get_summary_stats(Delta, type="mean_ci"))
mylab                      <- paste0((effect_perc*100),"% slowing")
res_sum2$RCT_trait         <- ordered(res_sum2$RCT_trait, levels=c("RCT_for_AXL","RCT_for_SER"))
res_sum2$RCT_trait         <- revalue(res_sum2$RCT_trait, c("RCT_for_AXL"="RCT outcome trait: AL", "RCT_for_SER"="RCT outcome trait: SER"))
res_sum2$Threshold         <- ordered(res_sum2$Threshold, levels=c("treatment_absolute","treatment_perc"))
res_sum2$Threshold         <- revalue(res_sum2$Threshold, c("treatment_absolute"="Absolute threshold", "treatment_perc"=mylab)) 
res_sum2$Selection_trait   <- ordered(res_sum2$Selection_trait, levels=c("unselected","select_on_axl","select_on_ser"))
res_sum2$Selection_trait   <- revalue(res_sum2$Selection_trait, c("unselected"="Unselected", "select_on_axl"="Fast AL prog.",
                                                                                           "select_on_ser"="Fast SER prog."))
res_sum2$Responder         <- ordered(res_sum2$Responder, levels=c("all","half","mixed"))
res_sum2$Responder         <- revalue(res_sum2$Responder, c("all"="All respond equally", "half"="Half respond; Half do not",
                                                                                           "mixed"="All respond, but variably"))

plot2 <- ggplot(res_sum2, aes(x=BaselineAge,y=mean, fill=Selection_trait))+
  theme(panel.grid.major = element_line(linetype="dotted", colour = "grey", linewidth = 0.5)) +
  theme(legend.key = element_rect(colour="white")) +
  labs(x="Age at baseline (years)",y="Treatment effect (mm or D)")+
  scale_fill_manual(values=c("dark red","blue","orange"),name="Selection scheme")+
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), position=position_dodge(width=0.8), linewidth=0.5, width=0.5)+
  geom_bar(stat = "identity", position=position_dodge(width=0.8), width=0.8)+
  theme(strip.background = element_rect(fill="#333333", colour="#333333"))+
  theme(strip.text = element_text(colour="white", face="bold"))+
  theme(legend.position = "bottom")+
  facet_grid2(Threshold ~ RCT_trait ~ Responder, scales = "free_y")+
  facetted_pos_scales(
    y = list(
      scale_y_continuous(limits=c(-0.2,0.05), breaks=c(-0.2,-0.15,-0.10,-0.05,0,0.05), labels=scaleFUN2),
      scale_y_continuous(limits=c(0,0.32), breaks=c(0,0.10,0.20,0.30), labels=scaleFUN2),
      scale_y_continuous(limits=c(-0.2,0.05), breaks=c(-0.2,-0.15,-0.10,-0.05,0,0.05), labels=scaleFUN2),
      scale_y_continuous(limits=c(0,0.32), breaks=c(0,0.10,0.20,0.30), labels=scaleFUN2)))
plot2

out_file=paste0(save_path,"scorm_fig2_v",save_version,".tiff")
ggsave(out_file, units="cm", width=18, height=21, dpi=300, compression = 'lzw')
# dev.off()


res_sum3               <- as.data.frame(results %>% group_by(Selection_trait,RCT_trait,Threshold,Responder,BaselineAge) %>% 
                          get_summary_stats(sd_deltaSER_POST, type="mean_ci"))
res_sum4               <- as.data.frame(results %>% group_by(Selection_trait,RCT_trait,Threshold,Responder,BaselineAge) %>% 
                          get_summary_stats(sd_deltaAXL_POST, type="mean_ci"))
res_sum5               <- rbind(res_sum3[which(res_sum3$RCT_trait=="RCT_for_SER"),],res_sum4[which(res_sum4$RCT_trait=="RCT_for_AXL"),])
mylab                      <- paste0((effect_perc*100),"% slowing")
res_sum5$RCT_trait         <- ordered(res_sum5$RCT_trait, levels=c("RCT_for_AXL","RCT_for_SER"))
res_sum5$RCT_trait         <- revalue(res_sum5$RCT_trait, c("RCT_for_AXL"="RCT outcome trait: AL", "RCT_for_SER"="RCT outcome trait: SER"))
res_sum5$Threshold         <- ordered(res_sum5$Threshold, levels=c("treatment_absolute","treatment_perc"))
res_sum5$Threshold         <- revalue(res_sum5$Threshold, c("treatment_absolute"="Absolute threshold", "treatment_perc"=mylab)) 
res_sum5$Selection_trait   <- ordered(res_sum5$Selection_trait, levels=c("unselected","select_on_axl","select_on_ser"))
res_sum5$Selection_trait   <- revalue(res_sum5$Selection_trait, c("unselected"="Unselected", "select_on_axl"="Fast AL prog.",
                                                                                           "select_on_ser"="Fast SER prog."))
res_sum5$Responder         <- ordered(res_sum5$Responder, levels=c("all","half","mixed"))
res_sum5$Responder         <- revalue(res_sum5$Responder, c("all"="All respond equally", "half"="Half respond; Half do not",
                                                                                           "mixed"="All respond, but variably"))

plot3 <- ggplot(res_sum5, aes(x=BaselineAge,y=mean, fill=Selection_trait))+
  theme(panel.grid.major = element_line(linetype="dotted", colour = "grey", linewidth = 0.5)) +
  theme(legend.key = element_rect(colour="white")) +
  theme(axis.title.y = element_markdown()) +
  labs(x="Age at baseline (years)",
       y="Standard deviation of &Delta;SER<sub>POST</sub> or &Delta;AL<sub>POST</sub>")+
  scale_fill_manual(values=c("dark red","blue","orange"),name="Selection scheme")+
  scale_y_continuous(limits=c(0,0.5), breaks=c(0,0.1,0.2,0.3,0.4,0.5))+
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), position=position_dodge(width=0.8), linewidth=0.5, width=0.5)+
  geom_bar(stat = "identity", position=position_dodge(width=0.8), width=0.8)+
  theme(strip.background = element_rect(fill="#333333", colour="#333333"))+
  theme(strip.text = element_text(colour="white", face="bold"))+
  theme(legend.position = "bottom")+
  facet_grid2(Threshold ~ RCT_trait ~ Responder, scales = "free_y")+
  facetted_pos_scales(
    y = list(
      scale_y_continuous(limits=c(0,0.5), breaks=c(0,0.1,0.2,0.3,0.4,0.5), labels=scaleFUN2),
      scale_y_continuous(limits=c(0,0.5), breaks=c(0,0.1,0.2,0.3,0.4,0.5), labels=scaleFUN2),
      scale_y_continuous(limits=c(0,0.5), breaks=c(0,0.1,0.2,0.3,0.4,0.5), labels=scaleFUN2),
      scale_y_continuous(limits=c(0,0.5), breaks=c(0,0.1,0.2,0.3,0.4,0.5), labels=scaleFUN2)))
plot3

out_file=paste0(save_path,"fig3_v",save_version,".tiff")
ggsave(out_file, units="cm", width=18, height=21, dpi=300, compression = 'lzw')
# dev.off()


res_sum6               <- as.data.frame(results %>% group_by(Selection_trait,RCT_trait,Threshold,Responder,BaselineAge) %>% 
                          get_summary_stats(deltaSER_POST, type="mean_ci"))
res_sum7               <- as.data.frame(results %>% group_by(Selection_trait,RCT_trait,Threshold,Responder,BaselineAge) %>% 
                          get_summary_stats(deltaAXL_POST, type="mean_ci"))
res_sum8               <- rbind(res_sum6[which(res_sum6$RCT_trait=="RCT_for_SER"),],res_sum7[which(res_sum7$RCT_trait=="RCT_for_AXL"),])
mylab                      <- paste0((effect_perc*100),"% slowing")
res_sum8$RCT_trait         <- ordered(res_sum8$RCT_trait, levels=c("RCT_for_AXL","RCT_for_SER"))
res_sum8$RCT_trait         <- revalue(res_sum8$RCT_trait, c("RCT_for_AXL"="RCT outcome trait: AL", "RCT_for_SER"="RCT outcome trait: SER"))
res_sum8$Threshold         <- ordered(res_sum8$Threshold, levels=c("treatment_absolute","treatment_perc"))
res_sum8$Threshold         <- revalue(res_sum8$Threshold, c("treatment_absolute"="Absolute threshold", "treatment_perc"=mylab)) 
res_sum8$Selection_trait   <- ordered(res_sum8$Selection_trait, levels=c("unselected","select_on_axl","select_on_ser"))
res_sum8$Selection_trait   <- revalue(res_sum8$Selection_trait, c("unselected"="Unselected", "select_on_axl"="Fast AL prog.",
                                                                                           "select_on_ser"="Fast SER prog."))
res_sum8$Responder         <- ordered(res_sum8$Responder, levels=c("all","half","mixed"))
res_sum8$Responder         <- revalue(res_sum8$Responder, c("all"="All respond equally", "half"="Half respond; Half do not",
                                                                                           "mixed"="All respond, but variably"))
plot4 <- ggplot(res_sum8, aes(x=BaselineAge,y=mean, fill=Selection_trait))+
  theme(panel.grid.major = element_line(linetype="dotted", colour = "grey", linewidth = 0.5)) +
  theme(legend.key = element_rect(colour="white")) +
  theme(axis.title.y = element_markdown()) +
  labs(x="Age at baseline (years)",
       y="&Delta;SER<sub>POST</sub> or &Delta;AL<sub>POST</sub>")+
  scale_fill_manual(values=c("dark red","blue","orange"),name="Selection scheme")+
  #scale_y_continuous(limits=c(0,0.5), breaks=c(0,0.2,0.4))+
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), position=position_dodge(width=0.8), linewidth=0.5, width=0.5)+
  geom_bar(stat = "identity", position=position_dodge(width=0.8), width=0.8)+
  theme(strip.background = element_rect(fill="#333333", colour="#333333"))+
  theme(strip.text = element_text(colour="white", face="bold"))+
  theme(legend.position = "bottom")+
  facet_grid2(Threshold ~ RCT_trait ~ Responder, scales = "free_y")+
  facetted_pos_scales(
    y = list(
      scale_y_continuous(limits=c(-0.1,0.65), breaks=c(0,0.2,0.4,0.6), labels=scaleFUN2),
      scale_y_continuous(limits=c(-1.10,0.05), breaks=c(-1.0,-0.75,-0.5,-0.25,0), labels=scaleFUN2),
      scale_y_continuous(limits=c(-0.1,0.65), breaks=c(0,0.2,0.4,0.6), labels=scaleFUN2),
      scale_y_continuous(limits=c(-1.10,0.05), breaks=c(-1.0,-0.75,-0.5,-0.25,0), labels=scaleFUN2)))
plot4

out_file=paste0(save_path,"scorm_fig4_v",save_version,".tiff")
ggsave(out_file, units="cm", width=18, height=21, dpi=300, compression = 'lzw')
# dev.off()

# Create results table (Table S2)

in_file                   <- paste0(save_path,"results_v",save_version,".csv")
results                   <- read.csv(file=in_file, header=TRUE)
res_sum                   <- as.data.frame(results %>% group_by(Selection_trait,RCT_trait,Threshold,Responder,BaselineAge) %>% 
                                           get_summary_stats(signif, type="mean_ci"))
res_sum$RCT_trait         <- ordered(res_sum$RCT_trait, levels=c("RCT_for_AXL","RCT_for_SER"))
res_sum$RCT_trait         <- revalue(res_sum$RCT_trait, c("RCT_for_AXL"="AL", "RCT_for_SER"="SER"))
res_sum$Threshold         <- ordered(res_sum$Threshold, levels=c("treatment_absolute","treatment_perc"))
res_sum$Threshold         <- revalue(res_sum$Threshold, c("treatment_absolute"="Absolute", "treatment_perc"="Relative")) 
res_sum$Selection_trait   <- ordered(res_sum$Selection_trait, levels=c("unselected","select_on_axl","select_on_ser"))
res_sum$Selection_trait   <- revalue(res_sum$Selection_trait, c("unselected"="Unselected", "select_on_axl"="Fast AL prog.",
                                                                                           "select_on_ser"="Fast SER prog."))
res_sum$Responder         <- ordered(res_sum$Responder, levels=c("all","half","mixed"))
res_sum$Responder         <- revalue(res_sum$Responder, c("all"="All respond equally", "half"="Half respond; Half do not",
                                                                                           "mixed"="All respond, but variably")) 
res_sum$power      <- paste0(res_sum$mean, " (", res_sum$mean - res_sum$ci, " to ", res_sum$mean + res_sum$ci, ")")
res_sum$n          <- NULL
res_sum$mean       <- NULL
res_sum$ci         <- NULL
res_sum$variable   <- NULL

out_file=paste0(save_path,"table S2_v",save_version,".csv")
write.csv(res_sum, file=out_file, row.names=FALSE)




