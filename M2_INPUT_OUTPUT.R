# M2 metabolite accumulation:
COSVEN4_B200 <- read.csv("/Users/justmendiola/ibsc/PROJECTS/M2/FULL_BDQ_PK_COVSVEN4.csv", header = TRUE, skip = 0)

filtered_COSVEN4_B200 <- subset(COSVEN4_B200, COSVEN4_B200$TIME == 0,)
filtered_COSVEN4_700 <- subset(filtered_COSVEN4_B200, AMT == 700,)

correct_IDs <- c(filtered_COSVEN4_700$ID)
filtered_data <- subset(COSVEN4_B200, ID %in% correct_IDs)
# write.csv(filtered_data, "/Users/justmendiola/ibsc/PROJECTS/Fed's_files/Simulations/M2 METABOLITE/filtered_COSVEN4_B200.csv", row.names = FALSE)

# SIMULATION 1 WAS DONE USING THE FULL_BDQ_PK_COVSVEN3new1.csv DATASET (28/2/24)
# ARM: DAY 1: 700mg QD, DAY 2: 500mg QD, DAY 3-14: 400mg QD, DAY 15-170: 200mg TIW, (every 2 days) 
#> 5 Simulations
suppressMessages(setAs("character","dummy.numeric", function(from) as.numeric(from)))
read.table.nm <- function(x,...){
  tmp <- suppressWarnings(read.table(x, fill=T, colClasses="dummy.numeric",...))
  return(tmp[complete.cases(tmp),])}

library(ggplot2)
library(plyr)

Output.DR <- read.table.nm("/Users/justmendiola/ibsc/PROJECTS/M2/simulation_1/sdtab1", header = TRUE,stringsAsFactors=FALSE,skip=1)
#Output.DR <- Output.DR[Output.DR$TIME>=1464,]
Output.DR$DOSE <- "BDQ - 200mg/tpw"
#Output.DR <- Output.DR[Output.DR$TIME>1510,]
Output.DR$DV <- exp(Output.DR$DV)
Output.DR$IPRED <- exp(Output.DR$IPRED)
Output.DR <- Output.DR[Output.DR$MDV==0,]
N <- 5
Output.DR$ITER <- rep(1:N, each = nrow(Output.DR)/N)
sim1<-ddply(Output.DR,.(TIME,DOSE),summarise, M=median(DV), L=quantile(DV,0.05) ,H=quantile(DV, 0.95))



#Output.DR <- Output.DR[Output.DR$TIME>1510,]
PK1 <- ggplot()+
  geom_line(data=sim1,aes(TIME/24, M), 
            linewidth=1,colour="#56B4E9")+
  labs(x="Time (Days)", y="M2 concentrations (mg/L)")+
  facet_wrap(~DOSE)+
  geom_hline(aes(yintercept=0.42), col = "red", linetype = "dashed",size=1)+
  theme_bw(base_size = 18)+
  geom_text(data=data.frame(sim1),aes(130,0.425,label="Safety target - 0.42mg/L"),colour = "black",size=6)+
  geom_ribbon(data=sim1,aes(TIME/24, ymin=L,ymax=H), fill = "grey",alpha=0.6)
PK1

sim1_exceed <- subset(Output.DR, Output.DR$DV >= 0.42)
sim1_exceed<- unique(sim1_exceed$ID)
sim1_exceed_grape <- unique(sim1_exceed_grape$ID)

# SIMULATION 3 WAS DONE USING THE FULL_BDQ_PK_COVSVEN3new1.csv DATASET (15/4/24)
# ARM: DAY 1-14: 400 QD, DAY 15-170: 200mg (every 2 days) 
#> 5 Simulations

COSVEN4_WHO <-COSVEN4_B200
COSVEN4_WHO$AMT[COSVEN4_WHO$AMT >= 200 & 0 <= COSVEN4_WHO$TIME & COSVEN4_WHO$TIME <= 312] <- 400
write.csv(COSVEN4_WHO, "/Users/justmendiola/ibsc/PROJECTS/M2/COSVEN4_WHO.csv")

Output.DR3 <- read.table.nm("/Users/justmendiola/ibsc/PROJECTS/M2/simulation_3/sdtab3", header = TRUE,stringsAsFactors=FALSE,skip=1)
#Output.DR <- Output.DR[Output.DR$TIME>=1464,]
Output.DR3$DOSE <- "M2 Exposure - WHO-Recommended BDQ Regimen"
#Output.DR <- Output.DR[Output.DR$TIME>1510,]
Output.DR3$DV <- exp(Output.DR3$DV)
Output.DR3$IPRED <- exp(Output.DR3$IPRED)
Output.DR3 <- Output.DR3[Output.DR3$MDV==0,]
N <- 5
Output.DR3$ITER <- rep(1:N, each = nrow(Output.DR3)/N)
sim3<-ddply(Output.DR3,.(TIME,DOSE),summarise, M=median(DV), L=quantile(DV,0.05) ,H=quantile(DV, 0.95))

#Output.DR <- Output.DR[Output.DR$TIME>1510,]
PK3 <- ggplot()+
  geom_line(data=sim3,aes(TIME/24, M), 
            linewidth=1,colour="#56B4E9")+
  labs(x="Time (Days)", y="M2 concentrations (mg/L)")+
  facet_wrap(~DOSE)+
  geom_hline(aes(yintercept=0.42), col = "red", linetype = "dashed",size=1)+
  theme_bw(base_size = 16)+
  geom_text(data=data.frame(sim3),aes(130,0.425,label="Safety target - 0.42mg/L"),colour = "black",size=6)+
  geom_ribbon(data=sim3,aes(TIME/24, ymin=L,ymax=H), fill = "grey",alpha=0.6)
PK3

sim3_exceed <- subset(Output.DR3, Output.DR3$DV >= 0.42)
sim3_exceed<- unique(sim3_exceed$ID)
sim3_exceed_grape <- unique(sim3_exceed_grape$ID)

# SIMULATION 12 WAS DONE USING THE FULL_BDQ_PK_COVSVEN3new1.csv DATASET (15/4/24)
# ARM: DAY 1: 700mg QD, DAY 2: 500mg QD, DAY 3-14: 400mg QD, DAY 15-170: 300mg TIW, (every 2 days)
#> 5 Simulations


COSVEN4_B300 <- COSVEN4_B200
COSVEN4_B300$AMT[COSVEN4_B300$AMT == 200] <- 300
write.csv(COSVEN4_B300, "/Users/justmendiola/ibsc/PROJECTS/M2/COSVEN4_B300.csv")

Output.DR12 <- read.table.nm("/Users/justmendiola/ibsc/PROJECTS/M2/simulation_12/sdtab12", header = TRUE,stringsAsFactors=FALSE,skip=1)
#Output.DR <- Output.DR[Output.DR$TIME>=1464,]
Output.DR12$DOSE <- "BDQ - 200mg/tpw"
#Output.DR <- Output.DR[Output.DR$TIME>1510,]
Output.DR12$DV <- exp(Output.DR12$DV)
Output.DR12$IPRED <- exp(Output.DR12$IPRED)
Output.DR12 <- Output.DR12[Output.DR12$MDV==0,]
N <- 5
Output.DR12$ITER <- rep(1:N, each = nrow(Output.DR12)/N)
sim12<-ddply(Output.DR12,.(TIME,DOSE),summarise, M=median(DV), L=quantile(DV,0.05) ,H=quantile(DV, 0.95))



#Output.DR <- Output.DR[Output.DR$TIME>1510,]
PK12 <- ggplot()+
  geom_line(data=sim12,aes(TIME/24, M), 
            linewidth=1,colour="#56B4E9")+
  labs(x="Time (Days)", y="M2 concentrations (mg/L)")+
  facet_wrap(~DOSE)+
  geom_hline(aes(yintercept=0.42), col = "red", linetype = "dashed",size=1)+
  theme_bw(base_size = 18)+
  geom_text(data=data.frame(sim12),aes(130,0.425,label="Safety target - 0.42mg/L"),colour = "black",size=6)+
  geom_ribbon(data=sim12,aes(TIME/24, ymin=L,ymax=H), fill = "grey",alpha=0.6)
PK12

sim12_exceed <- subset(Output.DR12, Output.DR12$DV >= 0.42)
sim12_exceed<- unique(sim12_exceed$ID)
sim12_exceed_grape <- unique(sim12_exceed_grape$ID)

