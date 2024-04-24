# M2 metabolite accumulation:
COSVEN4_B200_M2_grape <- read.csv("/Users/justmendiola/ibsc/PROJECTS/M2_GRAPE/FULL_BDQ_PK_COVSVEN4.csv", header = TRUE, skip = 0)

filtered_COSVEN4_B200_M2_grape <- subset(COSVEN4_B200_M2_grape, COSVEN4_B200_M2_grape$TIME == 0,)
filtered_COSVEN4_700 <- subset(COSVEN4_B200_M2_grape, AMT == 700,)
unique_IDs_M2 <- unique(COSVEN4_B200_M2_grape$ID)

correct_IDs_M2 <- c(filtered_COSVEN4_700$ID)
filtered_data <- subset(COSVEN3_B200_M2, ID %in% correct_IDs_M2)
# there are 60 subjects that follow the protocol


# SIMULATION 1 WAS DONE USING THE FULL_BDQ_PK_COVSVEN3new1.csv DATASET (28/3/24)
# ARM: DAY 1: 700mg QD, DAY 2: 500mg QD, DAY 3-14: 400mg QD, DAY 15-170: 200mg TIW, (every 2 days) 
#> 5 Simulations
suppressMessages(setAs("character","dummy.numeric", function(from) as.numeric(from)))
read.table.nm <- function(x,...){
  tmp <- suppressWarnings(read.table(x, fill=T, colClasses="dummy.numeric",...))
  return(tmp[complete.cases(tmp),])}

library(ggplot2)
library(plyr)

Output.DR_grape <- read.table.nm("/Users/justmendiola/ibsc/PROJECTS/M2_GRAPE/simulation_1/sdtab1grapefruit2", header = TRUE,stringsAsFactors=FALSE,skip=1)
#Output.DR <- Output.DR[Output.DR$TIME>=1464,]
Output.DR_grape$DOSE <- "Scenario 3"
#Output.DR <- Output.DR[Output.DR$TIME>1510,]
Output.DR_grape$DV <- exp(Output.DR_grape$DV)
Output.DR_grape$IPRED <- exp(Output.DR_grape$IPRED)
Output.DR_grape <- Output.DR_grape[Output.DR_grape$MDV==0,]
# N <- 4
# Output.DR_grape$ITER <- rep(1:N, each = nrow(Output.DR_grape)/N)
sim1_grape<-ddply(Output.DR_grape,.(TIME,DOSE),summarise, M=median(DV), L=quantile(DV,0.05) ,H=quantile(DV, 0.95))



#Output.DR <- Output.DR[Output.DR$TIME>1510,]
PK1_grape <- ggplot()+
  geom_line(data=sim1_grape,aes(TIME/24, M), 
            linewidth=1,colour="#56B4E9")+
  labs(x="Time (Days)", y="M2 concentrations (mg/L)")+
  facet_wrap(~DOSE)+
  geom_hline(aes(yintercept=0.42), col = "red", linetype = "dashed",size=1)+
  theme_bw(base_size = 16)+
  geom_text(data=data.frame(sim1_grape),aes(130,0.47,label="Safety target - 0.42mg/L"),colour = "black",size=6)+
  geom_ribbon(data=sim1_grape,aes(TIME/24, ymin=L,ymax=H), fill = "grey",alpha=0.6)
PK1_grape

sim1_exceed_grape <- subset(Output.DR_grape, Output.DR_grape$DV >= 0.42)
sim1_exceed_grape<- unique(sim1_exceed_grape)
length(sim1_exceed_grape)
sim1_exceed_grape_time <- sim1_grape$TIME[sim1_grape$M >= 0.42]
length(sim1_exceed_grape_time)


# SIMULATION 2 WAS DONE USING THE FULL_BDQ_PK_COVSVEN3new1.csv DATASET (28/3/24)
# ARM: DAY 1-14: 700mg QD, DAY 15-170: 200mg (every 2 days) 
#> 5 Simulations

COSVEN4_B200_L700_M2_grape <-COSVEN4_B200_M2_grape

COSVEN4_B200_L700_M2_grape[COSVEN4_B200_L700_M2_grape$AMT==500,]$AMT<-700
COSVEN4_B200_L700_M2_grape[COSVEN4_B200_L700_M2_grape$AMT==400,]$AMT<-700
write.csv(COSVEN4_B200_L700_M2_grape, "/Users/justmendiola/ibsc/PROJECTS/M2_GRAPE/COSVEN4_B200_L700_M2_grape.csv")


Output.DR1_grape <- read.table.nm("/Users/justmendiola/ibsc/PROJECTS/M2_GRAPE/simulation_2/sdtab2grapefruit2", header = TRUE,stringsAsFactors=FALSE,skip=1)
#Output.DR <- Output.DR[Output.DR$TIME>=1464,]
Output.DR1_grape$DOSE <- "BDQ w/GFJ - 700mg QD 2 week loading dose"
#Output.DR <- Output.DR[Output.DR$TIME>1510,]
Output.DR1_grape$DV <- exp(Output.DR1_grape$DV)
Output.DR1_grape$IPRED <- exp(Output.DR1_grape$IPRED)
Output.DR1_grape <- Output.DR1_grape[Output.DR1_grape$MDV==0,]
# N <- 5
# Output.DR1_grape$ITER <- rep(1:N, each = nrow(Output.DR1_grape)/N)
sim2_grape<-ddply(Output.DR1_grape,.(TIME,DOSE),summarise, M=median(DV), L=quantile(DV,0.05) ,H=quantile(DV, 0.95))



#Output.DR <- Output.DR[Output.DR$TIME>1510,]
PK2_grape <- ggplot()+
  geom_line(data=sim2_grape,aes(TIME/24, M), 
            linewidth=1,colour="#56B4E9")+
  labs(x="Time (Days)", y="M2 concentrations (mg/L)")+
  facet_wrap(~DOSE)+
  geom_hline(aes(yintercept=0.42), col = "red", linetype = "dashed",size=1)+
  theme_bw(base_size = 18)+
  geom_text(data=data.frame(sim2_grape),aes(130,0.47,label="Safety target - 0.42mg/L"),colour = "black",size=6)+
  geom_ribbon(data=sim2_grape,aes(TIME/24, ymin=L,ymax=H), fill = "grey",alpha=0.6)
PK2_grape

sim2_exceed_grape <- Output.DR1_grape$ID[Output.DR1_grape$DV >= 0.42]
sim2_exceed_grape<- unique(sim2_exceed_grape)
length(sim2_exceed_grape)

# SIMULATION 3 WAS DONE USING THE FULL_BDQ_PK_COVSVEN3new1.csv DATASET (16/4/24)
# ARM: DAY 1-14: 400 QD, DAY 15-170: 200mg (every 2 days) 
#> 5 Simulations

COSVEN4_WHO_M2_grape <-COSVEN4_B200_M2_grape

COSVEN4_WHO_M2_grape$AMT[COSVEN4_WHO_M2_grape$AMT >= 200 & 0 <= COSVEN4_WHO_M2_grape$TIME & COSVEN4_WHO_M2_grape$TIME <= 312] <- 400
write.csv(COSVEN4_WHO_M2_grape, "/Users/justmendiola/ibsc/PROJECTS/M2_GRAPE/COSVEN4_WHO_M2_grape.csv")

Output.DR2_grape <- read.table.nm("/Users/justmendiola/ibsc/PROJECTS/M2_GRAPE/simulation_3/sdtab3grapefruit2", header = TRUE,stringsAsFactors=FALSE,skip=1)
#Output.DR <- Output.DR[Output.DR$TIME>=1464,]
Output.DR2_grape$DOSE <- "Scenario 2"
#Output.DR <- Output.DR[Output.DR$TIME>1510,]
Output.DR2_grape$DV <- exp(Output.DR2_grape$DV)
Output.DR2_grape$IPRED <- exp(Output.DR2_grape$IPRED)
Output.DR2_grape <- Output.DR2_grape[Output.DR2_grape$MDV==0,]
# N <- 5
# Output.DR2_grape$ITER <- rep(1:N, each = nrow(Output.DR2_grape)/N)
sim3_grape<-ddply(Output.DR2_grape,.(TIME,DOSE),summarise, M=median(DV), L=quantile(DV,0.05) ,H=quantile(DV, 0.95))

# sim3a_grape<-ddply(Output.DR2_grape,.(TIME,DOSE, ID),summarise, M=median(IPRED), L=quantile(IPRED,0.05) ,H=quantile(IPRED, 0.95))


#Output.DR <- Output.DR[Output.DR$TIME>1510,]
PK3_grape <- ggplot()+
  geom_line(data=sim3_grape,aes(TIME/24, M), 
            linewidth=1,colour="#56B4E9")+
  labs(x="Time (Days)", y="M2 concentrations (mg/L)")+
  facet_wrap(~DOSE)+
  geom_hline(aes(yintercept=0.42), col = "red", linetype = "dashed",size=1)+
  theme_bw(base_size = 16)+
  geom_text(data=data.frame(sim3_grape),aes(130,0.47,label="Safety target - 0.42mg/L"),colour = "black",size=6)+
  geom_ribbon(data=sim3_grape,aes(TIME/24, ymin=L,ymax=H), fill = "grey",alpha=0.6)
PK3_grape



sim3_exceed_grape <- Output.DR2_grape$ID[Output.DR2_grape$DV >= 0.42]
length(sim3_exceed_grape)
sim3_exceed_grape_time <- sim3_grape$TIME[sim3_grape$M >= 0.42]
length(sim3_exceed_grape_time)
sim3_exceed_time <- sim3$TIME[sim3$M >= 0.42]
length(sim3_exceed_time) 

library(gridExtra)


combined_plot <- grid.arrange(PK3, PK3_grape, ncol = 2)

print(combined_plot)

combined_plot_2 <- grid.arrange(PK1_grape, PK3_grape, ncol = 2)
print(combined_plot_2)


# SIMULATION 4 WAS DONE USING THE FULL_BDQ_PK_COVSVEN3new1.csv DATASET (16/4/24)
# ARM: DAY 1-14: 600 QD, DAY 15-170: 200mg (every 2 days) 
#> 5 Simulations


COSVEN4_B200_L600_M2 <-COSVEN4_B200_M2_grape

COSVEN4_B200_L600_M2[COSVEN4_B200_L600_M2$AMT==500,]$AMT<-600
COSVEN4_B200_L600_M2[COSVEN4_B200_L600_M2$AMT==400,]$AMT<-600
write.csv(COSVEN4_B200_L600_M2, "/Users/justmendiola/ibsc/PROJECTS/M2_GRAPE/COSVEN4_B200_L600_M2.csv")


Output.DR3_grape <- read.table.nm("/Users/justmendiola/ibsc/PROJECTS/M2_GRAPE/simulation_4/sdtab4grapefruit", header = TRUE,stringsAsFactors=FALSE,skip=1)
#Output.DR <- Output.DR[Output.DR$TIME>=1464,]
Output.DR3_grape$DOSE <- "BDQ w/GFJ - 600mg QD 2 week loading dose"
#Output.DR <- Output.DR[Output.DR$TIME>1510,]
Output.DR3_grape$DV <- exp(Output.DR3_grape$DV)
Output.DR3_grape$IPRED <- exp(Output.DR3_grape$IPRED)
Output.DR3_grape <- Output.DR3_grape[Output.DR3_grape$MDV==0,]
# N <- 5
# Output.DR3_grape$ITER <- rep(1:N, each = nrow(Output.DR3_grape)/N)
sim4_grape<-ddply(Output.DR3_grape,.(TIME,DOSE),summarise, M=median(DV), L=quantile(DV,0.05) ,H=quantile(DV, 0.95))


#Output.DR <- Output.DR[Output.DR$TIME>1510,]
PK4_grape <- ggplot()+
  geom_line(data=sim4_grape,aes(TIME/24, M), 
            linewidth=1,colour="#56B4E9")+
  labs(x="Time (Days)", y="M2 concentrations (mg/L)")+
  facet_wrap(~DOSE)+
  geom_hline(aes(yintercept=0.42), col = "red", linetype = "dashed",size=1)+
  theme_bw(base_size = 18)+
  geom_text(data=data.frame(sim4_grape),aes(130,0.47,label="Safety target - 0.42mg/L"),colour = "black",size=6)+
  geom_ribbon(data=sim4_grape,aes(TIME/24, ymin=L,ymax=H), fill = "grey",alpha=0.6)
PK4_grape

sim4_exceed_grape <- Output.DR3_grape$ID[Output.DR3_grape$DV >= 0.42]
sim4_exceed_grape<- unique(sim4_exceed_grape)
length(sim4_exceed_grape)

# SIMULATION 5 WAS DONE USING THE FULL_BDQ_PK_COVSVEN3new1.csv DATASET (16/4/24)
# ARM: DAY 1-7: 700 QD, DAY 8-170: 200mg TIW (every 2 days)
#> 5 Simulations


COSVEN4_B200_L7002_M2 <-COSVEN4_B200_M2_grape

COSVEN4_B200_L7002_M2$AMT[COSVEN4_B200_L7002_M2$AMT >= 200 & 0 <= COSVEN4_B200_L7002_M2$TIME & COSVEN4_B200_L7002_M2$TIME <= 168] <- 700
COSVEN4_B200_L7002_M2$AMT[COSVEN4_B200_L7002_M2$AMT >= 200 & 168 < COSVEN4_B200_L7002_M2$TIME] <- 200
rows_to_modify <- COSVEN4_B200_L7002_M2$TIME %in% c(192, 240, 288)
COSVEN4_B200_L7002_M2$AMT[rows_to_modify] <- 0
COSVEN4_B200_L7002_M2$MDV[rows_to_modify] <- 0
COSVEN4_B200_L7002_M2$EVID[rows_to_modify] <- 0
COSVEN4_B200_L7002_M2$CMT[rows_to_modify] <- 5

write.csv(COSVEN4_B200_L7002_M2, "/Users/justmendiola/ibsc/PROJECTS/M2_GRAPE/COSVEN4_B200_L7002_M2.csv")


Output.DR4_grape <- read.table.nm("/Users/justmendiola/ibsc/PROJECTS/M2_GRAPE/simulation_5/sdtab5grapefruit", header = TRUE,stringsAsFactors=FALSE,skip=1)
#Output.DR <- Output.DR[Output.DR$TIME>=1464,]
Output.DR4_grape$DOSE <- "BDQ w/GFJ - 700mg QD 1 week loading dose"
#Output.DR <- Output.DR[Output.DR$TIME>1510,]
Output.DR4_grape$DV <- exp(Output.DR4_grape$DV)
Output.DR4_grape$IPRED <- exp(Output.DR4_grape$IPRED)
Output.DR4_grape <- Output.DR4_grape[Output.DR4_grape$MDV==0,]
N <- 5
# Output.DR4_grape$ITER <- rep(1:N, each = nrow(Output.DR4_grape)/N)
sim5_grape<-ddply(Output.DR4_grape,.(TIME,DOSE),summarise, M=median(DV), L=quantile(DV,0.05) ,H=quantile(DV, 0.95))



#Output.DR <- Output.DR[Output.DR$TIME>1510,]
PK5_grape <- ggplot()+
  geom_line(data=sim5_grape,aes(TIME/24, M), 
            linewidth=1,colour="#56B4E9")+
  labs(x="Time (Days)", y="M2 concentrations (mg/L)")+
  facet_wrap(~DOSE)+
  geom_hline(aes(yintercept=0.42), col = "red", linetype = "dashed",size=1)+
  theme_bw(base_size = 18)+
  geom_text(data=data.frame(sim5_grape),aes(130,0.47,label="Safety target - 0.42mg/L"),colour = "black",size=6)+
  geom_ribbon(data=sim5_grape,aes(TIME/24, ymin=L,ymax=H), fill = "grey",alpha=0.6)
PK5_grape

sim5_exceed_grape <- Output.DR4_grape$ID[Output.DR4_grape$DV >= 0.42]
sim5_exceed_grape<- unique(sim5_exceed_grape)
length(sim5_exceed_grape)

# SIMULATION 6 WAS DONE USING THE FULL_BDQ_PK_COVSVEN3new1.csv DATASET (16/4/24)
# ARM: DAY 1-7: 6 600 QD, DAY 8-170: 200mg TIW (every 2 days)
#> 5 Simulations


COSVEN4_B200_L6002_M2 <-COSVEN4_B200_M2_grape

COSVEN4_B200_L6002_M2$AMT[COSVEN4_B200_L6002_M2$AMT >= 200 & 0 <= COSVEN4_B200_L6002_M2$TIME & COSVEN4_B200_L6002_M2$TIME <= 168] <- 600
COSVEN4_B200_L6002_M2$AMT[COSVEN4_B200_L6002_M2$AMT >= 200 & 168 < COSVEN4_B200_L6002_M2$TIME] <- 200
rows_to_modify_1 <- COSVEN4_B200_L6002_M2$TIME %in% c(192, 240, 288)
COSVEN4_B200_L6002_M2$AMT[rows_to_modify_1] <- 0
COSVEN4_B200_L6002_M2$MDV[rows_to_modify_1] <- 0
COSVEN4_B200_L6002_M2$EVID[rows_to_modify_1] <- 0
COSVEN4_B200_L6002_M2$CMT[rows_to_modify_1] <- 5

write.csv(COSVEN4_B200_L6002_M2, "/Users/justmendiola/ibsc/PROJECTS/M2_GRAPE/simulation_6/COSVEN4_B200_L6002_M2.csv")


Output.DR5_grape <- read.table.nm("/Users/justmendiola/ibsc/PROJECTS/M2_GRAPE/simulation_6/sdtab6grapefruit2", header = TRUE,stringsAsFactors=FALSE,skip=1)
#Output.DR <- Output.DR[Output.DR$TIME>=1464,]
Output.DR5_grape$DOSE <- "BDQ w/GFJ - 600mg QD 1 week loading dose"
#Output.DR <- Output.DR[Output.DR$TIME>1510,]
Output.DR5_grape$DV <- exp(Output.DR5_grape$DV)
Output.DR5_grape$IPRED <- exp(Output.DR5_grape$IPRED)
Output.DR5_grape <- Output.DR5_grape[Output.DR5_grape$MDV==0,]
# N <- 5
# Output.DR5_grape$ITER <- rep(1:N, each = nrow(Output.DR5_grape)/N)
sim6_grape<-ddply(Output.DR5_grape,.(TIME,DOSE),summarise, M=median(DV), L=quantile(DV,0.05) ,H=quantile(DV, 0.95))



#Output.DR <- Output.DR[Output.DR$TIME>1510,]
PK6_grape <- ggplot()+
  geom_line(data=sim6_grape,aes(TIME/24, M), 
            linewidth=1,colour="#56B4E9")+
  labs(x="Time (Days)", y="M2 concentrations (mg/L)")+
  facet_wrap(~DOSE)+
  geom_hline(aes(yintercept=0.42), col = "red", linetype = "dashed",size=1)+
  theme_bw(base_size = 18)+
  geom_text(data=data.frame(sim6_grape),aes(130,0.47,label="Safety target - 0.42mg/L"),colour = "black",size=6)+
  geom_ribbon(data=sim6_grape,aes(TIME/24, ymin=L,ymax=H), fill = "grey",alpha=0.6)
PK6_grape

sim6_exceed_grape <- Output.DR5_grape$ID[Output.DR5_grape$DV >= 0.42]
sim6_exceed_grape<- unique(sim6_exceed_grape)
length(sim6_exceed_grape)


# SIMULATION 7 WAS DONE USING THE FULL_BDQ_PK_COVSVEN3new1.csv DATASET (17/4/24)
# ARM: DAY 1: 700mg QD, DAY 2: 500mg QD, DAY 3-14: 400mg QD, DAY 15-170: 100mg TIW, (every 2 days)
#> 5 Simulations


COSVEN4_B100 <- COSVEN4_B200_M2_grape
COSVEN4_B100$AMT[COSVEN4_B100$AMT == 200] <- 100
write.csv(COSVEN4_B100, "/Users/justmendiola/ibsc/PROJECTS/M2_GRAPE/COSVEN4_B100.csv")


Output.DR6_grape <- read.table.nm("/Users/justmendiola/ibsc/PROJECTS/M2_GRAPE/simulation_7/sdtab7grapefruit", header = TRUE,stringsAsFactors=FALSE,skip=1)
#Output.DR <- Output.DR[Output.DR$TIME>=1464,]
Output.DR6_grape$DOSE <- "BDQ w/GFJ - 100mg/tpw"
#Output.DR <- Output.DR[Output.DR$TIME>1510,]
Output.DR6_grape$DV <- exp(Output.DR6_grape$DV)
Output.DR6_grape$IPRED <- exp(Output.DR6_grape$IPRED)
Output.DR6_grape <- Output.DR6_grape[Output.DR6_grape$MDV==0,]
# N <- 5
# Output.DR6_grape$ITER <- rep(1:N, each = nrow(Output.DR6_grape)/N)
sim7_grape<-ddply(Output.DR6_grape,.(TIME,DOSE),summarise, M=median(DV), L=quantile(DV,0.05) ,H=quantile(DV, 0.95))



#Output.DR <- Output.DR[Output.DR$TIME>1510,]
PK7_grape <- ggplot()+
  geom_line(data=sim7_grape,aes(TIME/24, M), 
            linewidth=1,colour="#56B4E9")+
  labs(x="Time (Days)", y="M2 concentrations (mg/L)")+
  facet_wrap(~DOSE)+
  geom_hline(aes(yintercept=0.42), col = "red", linetype = "dashed",size=1)+
  theme_bw(base_size = 18)+
  geom_text(data=data.frame(sim7_grape),aes(130,0.47,label="Safety target - 0.42mg/L"),colour = "black",size=6)+
  geom_ribbon(data=sim7_grape,aes(TIME/24, ymin=L,ymax=H), fill = "grey",alpha=0.6)
PK7_grape

sim7_exceed_grape <- Output.DR6_grape$ID[Output.DR6_grape$DV >= 0.42]
sim7_exceed_grape<- unique(sim7_exceed_grape)
length(sim7_exceed_grape)

# SIMULATION 8 WAS DONE USING THE FULL_BDQ_PK_COVSVEN3new1.csv DATASET (16/4/24)
# ARM: DAY 1-14: 500 QD, DAY 15-170: 200mg (every 2 days)
#> 5 Simulations


COSVEN4_B200_L500_M2 <-COSVEN4_B200_M2_grape

COSVEN4_B200_L500_M2$AMT[COSVEN4_B200_L500_M2$AMT >= 200 & 0 <= COSVEN4_B200_L500_M2$TIME & COSVEN4_B200_L500_M2$TIME <= 312] <- 500
write.csv(COSVEN4_B200_L500_M2, "/Users/justmendiola/ibsc/PROJECTS/M2_GRAPE/COSVEN4_B200_L500_M2.csv")


Output.DR7_grape <- read.table.nm("/Users/justmendiola/ibsc/PROJECTS/M2_GRAPE/simulation_8/sdtab8grapefruit2", header = TRUE,stringsAsFactors=FALSE,skip=1)
#Output.DR <- Output.DR[Output.DR$TIME>=1464,]
Output.DR7_grape$DOSE <- "BDQ w/GFJ - 500mg QD 2 week loading dose"
#Output.DR <- Output.DR[Output.DR$TIME>1510,]
Output.DR7_grape$DV <- exp(Output.DR7_grape$DV)
Output.DR7_grape$IPRED <- exp(Output.DR7_grape$IPRED)
Output.DR7_grape <- Output.DR7_grape[Output.DR7_grape$MDV==0,]
# N <- 5
# Output.DR7_grape$ITER <- rep(1:N, each = nrow(Output.DR7_grape)/N)
sim8_grape<-ddply(Output.DR7_grape,.(TIME,DOSE),summarise, M=median(DV), L=quantile(DV,0.05) ,H=quantile(DV, 0.95))



#Output.DR <- Output.DR[Output.DR$TIME>1510,]
PK8_grape <- ggplot()+
  geom_line(data=sim8_grape,aes(TIME/24, M), 
            linewidth=1,colour="#56B4E9")+
  labs(x="Time (Days)", y="M2 concentrations (mg/L)")+
  facet_wrap(~DOSE)+
  geom_hline(aes(yintercept=0.42), col = "red", linetype = "dashed",size=1)+
  theme_bw(base_size = 18)+
  geom_text(data=data.frame(sim8_grape),aes(130,0.47,label="Safety target - 0.42mg/L"),colour = "black",size=6)+
  geom_ribbon(data=sim8_grape,aes(TIME/24, ymin=L,ymax=H), fill = "grey",alpha=0.6)
PK8_grape

sim8_exceed_grape <- Output.DR7_grape$ID[Output.DR7_grape$DV >= 0.42]
sim8_exceed_grape<- unique(sim8_exceed_grape)
length(sim8_exceed_grape)

# SIMULATION 9 WAS DONE USING THE FULL_BDQ_PK_COVSVEN3new1.csv DATASET (18/4/24)
# ARM: DAY 1-14: 200 QD, DAY 15-170: 200mg (every 2 days)
#> 5 Simulations


COSVEN4_B200_L200_M2 <-COSVEN4_B200_M2_grape

COSVEN4_B200_L200_M2$AMT[COSVEN4_B200_L200_M2$AMT >= 200 & 0 <= COSVEN4_B200_L200_M2$TIME & COSVEN4_B200_L200_M2$TIME <= 312] <- 200
write.csv(COSVEN4_B200_L200_M2, "/Users/justmendiola/ibsc/PROJECTS/M2_GRAPE/COSVEN4_B200_L200_M2.csv")


Output.DR8_grape <- read.table.nm("/Users/justmendiola/ibsc/PROJECTS/M2_GRAPE/simulation_9/sdtab9grapefruit", header = TRUE,stringsAsFactors=FALSE,skip=1)
#Output.DR <- Output.DR[Output.DR$TIME>=1464,]
Output.DR8_grape$DOSE <- "BDQ w/GFJ - 200mg QD 2 week loading dose"
#Output.DR <- Output.DR[Output.DR$TIME>1510,]
Output.DR8_grape$DV <- exp(Output.DR8_grape$DV)
Output.DR8_grape$IPRED <- exp(Output.DR8_grape$IPRED)
Output.DR8_grape <- Output.DR8_grape[Output.DR8_grape$MDV==0,]
# N <- 5
# Output.DR8_grape$ITER <- rep(1:N, each = nrow(Output.DR8_grape)/N)
sim9_grape<-ddply(Output.DR8_grape,.(TIME,DOSE),summarise, M=median(DV), L=quantile(DV,0.05) ,H=quantile(DV, 0.95))


#Output.DR <- Output.DR[Output.DR$TIME>1510,]
PK9_grape <- ggplot()+
  geom_line(data=sim9_grape,aes(TIME/24, M), 
            linewidth=1,colour="#56B4E9")+
  labs(x="Time (Days)", y="M2 concentrations (mg/L)")+
  facet_wrap(~DOSE)+
  geom_hline(aes(yintercept=0.42), col = "red", linetype = "dashed",size=1)+
  theme_bw(base_size = 18)+
  geom_text(data=data.frame(sim9_grape),aes(130,0.47,label="Safety target - 0.42mg/L"),colour = "black",size=6)+
  geom_ribbon(data=sim9_grape,aes(TIME/24, ymin=L,ymax=H), fill = "grey",alpha=0.6)
PK9_grape

sim9_exceed_grape <- Output.DR8_grape$ID[Output.DR8_grape$DV >= 0.42]
sim9_exceed_grape<- unique(sim9_exceed_grape)
length(sim9_exceed_grape)

# SIMULATION 10 WAS DONE USING THE FULL_BDQ_PK_COVSVEN3new1.csv DATASET (18/4/24)
# ARM: DAY 1-14: 300 QD, DAY 15-170: 200mg (every 2 days)
#> 5 Simulations


COSVEN4_B200_L300_M2 <-COSVEN4_B200_M2_grape

COSVEN4_B200_L300_M2$AMT[COSVEN4_B200_L300_M2$AMT >= 200 & 0 <= COSVEN4_B200_L300_M2$TIME & COSVEN4_B200_L300_M2$TIME <= 312] <- 300
write.csv(COSVEN4_B200_L300_M2, "/Users/justmendiola/ibsc/PROJECTS/M2_GRAPE/COSVEN4_B200_L300_M2.csv")


Output.DR9_grape <- read.table.nm("/Users/justmendiola/ibsc/PROJECTS/M2_GRAPE/simulation_10/sdtab10grapefruit", header = TRUE,stringsAsFactors=FALSE,skip=1)
#Output.DR <- Output.DR[Output.DR$TIME>=1464,]
Output.DR9_grape$DOSE <- "BDQ w/GFJ - 300mg QD 2 week loading dose"
#Output.DR <- Output.DR[Output.DR$TIME>1510,]
Output.DR9_grape$DV <- exp(Output.DR9_grape$DV)
Output.DR9_grape$IPRED <- exp(Output.DR9_grape$IPRED)
Output.DR9_grape <- Output.DR9_grape[Output.DR9_grape$MDV==0,]
# N <- 5
# Output.DR9_grape$ITER <- rep(1:N, each = nrow(Output.DR9_grape)/N)
sim10_grape<-ddply(Output.DR9_grape,.(TIME,DOSE),summarise, M=median(DV), L=quantile(DV,0.05) ,H=quantile(DV, 0.95))



#Output.DR <- Output.DR[Output.DR$TIME>1510,]
PK10_grape <- ggplot()+
  geom_line(data=sim10_grape,aes(TIME/24, M), 
            linewidth=1,colour="#56B4E9")+
  labs(x="Time (Days)", y="M2 concentrations (mg/L)")+
  facet_wrap(~DOSE)+
  geom_hline(aes(yintercept=0.42), col = "red", linetype = "dashed",size=1)+
  theme_bw(base_size = 18)+
  geom_text(data=data.frame(sim10_grape),aes(130,0.47,label="Safety target - 0.42mg/L"),colour = "black",size=6)+
  geom_ribbon(data=sim10_grape,aes(TIME/24, ymin=L,ymax=H), fill = "grey",alpha=0.6)
PK10_grape

sim10_exceed_grape <- Output.DR9_grape$ID[Output.DR9_grape$DV >= 0.42]
sim10_exceed_grape<- unique(sim10_exceed_grape)
length(sim10_exceed_grape)

# SIMULATION 11 WAS DONE USING THE FULL_BDQ_PK_COVSVEN3new1.csv DATASET (18/4/24)
# ARM: DAY 1-7: 6 500 QD, DAY 8-170: 200mg TIW (every 2 days)
#> 5 Simulations


COSVEN4_B200_L5002_M2 <-COSVEN4_B200_M2_grape

COSVEN4_B200_L5002_M2$AMT[COSVEN4_B200_L5002_M2$AMT >= 200 & 0 <= COSVEN4_B200_L5002_M2$TIME & COSVEN4_B200_L5002_M2$TIME <= 168] <- 500
COSVEN4_B200_L5002_M2$AMT[COSVEN4_B200_L5002_M2$AMT >= 200 & 168 < COSVEN4_B200_L5002_M2$TIME] <- 200
rows_to_modify_1 <- COSVEN4_B200_L5002_M2$TIME %in% c(192, 240, 288)
COSVEN4_B200_L5002_M2$AMT[rows_to_modify] <- 0
COSVEN4_B200_L5002_M2$MDV[rows_to_modify] <- 0
COSVEN4_B200_L5002_M2$EVID[rows_to_modify] <- 0
COSVEN4_B200_L5002_M2$CMT[rows_to_modify] <- 5

write.csv(COSVEN4_B200_L5002_M2, "/Users/justmendiola/ibsc/PROJECTS/M2_GRAPE/COSVEN4_B200_L5002_M2.csv")


Output.DR10_grape <- read.table.nm("/Users/justmendiola/ibsc/PROJECTS/M2_GRAPE/simulation_11/sdtab11grapefruit", header = TRUE,stringsAsFactors=FALSE,skip=1)
#Output.DR <- Output.DR[Output.DR$TIME>=1464,]
Output.DR10_grape$DOSE <- "BDQ w/GFJ - 500mg QD 1 week loading dose"
#Output.DR <- Output.DR[Output.DR$TIME>1510,]
Output.DR10_grape$DV <- exp(Output.DR10_grape$DV)
Output.DR10_grape$IPRED <- exp(Output.DR10_grape$IPRED)
Output.DR10_grape <- Output.DR10_grape[Output.DR10_grape$MDV==0,]
# N <- 5
# Output.DR10$ITER <- rep(1:N, each = nrow(Output.DR10)/N)
sim11_grape<-ddply(Output.DR10_grape,.(TIME,DOSE),summarise, M=median(DV), L=quantile(DV,0.05) ,H=quantile(DV, 0.95))



#Output.DR <- Output.DR[Output.DR$TIME>1510,]
PK11_grape <- ggplot()+
  geom_line(data=sim11_grape,aes(TIME/24, M), 
            linewidth=1,colour="#56B4E9")+
  labs(x="Time (Days)", y="M2 concentrations (mg/L)")+
  facet_wrap(~DOSE)+
  geom_hline(aes(yintercept=0.42), col = "red", linetype = "dashed",size=1)+
  theme_bw(base_size = 18)+
  geom_text(data=data.frame(sim11_grape),aes(130,0.47,label="Safety target - 0.42mg/L"),colour = "black",size=6)+
  geom_ribbon(data=sim11_grape,aes(TIME/24, ymin=L,ymax=H), fill = "grey",alpha=0.6)
PK11_grape

sim11_exceed_grape <- Output.DR10_grape$ID[Output.DR10_grape$DV >= 0.42]
sim11_exceed_grape<- unique(sim11_exceed_grape)
length(sim11_exceed_grape)

# SIMULATION 12 WAS DONE USING THE FULL_BDQ_PK_COVSVEN3new1.csv DATASET (15/4/24)
# ARM: DAY 1: 700mg QD, DAY 2: 500mg QD, DAY 3-14: 400mg QD, DAY 15-170: 300mg TIW, (every 2 days)
#> 5 Simulations


COSVEN4_B300 <- COSVEN4_B200_M2_grape
COSVEN4_B300$AMT[COSVEN4_B300$AMT == 200] <- 300
write.csv(COSVEN4_B300, "/Users/justmendiola/ibsc/PROJECTS/M2_GRAPE/COSVEN4_B300.csv")


Output.DR11_grape <- read.table.nm("/Users/justmendiola/ibsc/PROJECTS/M2_GRAPE/simulation_12/sdtab12grapefruit2", header = TRUE,stringsAsFactors=FALSE,skip=1)
#Output.DR <- Output.DR[Output.DR$TIME>=1464,]
Output.DR11_grape$DOSE <- "Scenario 5"
#Output.DR <- Output.DR[Output.DR$TIME>1510,]
Output.DR11_grape$DV <- exp(Output.DR11_grape$DV)
Output.DR11_grape$IPRED <- exp(Output.DR11_grape$IPRED)
Output.DR11_grape <- Output.DR11_grape[Output.DR11_grape$MDV==0,]
# N <- 5
# Output.DR6_grape$ITER <- rep(1:N, each = nrow(Output.DR6_grape)/N)
sim12_grape<-ddply(Output.DR11_grape,.(TIME,DOSE),summarise, M=median(DV), L=quantile(DV,0.05) ,H=quantile(DV, 0.95))



#Output.DR <- Output.DR[Output.DR$TIME>1510,]
PK12_grape <- ggplot()+
  geom_line(data=sim12_grape,aes(TIME/24, M), 
            linewidth=1,colour="#56B4E9")+
  labs(x="Time (Days)", y="M2 concentrations (mg/L)")+
  facet_wrap(~DOSE)+
  geom_hline(aes(yintercept=0.42), col = "red", linetype = "dashed",size=1)+
  theme_bw(base_size = 18)+
  geom_text(data=data.frame(sim12_grape),aes(130,0.47,label="Safety target - 0.42mg/L"),colour = "black",size=6)+
  geom_ribbon(data=sim12_grape,aes(TIME/24, ymin=L,ymax=H), fill = "grey",alpha=0.6)
PK12_grape

sim12_exceed_grape <- Output.DR11_grape$ID[Output.DR11_grape$DV >= 0.42]
sim12_exceed_grape<- unique(sim12_exceed_grape)
length(sim12_exceed_grape)

# SIMULATION 13 WAS DONE USING THE FULL_BDQ_PK_COVSVEN3new1.csv DATASET (15/4/24)
# ARM: DAY 1-14: 400mg QD, DAY 15-170: 300mg TIW, (every 2 days)
#> 5 Simulations

COSVEN4_B300_2 <- COSVEN4_B200_M2_grape
COSVEN4_B300_2$AMT[COSVEN4_B300_2$AMT >= 200 & 0 <= COSVEN4_B300_2$TIME & COSVEN4_B300_2$TIME <= 312] <- 400
COSVEN4_B300_2$AMT[COSVEN4_B300_2$AMT == 200] <- 300
write.csv(COSVEN4_B300_2, "/Users/justmendiola/ibsc/PROJECTS/M2_GRAPE/COSVEN4_B300_2.csv")


Output.DR12_grape <- read.table.nm("/Users/justmendiola/ibsc/PROJECTS/M2_GRAPE/simulation_13/sdtab13grapefruit2", header = TRUE,stringsAsFactors=FALSE,skip=1)
#Output.DR <- Output.DR[Output.DR$TIME>=1464,]
Output.DR12_grape$DOSE <- "Scenario 4"
#Output.DR <- Output.DR[Output.DR$TIME>1510,]
Output.DR12_grape$DV <- exp(Output.DR12_grape$DV)
Output.DR12_grape$IPRED <- exp(Output.DR12_grape$IPRED)
Output.DR12_grape <- Output.DR12_grape[Output.DR12_grape$MDV==0,]
# N <- 5
# Output.DR6_grape$ITER <- rep(1:N, each = nrow(Output.DR6_grape)/N)
sim13_grape<-ddply(Output.DR12_grape,.(TIME,DOSE),summarise, M=median(DV), L=quantile(DV,0.05) ,H=quantile(DV, 0.95))



#Output.DR <- Output.DR[Output.DR$TIME>1510,]
PK13_grape <- ggplot()+
  geom_line(data=sim13_grape,aes(TIME/24, M), 
            linewidth=1,colour="#56B4E9")+
  labs(x="Time (Days)", y="M2 concentrations (mg/L)")+
  facet_wrap(~DOSE)+
  geom_hline(aes(yintercept=0.42), col = "red", linetype = "dashed",size=1)+
  theme_bw(base_size = 18)+
  geom_text(data=data.frame(sim13_grape),aes(130,0.47,label="Safety target - 0.42mg/L"),colour = "black",size=6)+
  geom_ribbon(data=sim13_grape,aes(TIME/24, ymin=L,ymax=H), fill = "grey",alpha=0.6)
PK13_grape

sim13_exceed_grape_time <- sim13_grape$TIME[sim13_grape$M >= 0.42]
length(sim13_exceed_grape_time)


# combining the M2 grape plots using library(gridExtra)
combined_plots_fig <- grid.arrange(PK3_grape, PK1_grape, PK13_grape, PK12_grape, ncol = 2)
