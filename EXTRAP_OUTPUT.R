# BDQ PK and killing curves
EXTRAP_B200 <- read.csv("/Users/justmendiola/ibsc/PROJECTS/PK:PD_simulations/BDQ_MONO_GRAPE/EXTRAP_B200_WHO2CMT4.csv", header = TRUE, skip = 0)

filtered_EXTRAP_B200 <- subset(EXTRAP_B200, EXTRAP_B200$TIME == 0,)
filtered_EXTRAP_700 <- subset(filtered_EXTRAP_B200, AMT == 700,)
unique_IDs <- unique(EXTRAP_B200$ID)

correct_IDs <- c(filtered_EXTRAP_700$ID)
filtered_data <- subset(EXTRAP_B200, ID %in% correct_IDs)

# there are 28 subjects that follow the regimen



# SIMULATION 1 WAS DONE USING THE EXTRAP_B_WHO2.csv DATASET (7/3/24)
# ARM: DAY 1: 700mg QD, DAY 2: 500mg QD, DAY 3-14: 400mg QD, DAY 15-170: 200mg (every 2 days) 
#> 200 Simulations
suppressMessages(setAs("character","dummy.numeric", function(from) as.numeric(from)))
read.table.nm <- function(x,...){
  tmp <- suppressWarnings(read.table(x, fill=T, colClasses="dummy.numeric",...))
  return(tmp[complete.cases(tmp),])}

library(ggplot2)
library(plyr)

B200_EXTRAP_SIMULATION <- read.table.nm("/Users/justmendiola/ibsc/PROJECTS/PK:PD_simulations/BDQ_MONO/simulation_1/sdtabBDQ1", header = T, skip= 1)
B200_EXTRAP_SIMULATION_BUGS <- read.table.nm("/Users/justmendiola/ibsc/PROJECTS/PK:PD_simulations/BDQ_MONO/simulation_1/sdtabBDQ1k", header = T, skip= 1)
# B200_EXTRAP_SIMULATION_HIV <- B200_EXTRAP_SIMULATION_BUGS$ID %in% c(1080, 1081, 6037)

# compartment 2:
B200_EXTRAP_SIMULATION_CMT2 <- subset(B200_EXTRAP_SIMULATION_BUGS, CMT == 2,)
# B200_EXTRAP_SIMULATION_HIV_CMT2 <- subset(B200_EXTRAP_SIMULATION_BUGS, ID == 1080 | 1081 | 6037 & CMT == 2,)

# compartment 4:
B200_EXTRAP_SIMULATION_CMT4 <- subset(B200_EXTRAP_SIMULATION, CMT == 4,)

B200_EXTRAP_SIMULATION_AUC <- ddply(B200_EXTRAP_SIMULATION_CMT4, .(TIME), summarise,
                                    m = median(AUC), p025 = quantile(AUC, 0.025), p975 = quantile(AUC, 0.975))


# AVECONC - average concentration at time intervals

B200_EXTRAP_SIMULATION_AVECONC <- ddply(B200_EXTRAP_SIMULATION, .(TIME), summarise,
                                        m = median(AVECONC), p025 = quantile(AVECONC, 0.025), p975 = quantile(AVECONC, 0.975))

B200_EXTRAP_SIMULATION_AVECONC$DOSE <- "Scenario 3"

B200_AVECONC_plot <- ggplot(B200_EXTRAP_SIMULATION_AVECONC) + 
  geom_line(aes((TIME/24 - 60), m), 
            linewidth=1)+
  labs(x = 'Time (Days)',
       y = 'Avg. Plasma [BDQ] (mg/L)') +
  facet_wrap(~DOSE)+
  geom_hline(yintercept = (0.192), linetype = "dashed", color = "black") +
  geom_hline(yintercept = (3.04), linetype = "dashed", color = "black") +
  theme_bw(base_size = 16)+
  geom_text(data=data.frame(B200_EXTRAP_SIMULATION_AVECONC),aes(100,0.35,label="Fbugs IC50 - 0.192"),colour = "black",size=5) +
  geom_text(data=data.frame(B200_EXTRAP_SIMULATION_AVECONC),aes(100,3.19,label="Sbugs IC50 - 3.04"),colour = "black",size=5) +
  geom_ribbon(data= B200_EXTRAP_SIMULATION_AVECONC, aes(x=(TIME/24 - 60), ymin = p025, ymax = p975), fill = "grey1", alpha = 0.1) +
  xlim(0, 170)
B200_AVECONC_plot



# EBA
B200_EXTRAP_SIMULATION_FBUGS <- ddply(B200_EXTRAP_SIMULATION_CMT2, .(TIME), summarise, # fast growing bacteria
                                      m = median(FBUGS), p025 = quantile(FBUGS, 0.025), p975 = quantile(FBUGS, 0.975))
B200_EXTRAP_SIMULATION_FBUGS$DOSE <- "Scenario 3"

B200_EXTRAP_SIMULATION_SBUGS <- ddply(B200_EXTRAP_SIMULATION_CMT2, .(TIME), summarise, # slow growing bacteria
                                      m = median(SBUGS), p025 = quantile(SBUGS, 0.025), p975 = quantile(SBUGS, 0.975))
B200_EXTRAP_SIMULATION_SBUGS$DOSE <- "Scenario 3"

B200IDS_ZERO_SBUGS <- B200_EXTRAP_SIMULATION_CMT2$ID[B200_EXTRAP_SIMULATION_CMT2$SBUGS == 0 & B200_EXTRAP_SIMULATION_CMT2$TIME > 2640]
B200IDS_ZERO_SBUGS <- unique(B200IDS_ZERO_SBUGS)
print(B200IDS_ZERO_SBUGS)

B200IDS_ZERO_FBUGS <- B200_EXTRAP_SIMULATION_CMT2$ID[B200_EXTRAP_SIMULATION_CMT2$FBUGS == 0 & B200_EXTRAP_SIMULATION_CMT2$TIME > 2640]
B200IDS_ZERO_FBUGS <- unique(B200IDS_ZERO_FBUGS)
print(B200IDS_ZERO_FBUGS)


EBA_line_plot_B200 <- ggplot() + # plot for fast and slow growing populations
  geom_line(data = B200_EXTRAP_SIMULATION_FBUGS, 
            aes((TIME/24 - 60), m, color = "Fast Growing")) +
  geom_line(data = B200_EXTRAP_SIMULATION_SBUGS, 
            aes((TIME/24 - 60), m, color = "Slow Growing")) +
  scale_color_manual(values = c("blue", "red"), 
                     labels = c("Fast Growing", "Slow Growing")) +
  theme_bw(base_size = 18) +
  scale_y_log10() +
  facet_wrap(~DOSE)+
  labs(x = 'Time (Days)',
       y = 'Log10 mtb population (CFU/ml)') +
  xlim(0, 170)

EBA_line_plot_B200

FBUGS_plot_B200 <- ggplot() +
  geom_line(data = B200_EXTRAP_SIMULATION_FBUGS, 
          aes((TIME/24 - 60), m)) +
  geom_ribbon(data= B200_EXTRAP_SIMULATION_FBUGS, aes(x=(TIME/24 - 60), ymin = p025, ymax = p975), fill = "grey1", alpha = 0.1) +
  scale_y_log10() +
  labs(title = 'FBUGS Killing rate of ARM: B200',
       x = 'Time (Days)',
       y = 'Log10 mtb population (CFU/ml)') +
  xlim(0, 170)

FBUGS_plot_B200


SBUGS_plot_B200 <- ggplot() +
  geom_line(data = B200_EXTRAP_SIMULATION_SBUGS, 
            aes((TIME/24 - 60), m)) +
  geom_ribbon(data= B200_EXTRAP_SIMULATION_SBUGS, aes(x=(TIME/24 - 60), ymin = p025, ymax = p975), fill = "grey1", alpha = 0.1) +
  scale_y_log10() +
  labs(title = 'SBUGS Killing rate of ARM: B200',
       x = 'Time (Days)',
       y = 'Log10 mtb population (CFU/ml)') +
  xlim(0, 170)

SBUGS_plot_B200

# SIMULATION 3 WAS DONE USING THE EXTRAP_B.csv DATASET (15/4/24)
# ARM: DAY 1-14: 400 QD, DAY 15-170: 200mg (every 2 days) 
#> 200 Simulations

EXTRAP_WHO <-EXTRAP_B200
EXTRAP_WHO$AMT[EXTRAP_WHO$AMT >= 200 & 1440 <= EXTRAP_WHO$TIME & EXTRAP_WHO$TIME <= 1776] <- 400
write.csv(EXTRAP_WHO, "/Users/justmendiola/ibsc/PROJECTS/PK:PD_simulations/BDQ_MONO/EXTRAP_WHO.csv")

B200_WHO_EXTRAP_SIMULATION <- read.table.nm("/Users/justmendiola/ibsc/PROJECTS/PK:PD_simulations/BDQ_MONO/simulation_3/sdtabBDQ3", header = T, skip= 1)
B200_WHO_EXTRAP_SIMULATION_BUGS <- read.table.nm("/Users/justmendiola/ibsc/PROJECTS/PK:PD_simulations/BDQ_MONO/simulation_3/sdtabBDQ3k", header = T, skip= 1)


# compartment 2:
B200_WHO_EXTRAP_SIMULATION_CMT2 <- subset(B200_WHO_EXTRAP_SIMULATION_BUGS, CMT == 2,)

# compartment 4:
B200_WHO_EXTRAP_SIMULATION_CMT4 <- subset(B200_WHO_EXTRAP_SIMULATION, CMT == 4,)

B200_WHO_EXTRAP_SIMULATION_AUC <- ddply(B200_WHO_EXTRAP_SIMULATION_CMT4, .(TIME), summarise,
                                  m = median(AUC), p025 = quantile(AUC, 0.025), p975 = quantile(AUC, 0.975))


# AVECONC - average concentration at time intervals

B200_WHO_EXTRAP_SIMULATION_AVECONC <- ddply(B200_WHO_EXTRAP_SIMULATION, .(TIME), summarise,
                                             m = median(AVECONC), p025 = quantile(AVECONC, 0.025), p975 = quantile(AVECONC, 0.975))

B200_WHO_EXTRAP_SIMULATION_AVECONC$DOSE <- "Scenario 1/2"

B200_WHO_AVECONC_plot <- ggplot(B200_WHO_EXTRAP_SIMULATION_AVECONC) + 
  geom_line(aes((TIME/24 - 60), m), 
            linewidth=1)+
  labs(x = 'Time (Days)',
       y = 'Avg. Plasma [BDQ] (mg/L)') +
  facet_wrap(~DOSE)+
  geom_hline(yintercept = (0.192), linetype = "dashed", color = "black") +
  geom_hline(yintercept = (3.04), linetype = "dashed", color = "black") +
  theme_bw(base_size = 16)+
  geom_text(data=data.frame(B200_WHO_EXTRAP_SIMULATION_AVECONC),aes(100,0.35,label="Fbugs IC50 - 0.192"),colour = "black",size=5) +
  geom_text(data=data.frame(B200_WHO_EXTRAP_SIMULATION_AVECONC),aes(100,3.19,label="Sbugs IC50 - 3.04"),colour = "black",size=5) +
  geom_ribbon(data= B200_WHO_EXTRAP_SIMULATION_AVECONC, aes(x=(TIME/24 - 60), ymin = p025, ymax = p975), fill = "grey1", alpha = 0.1) +
  xlim(0, 170)
B200_WHO_AVECONC_plot


# EBA
B200_WHO_EXTRAP_SIMULATION_FBUGS <- ddply(B200_WHO_EXTRAP_SIMULATION_CMT2, .(TIME), summarise, # fast growing bacteria
                                           m = median(FBUGS), p025 = quantile(FBUGS, 0.025), p975 = quantile(FBUGS, 0.975))
B200_WHO_EXTRAP_SIMULATION_FBUGS$DOSE <- "Scenario 1/2"

B200_WHO_EXTRAP_SIMULATION_SBUGS <- ddply(B200_WHO_EXTRAP_SIMULATION_CMT2, .(TIME), summarise, # slow growing bacteria
                                           m = median(SBUGS), p025 = quantile(SBUGS, 0.025), p975 = quantile(SBUGS, 0.975))
B200_WHO_EXTRAP_SIMULATION_SBUGS$DOSE <- "Scenario 1/2"

B200_WHO_IDS_ZERO_SBUGS <- B200_WHO_EXTRAP_SIMULATION_CMT2$ID[B200_WHO_EXTRAP_SIMULATION_CMT2$SBUGS == 0 & B200_WHO_EXTRAP_SIMULATION_CMT2$TIME > 2640]
B200_WHO_IDS_ZERO_SBUGS <- unique(B200_WHO_IDS_ZERO_SBUGS)
print(B200_WHO_IDS_ZERO_SBUGS)

B200_WHO_IDS_ZERO_FBUGS <- B200_WHO_EXTRAP_SIMULATION_CMT2$ID[B200_WHO_EXTRAP_SIMULATION_CMT2$FBUGS == 0 & B200_WHO_EXTRAP_SIMULATION_CMT2$TIME > 2640]
B200_WHO_IDS_ZERO_FBUGS <- unique(B200_WHO_IDS_ZERO_FBUGS)
print(B200_WHO_IDS_ZERO_FBUGS)


EBA_line_plot_B200_WHO <- ggplot() + # plot for fast and slow growing populations
  geom_line(data = B200_WHO_EXTRAP_SIMULATION_FBUGS, 
            aes((TIME/24 - 60), m, color = "Fast Growing")) +
  geom_line(data = B200_WHO_EXTRAP_SIMULATION_SBUGS, 
            aes((TIME/24 - 60), m, color = "Slow Growing")) +
  scale_color_manual(values = c("blue", "red"), 
                     labels = c("Fast Growing", "Slow Growing")) +
  theme_bw(base_size = 18) +
  scale_y_log10() +
  facet_wrap(~DOSE)+
  labs(x = 'Time (Days)',
       y = 'Log10 mtb population (CFU/ml)') +
  xlim(0, 170)

EBA_line_plot_B200_WHO

FBUGS_plot_B200_WHO <- ggplot() +
  geom_line(data = B200_WHO_EXTRAP_SIMULATION_FBUGS, 
            aes((TIME/24 - 60), m)) +
  geom_ribbon(data= B200_WHO_EXTRAP_SIMULATION_FBUGS, aes(x=(TIME/24 - 60), ymin = p025, ymax = p975), fill = "grey1", alpha = 0.1) +
  scale_y_log10() +
  labs(title = 'FBUGS Killing rate of ARM: B200_WHO',
       x = 'Time (Days)',
       y = 'Log10 mtb population (CFU/ml)') +
  xlim(0, 170)

FBUGS_plot_B200_WHO


SBUGS_plot_B200_WHO <- ggplot() +
  geom_line(data = B200_WHO_EXTRAP_SIMULATION_SBUGS, 
            aes((TIME/24 - 60), m)) +
  geom_ribbon(data= B200_WHO_EXTRAP_SIMULATION_SBUGS, aes(x=(TIME/24 - 60), ymin = p025, ymax = p975), fill = "grey1", alpha = 0.1) +
  scale_y_log10() +
  labs(title = 'SBUGS Killing rate of ARM: B200_WHO',
       x = 'Time (Days)',
       y = 'Log10 mtb population (CFU/ml)') +
  xlim(0, 170)
SBUGS_plot_B200_WHO

# SIMULATION 6 WAS DONE USING THE FULL_BDQ_PK_COVSVEN3new1.csv DATASET (16/3/24)
# ARM: DAY 1-7: 6 600 QD, DAY 8-170: 200mg TIW (every 2 days)
#> 5 Simulations


EXTRAP_B200_L6002 <- EXTRAP_B200

EXTRAP_B200_L6002$AMT[EXTRAP_B200_L6002$AMT >= 200 & 1440 <= EXTRAP_B200_L6002$TIME & EXTRAP_B200_L6002$TIME <= 1608] <- 600
EXTRAP_B200_L6002$AMT[EXTRAP_B200_L6002$AMT >= 200 & 1608 < EXTRAP_B200_L6002$TIME] <- 200
rows_to_modify_1 <- EXTRAP_B200_L6002$TIME %in% c(1656, 1704, 1752)
EXTRAP_B200_L6002$AMT[rows_to_modify_1] <- 0
EXTRAP_B200_L6002$MDV[rows_to_modify_1] <- 0
EXTRAP_B200_L6002$EVID[rows_to_modify_1] <- 0
EXTRAP_B200_L6002$RATE[rows_to_modify_1] <- 0
EXTRAP_B200_L6002 <- EXTRAP_B200_L6002[!(EXTRAP_B200_L6002$CMT == 3 & EXTRAP_B200_L6002$TIME %in% c(1656, 1704, 1752)), ]
write.csv(EXTRAP_B200_L6002, "/Users/justmendiola/ibsc/PROJECTS/PK:PD_simulations/BDQ_MONO/EXTRAP_B200_L6002.csv")

B200_L6002_EXTRAP_SIMULATION <- read.table.nm("/Users/justmendiola/ibsc/PROJECTS/PK:PD_simulations/BDQ_MONO/simulation_6/sdtabBDQ6", header = T, skip= 1)
B200_L6002_EXTRAP_SIMULATION_BUGS <- read.table.nm("/Users/justmendiola/ibsc/PROJECTS/PK:PD_simulations/BDQ_MONO/simulation_6/sdtabBDQ6k", header = T, skip= 1)


# compartment 2:
B200_L6002_EXTRAP_SIMULATION_CMT2 <- subset(B200_L6002_EXTRAP_SIMULATION_BUGS, CMT == 2,)

# compartment 4:
B200_L6002_EXTRAP_SIMULATION_CMT4 <- subset(B200_L6002_EXTRAP_SIMULATION, CMT == 4,)

B200_L6002_EXTRAP_SIMULATION_AUC <- ddply(B200_L6002_EXTRAP_SIMULATION_CMT4, .(TIME), summarise,
                                        m = median(AUC), p025 = quantile(AUC, 0.025), p975 = quantile(AUC, 0.975))


# AVECONC - average concentration at time intervals

B200_L6002_EXTRAP_SIMULATION_AVECONC <- ddply(B200_L6002_EXTRAP_SIMULATION, .(TIME), summarise,
                                            m = median(AVECONC), p025 = quantile(AVECONC, 0.025), p975 = quantile(AVECONC, 0.975))

B200_L6002_EXTRAP_SIMULATION_AVECONC$DOSE <- "BDQ w/GFJ - Scenario 6"

B200_L6002_AVECONC_plot <- ggplot(B200_L6002_EXTRAP_SIMULATION_AVECONC) + 
  geom_line(aes((TIME/24 - 60), m), 
            linewidth=1)+
  labs(x = 'Time (Days)',
       y = 'Avg. Plasma [BDQ] (mg/L)') +
  facet_wrap(~DOSE)+
  geom_hline(yintercept = (0.192), linetype = "dashed", color = "black") +
  geom_hline(yintercept = (3.04), linetype = "dashed", color = "grey") +
  theme_bw(base_size = 18)+
  geom_text(data=data.frame(B200_L6002_EXTRAP_SIMULATION_AVECONC),aes(130,0.45,label="Fbugs IC50 - 0.192"),colour = "black",size=6) +
  geom_text(data=data.frame(B200_L6002_EXTRAP_SIMULATION_AVECONC),aes(130,3.19,label="Sbugs IC50 - 3.04"),colour = "black",size=6) +
  geom_ribbon(data= B200_L6002_EXTRAP_SIMULATION_AVECONC, aes(x=(TIME/24 - 60), ymin = p025, ymax = p975), fill = "grey1", alpha = 0.1) +
  xlim(0, 170)
B200_L6002_AVECONC_plot


# EBA
B200_L6002_EXTRAP_SIMULATION_FBUGS <- ddply(B200_L6002_EXTRAP_SIMULATION_CMT2, .(TIME), summarise, # fast growing bacteria
                                          m = median(FBUGS), p025 = quantile(FBUGS, 0.025), p975 = quantile(FBUGS, 0.975))
B200_L6002_EXTRAP_SIMULATION_FBUGS$DOSE <- "BDQ - Scenario 6"

B200_L6002_EXTRAP_SIMULATION_SBUGS <- ddply(B200_L6002_EXTRAP_SIMULATION_CMT2, .(TIME), summarise, # slow growing bacteria
                                          m = median(SBUGS), p025 = quantile(SBUGS, 0.025), p975 = quantile(SBUGS, 0.975))
B200_L6002_EXTRAP_SIMULATION_SBUGS$DOSE <- "BDQ - Scenario 6"

B200_L6002_IDS_ZERO_SBUGS <- B200_L6002_EXTRAP_SIMULATION_CMT2$ID[B200_L6002_EXTRAP_SIMULATION_CMT2$SBUGS == 0 & B200_L6002_EXTRAP_SIMULATION_CMT2$TIME > 2640]
B200_L6002_IDS_ZERO_SBUGS <- unique(B200_L6002_IDS_ZERO_SBUGS)
print(B200_L6002_IDS_ZERO_SBUGS)

B200_L6002_IDS_ZERO_FBUGS <- B200_L6002_EXTRAP_SIMULATION_CMT2$ID[B200_L6002_EXTRAP_SIMULATION_CMT2$FBUGS == 0 & B200_L6002_EXTRAP_SIMULATION_CMT2$TIME > 2640]
B200_L6002_IDS_ZERO_FBUGS <- unique(B200_L6002_IDS_ZERO_FBUGS)
print(B200_L6002_IDS_ZERO_FBUGS)


EBA_line_plot_B200_L6002 <- ggplot() + # plot for fast and slow growing populations
  geom_line(data = B200_L6002_EXTRAP_SIMULATION_FBUGS, 
            aes((TIME/24 - 60), m, color = "Fast Growing")) +
  geom_line(data = B200_L6002_EXTRAP_SIMULATION_SBUGS, 
            aes((TIME/24 - 60), m, color = "Slow Growing")) +
  scale_color_manual(values = c("blue", "red"), 
                     labels = c("Fast Growing", "Slow Growing")) +
  theme_bw(base_size = 18) +
  scale_y_log10() +
  facet_wrap(~DOSE)+
  labs(x = 'Time (Days)',
       y = 'Log10 mtb population (CFU/ml)') +
  xlim(0, 170)

EBA_line_plot_B200_L6002

FBUGS_plot_B200_L6002 <- ggplot() +
  geom_line(data = B200_L6002_EXTRAP_SIMULATION_FBUGS, 
            aes((TIME/24 - 60), m)) +
  geom_ribbon(data= B200_L6002_EXTRAP_SIMULATION_FBUGS, aes(x=(TIME/24 - 60), ymin = p025, ymax = p975), fill = "grey1", alpha = 0.1) +
  scale_y_log10() +
  labs(title = 'FBUGS Killing rate of ARM: B200_L6002',
       x = 'Time (Days)',
       y = 'Log10 mtb population (CFU/ml)') +
  xlim(0, 170)

FBUGS_plot_B200_L6002


SBUGS_plot_B200_L6002 <- ggplot() +
  geom_line(data = B200_L6002_EXTRAP_SIMULATION_SBUGS, 
            aes((TIME/24 - 60), m)) +
  geom_ribbon(data= B200_L6002_EXTRAP_SIMULATION_SBUGS, aes(x=(TIME/24 - 60), ymin = p025, ymax = p975), fill = "grey1", alpha = 0.1) +
  scale_y_log10() +
  labs(title = 'SBUGS Killing rate of ARM: B200_L6002',
       x = 'Time (Days)',
       y = 'Log10 mtb population (CFU/ml)') +
  xlim(0, 170)
SBUGS_plot_B200_L6002

# SIMULATION 12 WAS DONE USING THE EXTRAP_B.csv DATASET (15/4/24)
# ARM: DAY 1: 700mg, Day 2: 500mg, Day 3-14: 400mg QD DAY 15-170: 300mg (every 2 days) 
#> 200 Simulations

EXTRAP_B300 <- EXTRAP_B200
EXTRAP_B300$AMT[EXTRAP_B300$AMT == 200] <- 300
write.csv(EXTRAP_B300, "/Users/justmendiola/ibsc/PROJECTS/PK:PD_simulations/BDQ_MONO/EXTRAP_B300.csv")

B300_EXTRAP_SIMULATION <- read.table.nm("/Users/justmendiola/ibsc/PROJECTS/PK:PD_simulations/BDQ_MONO/simulation_12/sdtabBDQ12", header = T, skip= 1)
B300_EXTRAP_SIMULATION_BUGS <- read.table.nm("/Users/justmendiola/ibsc/PROJECTS/PK:PD_simulations/BDQ_MONO/simulation_12/sdtabBDQ12k", header = T, skip= 1)


# compartment 2:
B300_EXTRAP_SIMULATION_CMT2 <- subset(B300_EXTRAP_SIMULATION_BUGS, CMT == 2,)

# compartment 4:
B300_EXTRAP_SIMULATION_CMT4 <- subset(B300_EXTRAP_SIMULATION, CMT == 4,)

B300_EXTRAP_SIMULATION_AUC <- ddply(B300_EXTRAP_SIMULATION_CMT4, .(TIME), summarise,
                                          m = median(AUC), p025 = quantile(AUC, 0.025), p975 = quantile(AUC, 0.975))


# AVECONC - average concentration at time intervals

B300_EXTRAP_SIMULATION_AVECONC <- ddply(B300_EXTRAP_SIMULATION, .(TIME), summarise,
                                              m = median(AVECONC), p025 = quantile(AVECONC, 0.025), p975 = quantile(AVECONC, 0.975))

B300_EXTRAP_SIMULATION_AVECONC$DOSE <- "Scenario 5"

B300_AVECONC_plot <- ggplot(B300_EXTRAP_SIMULATION_AVECONC) + 
  geom_line(aes((TIME/24 - 60), m), 
            linewidth=1)+
  labs(x = 'Time (Days)',
       y = 'Avg. Plasma [BDQ] (mg/L)') +
  facet_wrap(~DOSE)+
  geom_hline(yintercept = (0.192), linetype = "dashed", color = "black") +
  geom_hline(yintercept = (3.04), linetype = "dashed", color = "black") +
  theme_bw(base_size = 16)+
  geom_text(data=data.frame(B300_EXTRAP_SIMULATION_AVECONC),aes(100,0.35,label="Fbugs IC50 - 0.192"),colour = "black",size=5) +
  geom_text(data=data.frame(B300_EXTRAP_SIMULATION_AVECONC),aes(100,3.19,label="Sbugs IC50 - 3.04"),colour = "black",size=5) +
  geom_ribbon(data= B300_EXTRAP_SIMULATION_AVECONC, aes(x=(TIME/24 - 60), ymin = p025, ymax = p975), fill = "grey1", alpha = 0.1) +
  xlim(0, 170)
B300_AVECONC_plot



# EBA
B300_EXTRAP_SIMULATION_FBUGS <- ddply(B300_EXTRAP_SIMULATION_CMT2, .(TIME), summarise, # fast growing bacteria
                                            m = median(FBUGS), p025 = quantile(FBUGS, 0.025), p975 = quantile(FBUGS, 0.975))
B300_EXTRAP_SIMULATION_FBUGS$DOSE <- "Scenario 5"

B300_EXTRAP_SIMULATION_SBUGS <- ddply(B300_EXTRAP_SIMULATION_CMT2, .(TIME), summarise, # slow growing bacteria
                                            m = median(SBUGS), p025 = quantile(SBUGS, 0.025), p975 = quantile(SBUGS, 0.975))
B300_EXTRAP_SIMULATION_SBUGS$DOSE <- "Scenario 5"

B300_IDS_ZERO_SBUGS <- B300_EXTRAP_SIMULATION_CMT2$ID[B300_EXTRAP_SIMULATION_CMT2$SBUGS == 0 & B300_EXTRAP_SIMULATION_CMT2$TIME > 2640]
B300_IDS_ZERO_SBUGS <- unique(B300_IDS_ZERO_SBUGS)
print(B300_IDS_ZERO_SBUGS)

B300_IDS_ZERO_FBUGS <- B300_EXTRAP_SIMULATION_CMT2$ID[B300_EXTRAP_SIMULATION_CMT2$FBUGS == 0 & B300_EXTRAP_SIMULATION_CMT2$TIME > 2640]
B300_IDS_ZERO_FBUGS <- unique(B300_IDS_ZERO_FBUGS)
print(B300_IDS_ZERO_FBUGS)


EBA_line_plot_B300 <- ggplot() + # plot for fast and slow growing populations
  geom_line(data = B300_EXTRAP_SIMULATION_FBUGS, 
            aes((TIME/24 - 60), m, color = "Fast Growing")) +
  geom_line(data = B300_EXTRAP_SIMULATION_SBUGS, 
            aes((TIME/24 - 60), m, color = "Slow Growing")) +
  scale_color_manual(values = c("blue", "red"), 
                     labels = c("Fast Growing", "Slow Growing")) +
  theme_bw(base_size = 18) +
  scale_y_log10() +
  facet_wrap(~DOSE)+
  labs(x = 'Time (Days)',
       y = 'Log10 mtb population (CFU/ml)') +
  xlim(0, 170)

EBA_line_plot_B300

FBUGS_plot_B300 <- ggplot() +
  geom_line(data = B300_EXTRAP_SIMULATION_FBUGS, 
            aes((TIME/24 - 60), m)) +
  geom_ribbon(data= B300_EXTRAP_SIMULATION_FBUGS, aes(x=(TIME/24 - 60), ymin = p025, ymax = p975), fill = "grey1", alpha = 0.1) +
  scale_y_log10() +
  labs(title = 'FBUGS Killing rate of ARM: B300',
       x = 'Time (Days)',
       y = 'Log10 mtb population (CFU/ml)') +
  xlim(0, 170)

FBUGS_plot_B300


SBUGS_plot_B300 <- ggplot() +
  geom_line(data = B300_EXTRAP_SIMULATION_SBUGS, 
            aes((TIME/24 - 60), m)) +
  geom_ribbon(data= B300_EXTRAP_SIMULATION_SBUGS, aes(x=(TIME/24 - 60), ymin = p025, ymax = p975), fill = "grey1", alpha = 0.1) +
  scale_y_log10() +
  labs(title = 'SBUGS Killing rate of ARM: B300',
       x = 'Time (Days)',
       y = 'Log10 mtb population (CFU/ml)') +
  xlim(0, 170)
SBUGS_plot_B300



# SIMULATION 13 WAS DONE USING THE EXTRAP_B.csv DATASET (15/4/24)
# ARM: DAY 1-14: 400 QD, DAY 15-170: 200mg (every 2 days) 
#> 200 Simulations
EXTRAP_B300_2 <- EXTRAP_B200
EXTRAP_B300_2$AMT[EXTRAP_B300_2$AMT >= 200 & 1440 <= EXTRAP_B300_2$TIME & EXTRAP_B300_2$TIME <= 1776] <- 400
EXTRAP_B300_2$AMT[EXTRAP_B300_2$AMT == 200] <- 300
write.csv(EXTRAP_B300_2, "/Users/justmendiola/ibsc/PROJECTS/PK:PD_simulations/BDQ_MONO//EXTRAP_B300_2.csv")

B300_2_EXTRAP_SIMULATION <- read.table.nm("/Users/justmendiola/ibsc/PROJECTS/PK:PD_simulations/BDQ_MONO/simulation_13/sdtabBDQ13", header = T, skip= 1)
B300_2_EXTRAP_SIMULATION_BUGS <- read.table.nm("/Users/justmendiola/ibsc/PROJECTS/PK:PD_simulations/BDQ_MONO/simulation_13/sdtabBDQ13k", header = T, skip= 1)


# compartment 2:
B300_2_EXTRAP_SIMULATION_CMT2 <- subset(B300_2_EXTRAP_SIMULATION_BUGS, CMT == 2,)

# compartment 4:
B300_2_EXTRAP_SIMULATION_CMT4 <- subset(B300_2_EXTRAP_SIMULATION, CMT == 4,)

B300_2_EXTRAP_SIMULATION_AUC <- ddply(B300_2_EXTRAP_SIMULATION_CMT4, .(TIME), summarise,
                                    m = median(AUC), p025 = quantile(AUC, 0.025), p975 = quantile(AUC, 0.975))


# AVECONC - average concentration at time intervals

B300_2_EXTRAP_SIMULATION_AVECONC <- ddply(B300_2_EXTRAP_SIMULATION, .(TIME), summarise,
                                        m = median(AVECONC), p025 = quantile(AVECONC, 0.025), p975 = quantile(AVECONC, 0.975))

B300_2_EXTRAP_SIMULATION_AVECONC$DOSE <- "Scenario 4"

B300_2_AVECONC_plot <- ggplot(B300_2_EXTRAP_SIMULATION_AVECONC) + 
  geom_line(aes((TIME/24 - 60), m), 
            linewidth=1)+
  labs(x = 'Time (Days)',
       y = 'Avg. Plasma [BDQ] (mg/L)') +
  facet_wrap(~DOSE)+
  geom_hline(yintercept = (0.192), linetype = "dashed", color = "black") +
  geom_hline(yintercept = (3.04), linetype = "dashed", color = "black") +
  theme_bw(base_size = 16)+
  geom_text(data=data.frame(B300_2_EXTRAP_SIMULATION_AVECONC),aes(100,0.35,label="Fbugs IC50 - 0.192"),colour = "black",size=5) +
  geom_text(data=data.frame(B300_2_EXTRAP_SIMULATION_AVECONC),aes(100,3.19,label="Sbugs IC50 - 3.04"),colour = "black",size=5) +
  geom_ribbon(data= B300_2_EXTRAP_SIMULATION_AVECONC, aes(x=(TIME/24 - 60), ymin = p025, ymax = p975), fill = "grey1", alpha = 0.1) +
  xlim(0, 170)
B300_2_AVECONC_plot



# EBA
B300_2_EXTRAP_SIMULATION_FBUGS <- ddply(B300_2_EXTRAP_SIMULATION_CMT2, .(TIME), summarise, # fast growing bacteria
                                      m = median(FBUGS), p025 = quantile(FBUGS, 0.025), p975 = quantile(FBUGS, 0.975))
B300_2_EXTRAP_SIMULATION_FBUGS$DOSE <- "Scenario 4"

B300_2_EXTRAP_SIMULATION_SBUGS <- ddply(B300_2_EXTRAP_SIMULATION_CMT2, .(TIME), summarise, # slow growing bacteria
                                      m = median(SBUGS), p025 = quantile(SBUGS, 0.025), p975 = quantile(SBUGS, 0.975))
B300_2_EXTRAP_SIMULATION_SBUGS$DOSE <- "Scenario 4"

B300_2_IDS_ZERO_SBUGS <- B300_2_EXTRAP_SIMULATION_CMT2$ID[B300_2_EXTRAP_SIMULATION_CMT2$SBUGS == 0 & B300_2_EXTRAP_SIMULATION_CMT2$TIME > 2640]
B300_2_IDS_ZERO_SBUGS <- unique(B300_2_IDS_ZERO_SBUGS)
print(B300_2_IDS_ZERO_SBUGS)

B300_2_IDS_ZERO_FBUGS <- B300_2_EXTRAP_SIMULATION_CMT2$ID[B300_2_EXTRAP_SIMULATION_CMT2$FBUGS == 0 & B300_2_EXTRAP_SIMULATION_CMT2$TIME > 2640]
B300_2_IDS_ZERO_FBUGS <- unique(B300_2_IDS_ZERO_FBUGS)
print(B300_2_IDS_ZERO_FBUGS)


EBA_line_plot_B300_2 <- ggplot() + # plot for fast and slow growing populations
  geom_line(data = B300_2_EXTRAP_SIMULATION_FBUGS, 
            aes((TIME/24 - 60), m, color = "Fast Growing")) +
  geom_line(data = B300_2_EXTRAP_SIMULATION_SBUGS, 
            aes((TIME/24 - 60), m, color = "Slow Growing")) +
  scale_color_manual(values = c("blue", "red"), 
                     labels = c("Fast Growing", "Slow Growing")) +
  theme_bw(base_size = 16) +
  scale_y_log10() +
  facet_wrap(~DOSE)+
  labs(x = 'Time (Days)',
       y = 'Log10 mtb population (CFU/ml)') +
  xlim(0, 170)

EBA_line_plot_B300_2

FBUGS_plot_B300_2 <- ggplot() +
  geom_line(data = B300_2_EXTRAP_SIMULATION_FBUGS, 
            aes((TIME/24 - 60), m)) +
  geom_ribbon(data= B300_2_EXTRAP_SIMULATION_FBUGS, aes(x=(TIME/24 - 60), ymin = p025, ymax = p975), fill = "grey1", alpha = 0.1) +
  scale_y_log10() +
  labs(title = 'FBUGS Killing rate of ARM: B300_2',
       x = 'Time (Days)',
       y = 'Log10 mtb population (CFU/ml)') +
  xlim(0, 170)

FBUGS_plot_B300_2


SBUGS_plot_B300_2 <- ggplot() +
  geom_line(data = B300_2_EXTRAP_SIMULATION_SBUGS, 
            aes((TIME/24 - 60), m)) +
  geom_ribbon(data= B300_2_EXTRAP_SIMULATION_SBUGS, aes(x=(TIME/24 - 60), ymin = p025, ymax = p975), fill = "grey1", alpha = 0.1) +
  scale_y_log10() +
  labs(title = 'SBUGS Killing rate of ARM: B300_2',
       x = 'Time (Days)',
       y = 'Log10 mtb population (CFU/ml)') +
  xlim(0, 170)
SBUGS_plot_B300_2



# Combine FBUGS and SBUGS datasets for TSS plots
combined_data <- rbind(
  cbind(B300_2_EXTRAP_SIMULATION_FBUGS, Phenotype = "Fast-growing Mtb"),
  cbind(B300_2_EXTRAP_SIMULATION_SBUGS, Phenotype = "Slow-growing Mtb"),
  cbind(B300_EXTRAP_SIMULATION_FBUGS, Phenotype = "Fast-growing Mtb"),
  cbind(B300_EXTRAP_SIMULATION_SBUGS, Phenotype = "Slow-growing Mtb"),
  cbind(B200_WHO_EXTRAP_SIMULATION_FBUGS, Phenotype = "Fast-growing Mtb"),
  cbind(B200_WHO_EXTRAP_SIMULATION_SBUGS, Phenotype = "Slow-growing Mtb"),
  cbind(B200_EXTRAP_SIMULATION_FBUGS, Phenotype = "Fast-growing Mtb"),
  cbind(B200_EXTRAP_SIMULATION_SBUGS, Phenotype = "Slow-growing Mtb")
)

# Plot
EBA_line_plot <- ggplot(combined_data, aes(x = (TIME / 24 - 60), y = m, color = Phenotype)) +
  geom_line() +
  scale_color_manual(values = c("blue", "red")) +
  facet_wrap(~ DOSE) +
  theme_bw(base_size = 18) +
  scale_y_log10() +
  labs(x = 'Time (Days)',
       y = 'Log10 mtb population (CFU/ml)',
       title = 'Time-to-Sterilisation Graphs') +
  xlim(0, 170)

EBA_line_plot


# Combine all dataframes for combined PK graph plot
combined_data <- rbind(
  transform(B300_2_EXTRAP_SIMULATION_AVECONC, Scenario = "Scenario 4"),
  transform(B300_EXTRAP_SIMULATION_AVECONC, Scenario = "Scenario 5"),
  transform(B200_WHO_EXTRAP_SIMULATION_AVECONC, Scenario = "Scenario 1/2"),
  transform(B200_EXTRAP_SIMULATION_AVECONC, Scenario = "Scenario 3")
)

# Plot
combined_AVECONC_plot <- ggplot(combined_data, aes(x = (TIME / 24 - 60), y = m)) + 
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = p025, ymax = p975), fill = "grey20", alpha = 0.1) +
  geom_hline(yintercept = 0.192, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 3.04, linetype = "dashed", color = "black") +
  geom_text(data = data.frame(combined_data), aes(x = 100, y = 0.38, label = "Fbugs IC50 - 0.192"), color = "black", size = 5) +
  geom_text(data = data.frame(combined_data), aes(x = 100, y = 3.24, label = "Sbugs IC50 - 3.04"), color = "black", size = 5) +
  labs(x = 'Time (Days)',
       y = 'Avg. Plasma [BDQ] (mg/L)',
       title = 'Concentration-Time Graphs') +
  xlim(0, 170) +
  facet_wrap(~ DOSE) +
  theme_bw(base_size = 16)

combined_AVECONC_plot
