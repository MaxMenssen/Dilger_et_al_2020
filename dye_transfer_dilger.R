#------------------------------------------------------------------------------
#--------------------------- Dye Transfer -------------------------------------
#------------------------------------------------------------------------------

# Loading required packages
library(hnp)
library(tidyverse)
library(lubridate)
library(plyr)
library(emmeans)

# Data import
dat <- read.csv2("https://raw.githubusercontent.com/MaxMenssen/Dilger_et_al_2020/master/dye_transfer_dilger.csv")

# Renaming variables, take Day as date
dat <- dat %>% 
        dplyr::rename(Day=Messtag,
                      Treatment=Behandlung,
                      Run=Repetition) %>% 
        mutate(Day=ymd(Day),
               Treatment=factor(Treatment, levels=c("Control", 
                                                    "NIM-1", 
                                                    "7d_NIM-2",
                                                    "30d_NIM-2",
                                                    "30d_NIM-2_1d_MAT")),
               Run=factor(Run))

head(dat)
str(dat)

levels(dat$Treatment)

#------------------------------------------------------------------------------

# Overview about the data
ggplot(dat, aes(x=Day, y=coupling_rate_injected))+
        theme_bw()+
        geom_point(aes(color=Run))+
        scale_x_date(date_breaks = "15 days", date_labels =  "%d %b")+
        facet_grid(~Treatment)+
        theme(axis.text.x=element_text(angle=90, hjust = 1))

# ATTENTION: Due to the experimental design the treatment effect might be somehow confonded 
# with a time effect.

#------------------------------------------------------------------------------

# Minimum coupling rate that is greater than 0
min_cr <- dat %>% 
        filter(coupling_rate_injected > 0) %>% 
        summarise(min(coupling_rate_injected)) %>% 
        unname() 

# coupling rate + min(coupling rate)
dat$cr_p_min <- dat$coupling_rate_injected + min_cr[,1]

# Fit the model
fit <- lm(log(cr_p_min)~Treatment, dat)

# Mean comparisons (already back transformed from log scale)
comp <- emmeans(fit, specs="Treatment", contr="trt.vs.ctrl", type="response")

# Extracting the means, backtransformation to original scale
ls_means <- comp$emmeans %>% 
        data.frame() %>% 
        group_by(Treatment) %>% 
        transmute(coupling_rate_injected=response-min_cr[,1],
                  lower=lower.CL-min_cr[,1],
                  upper=upper.CL-min_cr[,1])

# Mark significant differences
ls_means$sig_star <- c("", "", "*", "*", "*")      
  
# Save the least square means as csv file
write.csv2(ls_means, "dye_transfer_means_dilger.csv", row.names=FALSE)      


# Save the contrasts 
write.csv2(data.frame(comp$contrasts), "dye_transfer_contrasts_dilger.csv", 
           row.names=FALSE)

#------------------------------------------------------------------------------

# Grafic 2b
ggplot(dat, aes(x=Treatment, y=coupling_rate_injected))+
        theme_bw()+
        geom_point(aes(color=Treatment),
                   alpha=0.6,
                   position=position_jitter(height=0, width=0.1))+
        geom_linerange(data=ls_means, aes(x=Treatment,
                                         ymin=lower,
                                         ymax=upper))+
        geom_point(data=ls_means, 
                   aes(x=Treatment, y=coupling_rate_injected),
                   shape="-", size=6)+
        geom_text(data=ls_means, aes(y=upper+2.6,
                                    x=Treatment,
                                    label=sig_star), size=6)+
        scale_y_continuous(trans="log",
                           breaks=c(0.001, 0.01, 0.1, 1, 10, 25),
                           limits=c(0.0005, 26))+
  scale_x_discrete(breaks=c("Control", 
                            "NIM-1", 
                            "7d_NIM-2",
                            "30d_NIM-2",
                            "30d_NIM-2_1d_MAT"))+
  scale_color_discrete(breaks=c("Control", 
                                "NIM-1", 
                                "7d_NIM-2",
                                "30d_NIM-2",
                                "30d_NIM-2_1d_MAT"),
                       labels=c("Control", 
                                "NIM-1", 
                                "7d NIM-2",
                                "30d NIM-2",
                                "30d NIM-2 1d MAT"))+
  scale_color_manual(values=c("#FF6666",
                              "#CCCC33",
                              "#33CC33",
                              "#0099CC",
                              "#FF99FF"))+
  theme(axis.text.x=element_blank(),
        legend.title.align=0.5,
        legend.position="bottom",
        legend.box.margin=margin(-10, -10, -10, -10),
        legend.title = element_blank(),
        legend.text = element_text(face="bold"),
        strip.text.x = element_text(face="bold.italic"),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(color="black"),
        axis.title.y = element_text(color="black", size=15))+
        ylab("LY dye coupling rate")+
        xlab("")

ggsave("dye_transfer_fig_2.png", width=10, height=9*1.5, units="cm", dpi=600)


