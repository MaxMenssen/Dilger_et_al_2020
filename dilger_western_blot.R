#------------------------------------------------------------------------------
#----------------- Western Blot Nadine Dilger ---------------------------------
#------------------------------------------------------------------------------

# Loading required packages
library(hnp)
library(tidyverse)
library(plyr)
library(emmeans)

# import the data set
dat <- read.csv2("https://raw.githubusercontent.com/MaxMenssen/Dilger_et_al_2020/master/western_blot_dilger.csv", sep=",")

dat$Run <- as.factor(dat$Run)

dat <- dat %>% 
        dplyr::rename(Treatment=Behandlung)

dat$Treatment <- factor(dat$Treatment, levels=c("Control", 
                                                "NIM-1", 
                                                "7d_NIM-2",
                                                "30d_NIM-2",
                                                "30d_NIM-2_1d_MAT"))

# Overview about the data
head(dat)
str(dat)

#------------------------------------------------------------------------------

# Graphical overview about the exp. design
# Show which proteinse were on the same run 
ggplot(dat, aes(x=Treatment, y=Int_House_Cont))+
        theme_bw()+
        geom_point()+
        facet_grid(Run~Protein)+
        scale_x_discrete(breaks=c("Control", "30d_NIM-2",
                                  "30d_NIM-2_1d_MAT", "7d_NIM-2",
                                  "NIM-1"),
                         labels=c("Contr.", "A", "B", "C", "D"))

#------------------------------------------------------------------------------

# Fit a linear model such that log(delta_ct)~Run+Treatment for each protein

# Vector with names of the proteins
proteins <- levels(dat$Protein)

# Open vectors to store results
res_vs_fitted <- vector("list", length=length(proteins))
names(res_vs_fitted) <- proteins

hnp_plot <- vector("list", length=length(proteins))
names(hnp_plot) <- proteins

means <- vector("list", length=length(proteins))
names(means) <- proteins

contrasts <- vector("list", length=length(proteins))
names(contrasts) <- proteins

for(i in 1:length(proteins)){
        
        # Filtering for proteins (without control group)
        dat_fitered <- dat %>% 
                filter(Protein==proteins[i] & Treatment!="Control")
        
        # Fit the linear model 
        fit <- lm(log(Int_House_Cont)~Run+Treatment, dat_fitered)
        
        # Store the residual vs. fitted plot
        res_vs_fitted[[i]] <- plot(fit, which=1)
        
        # Normal plot
        hnp_plot[[i]] <- hnp(fit, halfnormal=FALSE,
                             print.on=TRUE, paint.out=TRUE)
        
        # Mean Comparisons
        comp <- emmeans(fit, 
                        specs="Treatment",  
                        # contr="trt.vs.ctrl", # Kein contr angeben
                        type="response")
        
        # Data frame with ls means
        means_df <- comp %>% 
                data.frame() %>% 
                dplyr::rename(Int_House_Cont=response)
        
        # print(means_df)
        
        # Differ the treatments significantly from 1?
        contr_df <- test(comp,
                         delta=(log(1)),
                         side="both") %>% 
                data.frame()
        
        # print(contr_df)
        # stop()
        
        contrasts[[i]] <- contr_df
        
        # Mark the significant comparisons
        sig_star <- vector(length=length(contr_df$p.value))
        
        for(k in 1:length(contr_df$p.value)){
                
                # print(contr_df$p.value[k])
                
                if(contr_df$p.value[k]<0.05){
                        sig_star[k] <- "*"
                } # else
                        
                        # if(between(contr_df$p.value[k], 0.01, 0.001)){
                        #         sig_star[k] <- "**"
                        # } # else
                        #         
                        #         if(contr_df$p.value[k]<0.001){
                        #                 sig_star[k] <- "***"
                        #         } 
                
                else{
                        sig_star[k] <- ""
                }
                
        }
        
        means_df$sig_star <- c(sig_star)
        
        means[[i]] <- means_df
}

#------------------------------------------------------------------------------
#--------------------- Table with means and CI --------------------------------
#------------------------------------------------------------------------------


means_data_frame <- ldply(means, .fun=data.frame) %>% 
        dplyr::rename(Protein=.id, Treatment=Treatment) %>% 
        data.frame()

write.csv2(means_data_frame, "means_western_blot_dilger.csv", 
           row.names=FALSE)


#------------------------------------------------------------------------------
#--------------- Table with Contrasts and p.values ----------------------------
#------------------------------------------------------------------------------

contrasts_data_frame <- ldply(contrasts, .fun=data.frame) %>% 
        dplyr::rename(Gene=.id) %>% 
        data.frame()

write.csv2(contrasts_data_frame, "contrasts_western_blot_dilger.csv", 
           row.names=FALSE)

#------------------------------------------------------------------------------
#-------------------- Connexin Expression -------------------------------------
#------------------------------------------------------------------------------
levels(dat$Protein)

dat_con <- dat %>% 
        filter(Protein=="Cx26"|
                       Protein=="Cx43"|
                       Protein=="Cx45",
               Treatment!="Control") %>% 
        droplevels()


means_con <- means_data_frame %>% 
        filter(Protein=="Cx26"|
                       Protein=="Cx43"|
                       Protein=="Cx45") %>% 
        droplevels()



ggplot(dat_con, aes(x=Treatment, y=Int_House_Cont))+
        theme_bw()+
        geom_point(aes(color=Treatment),
                   alpha=0.6,
                   position=position_jitter(height=0, width=0.1))+
        geom_hline(yintercept=1, alpha=0.5)+
        geom_linerange(data=means_con, aes(x=Treatment,
                                         ymin=lower.CL,
                                         ymax=upper.CL))+
        geom_point(data=means_con, aes(x=Treatment, y=Int_House_Cont))+
        geom_text(data=means_con, aes(y=upper.CL+0.00001,
                                    x=Treatment,
                                    label=sig_star), size=8)+
        scale_y_continuous(breaks=c(0.025, 0.05, 0.1, 1, 10, 15),
                           limits=c(0.025, 16),
                           trans="log10")+
        scale_x_discrete(breaks=c(# "Control",
                                  "NIM-1",
                                  "7d_NIM-2",
                                  "30d_NIM-2",
                                  "30d_NIM-2_1d_MAT"))+
        scale_color_discrete(breaks=c(#"Control",
                                      "NIM-1",
                                      "7d_NIM-2",
                                      "30d_NIM-2",
                                      "30d_NIM-2_1d_MAT"),
                             labels=c(#"Control",
                                      "NIM-1",
                                      "7d NIM-2",
                                      "30d NIM-2",
                                      "30d NIM-2 1d MAT"))+
        theme(axis.text.x=element_blank(),
              legend.title.align=0.5,
              legend.position="top",
              legend.box.margin=margin(-10, -10, -10, -10),
              legend.title = element_blank(),
              # strip.text.x = element_text(face="italic"),
              axis.ticks.x=element_blank())+
        ylab("Relative protein expression level\n scaled on MSC")+
        xlab("")+
        facet_grid(~Protein)

ggsave("western_blot_connexin_fig_3_b.png", width=13*1.5, height=9*1.5, units="cm")


#------------------------------------------------------------------------------
#----------------------- Neuronal Markers -------------------------------------
#------------------------------------------------------------------------------
levels(dat$Protein)

dat_neu <- dat %>% 
        filter(Protein=="NeuN"|
                       Protein=="Tuj1",
               Treatment!="Control") %>% 
        droplevels()


means_neu <- means_data_frame %>% 
        filter(Protein=="NeuN"|
                       Protein=="Tuj1") %>% 
        droplevels()



ggplot(dat_neu, aes(x=Treatment, y=Int_House_Cont))+
        theme_bw()+
        geom_point(aes(color=Treatment),
                   alpha=0.6,
                   position=position_jitter(height=0, width=0.1))+
        geom_hline(yintercept=1, alpha=0.5)+
        geom_linerange(data=means_neu, aes(x=Treatment,
                                           ymin=lower.CL,
                                           ymax=upper.CL))+
        geom_point(data=means_neu, aes(x=Treatment, y=Int_House_Cont))+
        geom_text(data=means_neu, aes(y=upper.CL+0.00001,
                                      x=Treatment,
                                      label=sig_star), size=8)+
        scale_y_continuous(breaks=c(0.025, 0.05, 0.1, 1, 10, 15),
                           limits=c(0.025, 16),
                           trans="log10")+
        scale_x_discrete(breaks=c(# "Control",
                                  "NIM-1",
                                  "7d_NIM-2",
                                  "30d_NIM-2",
                                  "30d_NIM-2_1d_MAT"))+
        scale_color_discrete(breaks=c(#"Control",
                "NIM-1",
                "7d_NIM-2",
                "30d_NIM-2",
                "30d_NIM-2_1d_MAT"),
                labels=c(#"Control",
                        "NIM-1",
                        "7d NIM-2",
                        "30d NIM-2",
                        "30d NIM-2 1d MAT"))+
        theme(axis.text.x=element_blank(),
              legend.title.align=0.5,
              legend.position="top",
              legend.box.margin=margin(-10, -10, -10, -10),
              legend.title = element_blank(),
              # strip.text.x = element_text(face="italic"),
              axis.ticks.x=element_blank())+
        ylab("Relative protein expression level\n scaled on MSC")+
        xlab("")+
        facet_grid(~Protein)

ggsave("western_blot_neuronal_fig_1_c.png", width=13*1.5, height=9*1.5, units="cm")












