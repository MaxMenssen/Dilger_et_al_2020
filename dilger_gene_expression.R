#------------------------------------------------------------------------------
#---------------------- qPCR Nadine Dilger ------------------------------------
#------------------------------------------------------------------------------

# Loading required packages
library(hnp)
library(tidyverse)
library(plyr)
library(emmeans)

# import the data set directly from github
dat <- read.csv2("https://raw.githubusercontent.com/MaxMenssen/Dilger_et_al_2020/master/gene_expression_dilger.csv")
dat$Run <- as.factor(dat$Run)

# Overview about the data
head(dat)
str(dat)

# Rename the treatment variable
dat <- dat %>% 
        dplyr::rename(Treatment=BehaNAlung)

# Reorder the Treatment levels
dat$Treatment <- factor(dat$Treatment, levels=c("Control", 
                                            "NIM-1", 
                                            "7d_NIM-2",
                                            "30d_NIM-2",
                                            "30d_NIM-2_1d_MAT"))
levels(dat$Treatment)


#------------------------------------------------------------------------------

# Graphical overview about the exp. design
# Show which Genes were on the same plate 
ggplot(dat, aes(x=Treatment, y=delta_ct))+
        theme_bw()+
        geom_point()+
        facet_grid(Run~Gen)+
        scale_x_discrete(breaks=c("Control", "30d_NIM-2",
                                  "30d_NIM-2_1d_MAT", "7d_NIM-2", 
                                  "NIM-1"),
                         labels=c("Contr.", "A", "B", "C", "D"))

#------------------------------------------------------------------------------

# Fit a linear model such that log(delta_ct)~Run+Treatment for each gene

# Vector with names of the genes
genes <- levels(dat$Gen)

# Open vectors to store results
res_vs_fitted <- vector("list", length=length(genes))
names(res_vs_fitted) <- genes

hnp_plot <- vector("list", length=length(genes))
names(hnp_plot) <- genes

means <- vector("list", length=length(genes))
names(means) <- genes

contrasts <- vector("list", length=length(genes))
names(contrasts) <- genes

for(i in 1:length(genes)){
        
        # Filtern für TAGLIN
        dat_fitered <- dat %>% 
                filter(Gen==genes[i])
        
        # Fit the linear model 
        fit <- lm(log(delta_ct)~Run+Treatment, dat_fitered)
        
        # Store the residual vs. fitted plot
        res_vs_fitted[[i]] <- plot(fit, which=1)
        
        # Normal plot
        hnp_plot[[i]] <- hnp(fit, halfnormal=FALSE,
                             print.on=TRUE, paint.out=TRUE)
        
        # Mean Comparisons
        comp <- emmeans(fit, 
                        specs="Treatment",  
                        contr="trt.vs.ctrl", 
                        type="response")
        
        
        # Data frame with ls means
        means_df <- comp$emmeans %>% 
                data.frame() %>% 
                dplyr::rename(delta_ct=response)
        
        # Data frame with comparisons
        contr_df <- comp$contrasts %>% 
                data.frame() %>% 
                mutate(p=round(p.value, 3))
        
        contrasts[[i]] <- contr_df
        
        # Mark the significant comparisons
        sig_star <- vector(length=length(contr_df$p.value))
        
        for(k in 1:length(contr_df$p.value)){
                
                # print(contr_df$p.value[k])
                
                if(contr_df$p.value[k]<0.05){
                        sig_star[k] <- "*"
                } 
                
                else{
                        sig_star[k] <- ""
                }
                
        }
        
        means_df$sig_star <- c("",sig_star)
        
        means[[i]] <- means_df
}


#------------------------------------------------------------------------------
#--------------------- Tables with means and CI -------------------------------
#------------------------------------------------------------------------------


means_data_frame <- ldply(means, .fun=data.frame) %>% 
        dplyr::rename(Gen=.id, Treatment=Treatment) %>% 
        data.frame()

write.csv2(means_data_frame, "gene_expression_means_dilger.csv", row.names=FALSE)


#------------------------------------------------------------------------------
#--------------- Table with Contrasts and p.values ----------------------------
#------------------------------------------------------------------------------

contrasts_data_frame <- ldply(contrasts, .fun=data.frame) %>% 
        dplyr::rename(Gene=.id) %>% 
        data.frame()

write.csv2(contrasts_data_frame, "gene_expression_contrasts_dilger.csv", row.names=FALSE)


#------------------------------------------------------------------------------
#-------------------------- Grafic 1 B ----------------------------------------
#------------------------------------------------------------------------------

# Filter original data for the genes
dat_1 <- dat %>% 
        dplyr::filter(Gen=="NT5E"| 
                              Gen=="THY1"|
                              Gen=="ENG"|
                              Gen=="ALCAM"|
                              Gen=="NES"|
                              Gen=="TUBB3"|
                              Gen=="RBFOX3"|
                              Gen=="MAP2"|
                              Gen=="POU3F2"|
                              Gen=="SOX2"|
                              Gen=="MYT1L") %>% 
        dplyr::rename(Treatment=Treatment)

# Relevel the genes (for facet order)
dat_1$Gen <- factor(dat_1$Gen, 
                      levels=c("NT5E", "THY1", "ENG", "ALCAM", "NES",
                               "TUBB3", "RBFOX3", "MAP2", "POU3F2",
                               "SOX2", "MYT1L"))

# Filter the means for the genes
means_1 <- means_data_frame %>% 
        dplyr::filter(Gen=="NT5E"| 
                              Gen=="THY1"|
                              Gen=="ENG"|
                              Gen=="ALCAM"|
                              Gen=="NES"|
                              Gen=="TUBB3"|
                              Gen=="RBFOX3"|
                              Gen=="MAP2"|
                              Gen=="POU3F2"|
                              Gen=="SOX2"|
                              Gen=="MYT1L")

# Find the max. delta ct values (for aligning the stars in the grafic)
max_delta_1 <- dat_1 %>% 
        group_by(Gen, Treatment) %>% 
        dplyr::summarize(max_delta=max(delta_ct, na.rm=TRUE)) 

# Join the means and the max_delta data setws
means_1 <- full_join(means_1, max_delta_1) 

# Relevel the genes (for facet order)
means_1$Gen <- factor(means_1$Gen, 
                      levels=c("NT5E", "THY1", "ENG", "ALCAM", "NES",
                               "TUBB3", "RBFOX3", "MAP2", "POU3F2",
                               "SOX2", "MYT1L"))

# Plot the grafic
ggplot(dat_1, aes(x=Treatment, y=delta_ct))+
        theme_bw()+
        geom_point(aes(color=Treatment),
                   alpha=0.6,
                   position=position_jitter(height=0, width=0.1))+
        geom_linerange(data=means_1, aes(x=Treatment,
                                           ymin=lower.CL,
                                           ymax=upper.CL))+
        geom_point(data=means_1, aes(x=Treatment, y=delta_ct), shape="-", size=6)+
        geom_text(data=means_1, aes(y=max_delta,
                                    x=Treatment,
                                    label=sig_star), size=6)+
        scale_y_continuous(breaks=c(1e-8, 1e-06, 1e-04, 1e-02, 1), 
                           trans="log10")+
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
        ylab(bquote("Relative gene expression ("~2^{-Delta~"Ct"}~")"))+
        xlab("")+
        facet_grid(~Gen)

# Save the grafic
ggsave("gene_expression_fig_1.png", width=13*2, height=9*1.5, units="cm", dpi=600)

#------------------------------------------------------------------------------
#----------------------------- Grafic 3 A -------------------------------------
#------------------------------------------------------------------------------

# Filter original data for the genes
dat_3 <- dat %>% 
        dplyr::filter(Gen=="GJB2"| 
                              Gen=="GJD2"|
                              Gen=="GJA4"|
                              Gen=="GJA5"|
                              Gen=="GJA1"|
                              Gen=="GJC1") %>% 
        dplyr::rename(Treatment=Treatment)

# Relevel the genes (for facet order)
dat_3$Gen <- factor(dat_3$Gen, 
                    levels=c("GJB2", "GJD2", "GJA4", "GJA5", "GJA1", "GJC1"))

# Filter the means for the genes
means_3 <- means_data_frame %>% 
        dplyr::filter(Gen=="GJB2"| 
                              Gen=="GJD2"|
                              Gen=="GJA4"|
                              Gen=="GJA5"|
                              Gen=="GJA1"|
                              Gen=="GJC1")

# Find the max. delta ct values (for aligning the stars in the grafic)
max_delta_3 <- dat_3 %>% 
        group_by(Gen, Treatment) %>% 
        dplyr::summarize(max_delta=max(delta_ct, na.rm=TRUE)) 

# Join the means and the max_delta data setws
means_3 <- full_join(means_3, max_delta_3) 

# Relevel the genes (for facet order)
means_3$Gen <- factor(means_3$Gen, 
                      levels=c("GJB2", "GJD2", "GJA4", "GJA5", "GJA1", "GJC1"))

# Plot the grafic
ggplot(dat_3, aes(x=Treatment, y=delta_ct))+
        theme_bw()+
        geom_point(aes(color=Treatment),
                   alpha=0.6,
                   position=position_jitter(height=0, width=0.1))+
        geom_linerange(data=means_3, aes(x=Treatment,
                                         ymin=lower.CL,
                                         ymax=upper.CL))+
        geom_point(data=means_3, aes(x=Treatment, y=delta_ct), shape="-", size=6)+
        geom_text(data=means_3, aes(y=max_delta,
                                    x=Treatment,
                                    label=sig_star), size=6)+
        scale_y_continuous(breaks=c(1e-8, 1e-06, 1e-04, 1e-02, 1), 
                           trans="log10")+
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
        ylab(bquote("Relative gene expression ("~2^{-Delta~"Ct"}~")"))+
        xlab("")+
        facet_grid(~Gen)

# Save the grafic
ggsave("gene_expression_fig_3.png", width=13*2, height=9*1.5, units="cm", dpi=600)


#------------------------------------------------------------------------------
#--------------- Gene expression Appendix -------------------------------------
#------------------------------------------------------------------------------

# Filter original data for the genes
dat_A <- dat %>% 
        dplyr::filter(Gen=="ACTA2"| 
                              Gen=="TAGLN"|
                              Gen=="MYL2"|
                              Gen=="SMYD1") %>% 
        dplyr::rename(Treatment=Treatment)

# Relevel the genes (for facet order)
dat_A$Gen <- factor(dat_A$Gen, 
                    levels=c("ACTA2", "TAGLN", "MYL2", "SMYD1"))

# Filter the means for the genes
means_A <- means_data_frame %>% 
        dplyr::filter(Gen=="ACTA2"| 
                              Gen=="TAGLN"|
                              Gen=="MYL2"|
                              Gen=="SMYD1")

# Find the max. delta ct values (for aligning the stars in the grafic)
max_delta_A <- dat_A %>% 
        group_by(Gen, Treatment) %>% 
        dplyr::summarize(max_delta=max(delta_ct, na.rm=TRUE)) 

# Join the means and the max_delta data setws
means_A <- full_join(means_A, max_delta_A) 

# Relevel the genes (for facet order)
means_A$Gen <- factor(means_A$Gen, 
                      levels=c("ACTA2", "TAGLN", "MYL2", "SMYD1"))

# Plot the grafic
ggplot(dat_A, aes(x=Treatment, y=delta_ct))+
        theme_bw()+
        geom_point(aes(color=Treatment),
                   alpha=0.6,
                   position=position_jitter(height=0, width=0.1))+
        geom_linerange(data=means_A, aes(x=Treatment,
                                         ymin=lower.CL,
                                         ymax=upper.CL))+
        geom_point(data=means_A, aes(x=Treatment, y=delta_ct), shape="-", size=6)+
        geom_text(data=means_A, aes(y=max_delta,
                                    x=Treatment,
                                    label=sig_star), size=6)+
        scale_y_continuous(breaks=c(1e-8, 1e-06, 1e-04, 1e-02, 1), 
                           limits=c(1e-08, 1),
                           trans="log10")+
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
        ylab(bquote("Relative gene expression ("~2^{-Delta~"Ct"}~")"))+
        xlab("")+
        facet_grid(~Gen)

# Save the grafic
ggsave("gene_expression_fig_appendix.png", width=13*1.5, height=9*1.5, units="cm", dpi=600)


























