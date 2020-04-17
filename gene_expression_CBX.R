#------------------------------------------------------------------------------
#---------------------- qPCR DCB Nadine Dilger ------------------------------------
#------------------------------------------------------------------------------

# Loading required packages
library(hnp)
library(tidyverse)
library(plyr)
library(emmeans)

# import the data set
dat <- read.csv2("https://raw.githubusercontent.com/MaxMenssen/Dilger_et_al_2020/master/gene_expression_CBX_dilger.csv")
str(dat)

dat$Run <- as.factor(dat$Run)

# Separate Tretment and Concentration
dat <- separate(dat, Behandlung, c("Treat", "Conc"), sep="C") %>% 
        mutate(Treat=factor(Treat),
               Conc=factor(Conc, levels=c("0", "50", "100")))


# Overview about the data
head(dat)
str(dat)




#------------------------------------------------------------------------------

# Graphical overview about the exp. design
# Show which Genes were on the same plate 

ggplot(dat, aes(x=Treat, y=delta_ct))+
        theme_bw()+
        geom_point()+
        facet_grid(Run~Gen)+
        xlab("Treatment")


#------------------------------------------------------------------------------

# Contrast matrix for user defined contrasts 
# (all ttreatments against MSC_0 and AC_50, AC_100 vs. AC_0 )

AC_0 <- c(1,0,0,0,0,0)
MSC_0 <- c(0,1,0,0,0,0)
AC_50 <- c(0,0,1,0,0,0)
MSC_50 <- c(0,0,0,1,0,0) 
AC_100 <- c(0,0,0,0,1,0)
MSC_100 <- c(0,0,0,0,0,1)

contrast_list <- list("MSC_0 / MSC_50"=MSC_0-MSC_50,
                      "MSC_0 / MSC_100"=MSC_0-MSC_100,
                      "MSC_0 / AC_0"=MSC_0-AC_0,
                      "MSC_0 / AC_50"=MSC_0-AC_50,
                      "MSC_0 / AC_100"=MSC_0-AC_100,
                      "AC_0 / AC_50"=AC_0-AC_50,
                      "AC_0 / AC_100"=AC_0-AC_100)

contrast_list


#------------------------------------------------------------------------------

# Vector with names of the genes
genes <- levels(dat$Gen)

# Open vectors to store results
anova_results <- vector("list", length=length(genes))
names(anova_results) <- genes

res_vs_fitted <- vector("list", length=length(genes))
names(res_vs_fitted) <- genes

hnp_plot <- vector("list", length=length(genes))
names(hnp_plot) <- genes

means <- vector("list", length=length(genes))
names(means) <- genes

contrasts <- vector("list", length=length(genes))
names(contrasts) <- genes

for(i in 1:length(genes)){
        
        # Filter for the genes
        dat_fitered <- dat %>% 
                filter(Gen==genes[i])
        
        # Fit the linear model 
        fit <- lm(log(delta_ct)~Run+Treat*Conc, dat_fitered)
        
        # ANOVA Results
        # print(anova_results[[i]])
        
        anova_res <- anova(fit) %>% 
                data.frame() %>% 
                mutate(p=round(Pr..F., 3))
        
        anova_res$predictor <- row.names(anova_res)
        
        anova_res$p[anova_res$p < 0.001] <- "<0.001"
        
        anova_results[[i]] <- anova_res
        
        # print(anova_results[[i]])
        
        # Store the residual vs. fitted plot
        res_vs_fitted[[i]] <- plot(fit, which=1)
        
        # Normal plot
        hnp_plot[[i]] <- hnp(fit, halfnormal=FALSE,
                             print.on=TRUE, paint.out=TRUE)
        
        # Extracting the means
        comp <- emmeans(fit, ~Treat*Conc, type="response")
        
        # Data frame with ls means
        means_df <- comp %>% 
                data.frame() %>% 
                dplyr::rename(delta_ct=response)
        
        # print(means_df)
        
        # Mean Comparisons
        contr <- contrast(comp, method=contrast_list, type="response")
        
        # Data frame with comparisons
        contr_df <- contr %>% 
                as.data.frame() %>% 
                mutate(p=round(p.value, 3))
        
        # print(contr_df)
        
        contrasts[[i]] <- contr_df
        
        
        # Mark the significant comparisons for msc 0
        
        contr_df_msc <- contr_df[1:5,] 
        
        sig_star_msc <- vector(length=5)
        
        for(k in 1:length(sig_star_msc)){
                
                # print(contr_df$p.value[k])
                
                if(contr_df_msc$p.value[k]<0.05){
                        sig_star_msc[k] <- "*"
                                } 
                
                else{
                        sig_star_msc[k] <- ""
                }
                
                # print(contr_df_msc$p.value[k])
                # print(sig_star_msc[k])
        }
        
        
        # Significant comparisons ac 0
        contr_df_ac <- contr_df[6:7,] 
        
        # print(contr_df_ac)
        
        sig_star_ac <- vector(length=2)
        
        for(l in 1:length(sig_star_ac)){
                
                if(contr_df_ac$p.value[l]<0.05){
                        sig_star_ac[l] <- "#"
                } 
                
                else{
                        sig_star_ac[l] <- ""
                }
                
                # print(contr_df_ac$p.value[l])
                # print(sig_star_ac[l])
        }
        
        # ACHTUNG: Die Zuordnung klappt noch nicht!!!
        # means_df muss nach MS 0 geordnet werden
        
        # Labels für ANOVA-Ergebnisse einbauen!
        
        # means_df$Conc_sort <- c("A", "C", "B", "D", "F", "E")
        
        means_df <- means_df %>% 
                dplyr::arrange(Treat=desc(Treat))
        
        means_df$sig_star_msc <- c("", sig_star_msc)
        means_df$sig_star_ac <- c("", "", "", "", sig_star_ac)
        
        # print(means_df)
        
        means[[i]] <- means_df
}

#------------------------------------------------------------------------------

anova_data_frame <- ldply(anova_results, .fun=data.frame) %>% 
        dplyr::rename(Gen=.id) %>% 
        select(Gen, Df, predictor, everything()) %>% 
        data.frame()

write.csv2(anova_data_frame, "genes_anova_CBX_dilger.csv", 
           row.names=FALSE)

#------------------------------------------------------------------------------


means_data_frame <- ldply(means, .fun=data.frame) %>% 
        dplyr::rename(Gen=.id) %>% 
        data.frame()

write.csv2(means_data_frame, "genes_means_CBX_dilger.csv", 
           row.names=FALSE)

#------------------------------------------------------------------------------


contrasts_data_frame <- ldply(contrasts, .fun=data.frame) %>% 
        dplyr::rename(Gene=.id) %>% 
        data.frame()

write.csv2(contrasts_data_frame, "genes_contrasts_CBX_dilger.csv", 
           row.names=FALSE)


#------------------------------------------------------------------------------
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
        dplyr::rename(Treatment=Treat,
                      Concentration=Conc)

# Relevel the genes (for facet order)
dat_1$Gen <- factor(dat_1$Gen, 
                    levels=c("NT5E", "THY1", "ENG", "ALCAM", "NES",
                             "TUBB3", "RBFOX3", "MAP2", "POU3F2",
                             "SOX2", "MYT1L"))
str(dat_1)

dat_1$Interaction <- dat_1$Treat : dat_1$Conc
dat_1$Interaction <- factor(dat_1$Interaction, 
                            levels=c("MS:0", "MS:50", "MS:100",
                                     "A:0", "A:50" , "A:100"))


# dat_1$inter <- droplevels(dat_1$Treatment : dat_1$Concentration)

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
                              Gen=="MYT1L")%>% 
        dplyr::rename(Treatment=Treat,
                      Concentration=Conc)

# Find the max. delta ct values (for aligning the stars in the grafic)
max_delta_1 <- dat_1 %>% 
        group_by(Gen, Treatment, Concentration) %>% 
        dplyr::summarize(max_delta=max(delta_ct, na.rm=TRUE)) 

# Join the means and the max_delta data setws
means_1 <- full_join(means_1, max_delta_1) 

# Relevel the genes (for facet order)
means_1$Gen <- factor(means_1$Gen, 
                      levels=c("NT5E", "THY1", "ENG", "ALCAM", "NES",
                               "TUBB3", "RBFOX3", "MAP2", "POU3F2",
                               "SOX2", "MYT1L"))

means_1$Interaction <- means_1$Treat : means_1$Conc
means_1$Interaction <- factor(means_1$Interaction, 
                            levels=c("MS:0", "MS:50", "MS:100",
                                     "A:0", "A:50" , "A:100"))

# Labels for the anova results
anova_labs_1 <- anova_data_frame %>% 
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

# Plot the grafic
ggplot(dat_1, aes(x=Interaction, y=delta_ct))+
        theme_bw()+
        geom_point(aes(color=Interaction),
                   alpha=0.6,
                   position=position_jitter(height=0, width=0.1))+
        geom_linerange(data=means_1, aes(x=Interaction,
                                         ymin=lower.CL,
                                         ymax=upper.CL))+
        geom_point(data=means_1, aes(x=Interaction,
                                     y=delta_ct),
                   shape="-", size=6)+
        geom_text(data=means_1, aes(y=1.2*max_delta,
                                    x=Interaction,
                                    label=sig_star_msc), size=6)+
        geom_text(data=means_1, aes(y=2.7*max_delta,
                                    x=Interaction,
                                    label=sig_star_ac), size=3)+
        scale_y_continuous(breaks=c(1e-8, 1e-06, 1e-04, 1e-02, 1, 10), 
                           trans="log10")+
        scale_color_discrete(breaks=c("MS:0",
                                      "MS:50",
                                      "MS:100",
                                      "A:0",
                                      "A:50",
                                      "A:100"),
                             labels=c(bquote("MSC, 0"~mu*"M"),
                                      bquote("MSC, 50"~mu*"M"),
                                      bquote("MSC, 100"~mu*"M"),
                                      bquote("AC, 0"~mu*"M"),
                                      bquote("AC, 50"~mu*"M"),
                                      bquote("AC, 100"~mu*"M")))+
        scale_color_manual(values=c("#FFCC00",
                                    "#FF9900",
                                    "#FF6600",
                                    "#0099FF",
                                    "#66CCFF",
                                    "#3399CC"))+
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
ggsave("gene_expression_CBX_fig_1.png", width=13*2, height=9*1.5, units="cm", dpi=600)











