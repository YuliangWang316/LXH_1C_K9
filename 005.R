library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggsci)
library(Rmisc)
library(data.table)

setwd("e:/K9homer/")
KO <- fread('merged_ALL_KO_WT_USE2_refpoint_binsize10_5k',header = F)

clean_scaler_KO <- KO %>% select(c(-1:-3,-5,-6))

long_scaler_KO <- melt(clean_scaler_KO,id.vars = 'V4',value.name = 'signal') %>%
  select(-variable)


long_scaler_KO$sample <- rep(c("KO","WT"),
                          each = nrow(KO)*1000)


# add x position
long_scaler_KO$pos <- rep(c(1:1000),each = nrow(KO),times = 2)

# calculate means
filnal_scaler_KO <- long_scaler_KO %>%
  # remove na
  tidyr::drop_na() %>%
  dplyr::group_by(sample,pos) %>%
  # mean and 95% interval confidence
  dplyr::summarise(mean_signal = mean(signal),
                   sd = sd(signal),
                   upper = CI(signal,ci = 0.95)[1],
                   lower = CI(signal,ci = 0.95)[3])



filnal_scaler_KO_WT<-filnal_scaler_KO[which(filnal_scaler_KO$sample == "WT"),]
filnal_scaler_KO_KO<-filnal_scaler_KO[which(filnal_scaler_KO$sample == "KO"),]
# filnal_scaler_KO_KO$mean_signal<-filnal_scaler_KO_KO$mean_signal +2
filnal_scaler_KO_WT$mean_signal<-filnal_scaler_KO_WT$mean_signal +19
filnal_scaler_KO<-rbind(filnal_scaler_KO_KO,filnal_scaler_KO_WT)
filnal_scaler_KO$sample<-factor(filnal_scaler_KO$sample,levels = c("KO","WT"))


p <- ggplot(filnal_scaler_KO,aes(x = pos,y = mean_signal)) +
  geom_line(aes(color = sample),size = 1) +
  theme_classic(base_size = 20) +
  scale_color_manual(values = c("#660303","#031d66")) +
  # x label
  scale_x_continuous(breaks = c(0,500,1000),
                     labels = c('-5 kb','center','+5 kb')) +
  xlab('') + ylab('Normalized signal') +
  theme(aspect.ratio = 0.8)

p

