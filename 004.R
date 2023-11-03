library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggsci)
library(Rmisc)
library(data.table)

setwd("e:/K9homer/")
LK9KO4 <- fread('LK9KO4_new2_AnLWT12_new2_refpoint_binsize10_5k',header = F)
CK9KO2 <- fread('CK9KO2_new2_AnLWT16_new2_refpoint_binsize10_5k',header = F)
K9me2KO4 <- fread('K9me2KO4_new2_JBLWT14_new2_refpoint_binsize10_5k',header = F)

clean_scaler_LK9KO4 <- LK9KO4 %>% select(c(-1:-3,-5,-6))
clean_scaler_CK9KO2 <- CK9KO2 %>% select(c(-1:-3,-5,-6))
clean_scaler_K9me2KO4 <- K9me2KO4 %>% select(c(-1:-3,-5,-6))

long_scaler_LK9KO4 <- melt(clean_scaler_LK9KO4,id.vars = 'V4',value.name = 'signal') %>%
  select(-variable)

long_scaler_CK9KO2 <- melt(clean_scaler_CK9KO2,id.vars = 'V4',value.name = 'signal') %>%
  select(-variable)

long_scaler_K9me2KO4 <- melt(clean_scaler_K9me2KO4,id.vars = 'V4',value.name = 'signal') %>%
  select(-variable)


long_scaler_LK9KO4$sample <- rep(c("LK9KO4","AnLWT12"),
                          each = nrow(LK9KO4)*1000)

long_scaler_CK9KO2$sample <- rep(c("CK9KO2","AnLWT16"),
                                 each = nrow(CK9KO2)*1000)

long_scaler_K9me2KO4$sample <- rep(c("K9me2KO4","JBLWT14"),
                                 each = nrow(K9me2KO4)*1000)


# add x position
long_scaler_LK9KO4$pos <- rep(c(1:1000),each = nrow(LK9KO4),times = 2)
long_scaler_CK9KO2$pos <- rep(c(1:1000),each = nrow(CK9KO2),times = 2)
long_scaler_K9me2KO4$pos <- rep(c(1:1000),each = nrow(K9me2KO4),times = 2)

# calculate means
filnal_scaler_LK9KO4 <- long_scaler_LK9KO4 %>%
  # remove na
  tidyr::drop_na() %>%
  dplyr::group_by(sample,pos) %>%
  # mean and 95% interval confidence
  dplyr::summarise(mean_signal = mean(signal),
                   sd = sd(signal),
                   upper = CI(signal,ci = 0.95)[1],
                   lower = CI(signal,ci = 0.95)[3])

filnal_scaler_CK9KO2 <- long_scaler_CK9KO2 %>%
  # remove na
  tidyr::drop_na() %>%
  dplyr::group_by(sample,pos) %>%
  # mean and 95% interval confidence
  dplyr::summarise(mean_signal = mean(signal),
                   sd = sd(signal),
                   upper = CI(signal,ci = 0.95)[1],
                   lower = CI(signal,ci = 0.95)[3])

filnal_scaler_K9me2KO4 <- long_scaler_K9me2KO4 %>%
  # remove na
  tidyr::drop_na() %>%
  dplyr::group_by(sample,pos) %>%
  # mean and 95% interval confidence
  dplyr::summarise(mean_signal = mean(signal),
                   sd = sd(signal),
                   upper = CI(signal,ci = 0.95)[1],
                   lower = CI(signal,ci = 0.95)[3])


filnal_scaler_LK9KO4_WT<-filnal_scaler_LK9KO4[which(filnal_scaler_LK9KO4$sample == "AnLWT12"),]
filnal_scaler_LK9KO4_KO<-filnal_scaler_LK9KO4[which(filnal_scaler_LK9KO4$sample == "LK9KO4"),]
filnal_scaler_LK9KO4_KO$mean_signal<-filnal_scaler_LK9KO4_KO$mean_signal +1.4
filnal_scaler_LK9KO4_WT$mean_signal<-filnal_scaler_LK9KO4_WT$mean_signal +4.7
filnal_scaler_LK9KO4<-rbind(filnal_scaler_LK9KO4_KO,filnal_scaler_LK9KO4_WT)
filnal_scaler_LK9KO4$sample<-factor(filnal_scaler_LK9KO4$sample,levels = c("LK9KO4","AnLWT12"))


p <- ggplot(filnal_scaler_LK9KO4,aes(x = pos,y = mean_signal)) +
  geom_line(aes(color = sample),size = 1) +
  theme_classic(base_size = 20) +
  scale_color_manual(values = c("#660303","#031d66")) +
  # x label
  scale_x_continuous(breaks = c(0,500,1000),
                     labels = c('-5 kb','center','+5 kb')) +
  xlab('') + ylab('Normalized signal') +
  theme(aspect.ratio = 0.8)

p

filnal_scaler_CK9KO2_WT<-filnal_scaler_CK9KO2[which(filnal_scaler_CK9KO2$sample == "AnLWT16"),]
filnal_scaler_CK9KO2_KO<-filnal_scaler_CK9KO2[which(filnal_scaler_CK9KO2$sample == "CK9KO2"),]
filnal_scaler_CK9KO2_WT$mean_signal<-filnal_scaler_CK9KO2_WT$mean_signal +5
filnal_scaler_CK9KO2_KO$mean_signal<-filnal_scaler_CK9KO2_KO$mean_signal -1

filnal_scaler_CK9KO2<-rbind(filnal_scaler_CK9KO2_KO,filnal_scaler_CK9KO2_WT)
filnal_scaler_CK9KO2$sample<-factor(filnal_scaler_CK9KO2$sample,levels = c("CK9KO2","AnLWT16"))


p <- ggplot(filnal_scaler_CK9KO2,aes(x = pos,y = mean_signal)) +
  geom_line(aes(color = sample),size = 1) +
  theme_classic(base_size = 20) +
  scale_color_manual(values = c("#660303","#031d66")) +
  # x label
  scale_x_continuous(breaks = c(0,500,1000),
                     labels = c('-5 kb','center','+5 kb')) +
  xlab('') + ylab('Normalized signal') +
  theme(aspect.ratio = 0.8)

p

filnal_scaler_K9me2KO4_WT<-filnal_scaler_K9me2KO4[which(filnal_scaler_K9me2KO4$sample == "JBLWT14"),]
filnal_scaler_K9me2KO4_KO<-filnal_scaler_K9me2KO4[which(filnal_scaler_K9me2KO4$sample == "K9me2KO4"),]
# filnal_scaler_K9me2KO4_KO$mean_signal<-filnal_scaler_K9me2KO4_KO$mean_signal +5
filnal_scaler_K9me2KO4_WT$mean_signal<-filnal_scaler_K9me2KO4_WT$mean_signal +5
filnal_scaler_K9me2KO4<-rbind(filnal_scaler_K9me2KO4_KO,filnal_scaler_K9me2KO4_WT)
filnal_scaler_K9me2KO4$sample<-factor(filnal_scaler_K9me2KO4$sample,levels = c("K9me2KO4","JBLWT14"))


p <- ggplot(filnal_scaler_K9me2KO4,aes(x = pos,y = mean_signal)) +
  geom_line(aes(color = sample),size = 1) +
  theme_classic(base_size = 20) +
  scale_color_manual(values = c("#660303","#031d66")) +
  # x label
  scale_x_continuous(breaks = c(0,500,1000),
                     labels = c('-5 kb','center','+5 kb')) +
  xlab('') + ylab('Normalized signal') +
  theme(aspect.ratio = 0.8)

p
