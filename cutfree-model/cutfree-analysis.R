library(ggplot2)

plota <- ggplot(runtime_data)+
  geom_point(aes(x=OligoLength, y=CutFree_Time))+
  theme_classic()+
  xlab("Oligo Length")+
  ylab("Time (s)")

plota

modela <- lm(CutFree_Time ~ OligoLength, data=runtime_data)
summary(modela)

plota_RL <- ggplot(runtime_data)+
  geom_point(aes(x=OligoLength, y=CutFreeRL_Time))+
  theme_classic()+
  xlab("Oligo Length")+
  ylab("Time (s)")

plota_RL

modela_RL <- lm(CutFreeRL_Time ~ OligoLength, data=runtime_data)
summary(modela_RL)

plotb <- ggplot(runtime_data)+
  geom_point(aes(x=SiteLength, y=CutFreeRL_Time))+
  theme_classic()+
  xlab("Site Length")+
  ylab("Time (s)")

plotb

modelb <- lm(CutFreeRL_Time ~ SiteLength, data=runtime_data)
summary(modelb)

plotc <- ggplot(runtime_data)+
  geom_point(aes(x=TotalSites, y=CutFreeRL_Time))+
  theme_classic()+
  xlab("Total Sites")+
  ylab("Time (s)")

plotc

modelc <- lm(CutFreeRL_Time ~ TotalSites, data=runtime_data)
summary(modelc)

plot1 <- ggplot(runtime_data)+
  geom_point(aes(x=OligoLength+SiteLength+TotalSites, y=CutFreeRL_Time))+
  theme_classic()+
  xlab("Oligo Length")+
  ylab("Time (s)")

plot1

model1 <- lm(CutFreeRL_Time ~ OligoLength+SiteLength+TotalSites, data=runtime_data)
summary(model1)

plot2 <- ggplot(runtime_data)+
  geom_point(aes(x=log(SiteLength)+OligoLength*SiteLength*TotalSites, y=CutFree_Time))+
  theme_classic()+
  xlab("Oligo Length")+
  ylab("Time (s)")

plot2

model2 <- lm(CutFree_Time ~ log(SiteLength)+OligoLength*SiteLength*TotalSites, data=runtime_data)
summary(model2)

plot2_RL <- ggplot(runtime_data)+
  geom_point(aes(x=log(SiteLength)+OligoLength*SiteLength*TotalSites, y=CutFreeRL_Time))+
  theme_classic()+
  xlab("Oligo Length")+
  ylab("Time (s)")

plot2_RL

model2_RL <- lm(CutFreeRL_Time ~ log(SiteLength)+OligoLength*SiteLength*TotalSites, data=runtime_data)
summary(model2_RL)

plot3 <- ggplot(runtime_data)+
  geom_point(aes(x=log(SiteLength)+log(OligoLength)*log(SiteLength)*log(TotalSites), y=CutFreeRL_Time))+
  theme_classic()+
  xlab("Oligo Length")+
  ylab("Time (s)")

plot3

model3 <- lm(CutFreeRL_Time ~ log(SiteLength)+log(OligoLength)*log(SiteLength)*log(TotalSites), data=runtime_data)
summary(model3)

plot4 <- ggplot(runtime_data)+
  geom_point(aes(x=log(SiteLength)+log(OligoLength)*log(SiteLength)*log(TotalSites), y=CutFree_Time))+
  theme_classic()+
  xlab("Oligo Length")+
  ylab("Time (s)")

plot4

model4 <- lm(CutFree_Time ~ log(SiteLength)+log(OligoLength)*log(SiteLength)*log(TotalSites), data=runtime_data)
summary(model4)


