
### Analysis associated with Weber et al. publication: “Isotopic niche overlap among foraging marine turtle species in the Gulf of Mexico”

#codes for ranges, max and mins, means & SDs, and boxplots #

#load required packages
library(tidyverse)
library(gridExtra)
library(viridis)
library(hrbrthemes)
library(ggplot2)


#load data
CR_CNS <- read.csv("SIA_Values.csv")

summary(CR_CNS)

CR_CNS_factors<-c("Species") #select appropriate vars to convert into factors for CR_CNS data
CR_CNS[CR_CNS_factors]<-lapply(CR_CNS[CR_CNS_factors], as.factor) #convert selected vars to factors

CR_CNS_nums<-c("d15N",
             "d13C",
             "d34S")
CR_CNS[CR_CNS_nums]<-lapply(CR_CNS[CR_CNS_nums], as.numeric)

summary(CR_CNS) #check conversion

sd(CR_CNS$d13C)
sd(CR_CNS$d15N)
sd(CR_CNS$d34S)

CR_CNS %>% group_by(Species) %>%
  summarize(count = n(),
            mC = mean(d13C),
            sdC = sd(d13C), 
            mN = mean(d15N), 
            sdN = sd(d15N),
            mS = mean(d34S),
            sdS = sd(d34S))

CR_CNS %>% group_by(Species) %>%
  summarize(count = n(),
            rC = range(d13C),
            rN = range(d15N),
            rS = range(d34S))



# Boxplots #

NbpC <- ggplot(CR_CNS, aes(x = Species, y = d13C, fill = Species, color = Species)) +
  stat_boxplot(geom ="errorbar", width = 0.3, lwd = 1.3) + 
  geom_boxplot(fill = "white", outlier.size = 3, lwd = 1.3, fatten = 1) + 
  scale_color_manual(values = c("Cc" = "#E31A1C",
                                "Cm" = "#33A02C", 
                                "Lk" = "#0000FF")) +
  ###line below adds the mean, shown as a circle with a cross through it
  stat_summary(fun=mean, geom="point", shape=10, size=5, color="black") +
  theme_bw() + theme(legend.position="none", axis.line = element_line(color='black'),
                     plot.background = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.grid.major = element_blank(),
                     axis.text.x = element_text(color="black", size=16),
                     axis.text.y = element_text(size=16),
                     axis.title.x = element_text(size=20),
                     axis.title.y = element_text(size=20)) +
  xlab (NULL) + ylab(expression(delta^13* C * "  \u2030"))

NbpN <- ggplot(CR_CNS, aes(x = Species, y = d15N, fill = Species, color = Species)) +
  stat_boxplot(geom ="errorbar", width = 0.3, lwd = 1.3) + 
  geom_boxplot(fill = "white", outlier.size = 3, lwd = 1.3, fatten = 1) + 
  scale_color_manual(values = c("Cc" = "#E31A1C",
                                "Cm" = "#33A02C", 
                                "Lk" = "#0000FF")) +
  stat_summary(fun=mean, geom="point", shape=10, size=5, color="black") + 
  theme_bw() + theme(legend.position="none", axis.line = element_line(color='black'),
                     plot.background = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.grid.major = element_blank(),
                     axis.text.x = element_text(color="black", size=16),
                     axis.text.y = element_text(size=16),
                     axis.title.x = element_text(size=20),
                     axis.title.y = element_text(size=20)) +
  xlab (NULL) + ylab(expression(delta^15* N * "  \u2030"))

NbpS <- ggplot(CR_CNS, aes(x = Species, y = d34S, fill = Species, color = Species)) +
  stat_boxplot(geom ="errorbar", width = 0.3, lwd = 1.3) + 
  geom_boxplot(fill = "white", outlier.size = 3, lwd = 1.3, fatten = 1) + 
  scale_color_manual(values = c("Cc" = "#E31A1C",
                                "Cm" = "#33A02C", 
                                "Lk" = "#0000FF")) +
  stat_summary(fun=mean, geom="point", shape=10, size=5, color="black") + 
  theme_bw() + theme(legend.position="none", axis.line = element_line(color='black'),
                     plot.background = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.grid.major = element_blank(),
                     axis.text.x = element_text(color="black", size=16),
                     axis.text.y = element_text(size=16),
                     axis.title.x = element_text(size=20),
                     axis.title.y = element_text(size=20)) +
  xlab (NULL) + ylab(expression(delta^34* S * "  \u2030"))



library(grid)

bottom<-textGrob("Species", gp = gpar(fontsize=16))

grid.arrange(NbpC, NbpN, NbpS, nrow=1, ncol=3, bottom = bottom)



