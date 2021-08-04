library(ggstatsplot)
library(tidyverse)
library(dplyr)
library(forcats)
library(PairedData)
library(ggpubr)
library(rstatix)

mmdata <- read.csv("D:\\Dropbox\\Conservation\\Shark Stewards\\MethylMercury\\data.csv", header=TRUE, stringsAsFactors=FALSE) %>%
  filter(market_category != "Unknown" & market_category != "Whole Fish" & market_category != "File Fish")

fda <- read.csv("D:\\Dropbox\\Conservation\\Shark Stewards\\MethylMercury\\fda.csv", header=TRUE, stringsAsFactors=FALSE) %>%
  group_by(market_category) %>%
  summarise(mean_ppm = mean(ppm), median_ppm = median(ppm))

# Collapse MMData into the market categories, summarizing by SD, Mean, and N
mkcat <- mmdata %>%
  group_by(market_category) %>% 
  summarise(sd_ppm = sd(mercury_ppm), mean_ppm = mean(mercury_ppm), n = n())

#Add the mean in each Market Category for the US Samples
mkcat$us_ppm <- mmdata %>%
  group_by(market_category, country) %>%
  summarise(country_mean = mean(mercury_ppm)) %>%
  filter(country == "US") %>%
  pull(country_mean)

#Add the mean in each Market Category for the UK Samples
mkcat$hk_ppm <- mmdata %>%
  group_by(market_category, country) %>%
  summarise(country_mean = mean(mercury_ppm)) %>%
  filter(country == "HK") %>%
  pull(country_mean)

#Calculate the 95% confidence Intervals
mkcat$interval <- qt(0.975,df=mkcat$n-1) * ( mkcat$sd_ppm / sqrt(mkcat$n))
mkcat$lower <- mkcat$mean_ppm - mkcat$interval
mkcat$upper <- mkcat$mean_ppm + mkcat$interval

#Summarize by Trophic Levels
trophic_level <- mmdata %>%
  filter(trophic_level != "") %>%
  group_by(trophic_level) %>%
  summarise(mean_ppm = mean(mercury_ppm), median_ppm = median(mercury_ppm), n = n())

#Summarize by Habitat Level
habitat <- mmdata %>%
  filter(habitat != "") %>%
  group_by(habitat) %>%
  summarise(mean_ppm = mean(mercury_ppm), median_ppm = median(mercury_ppm), n = n())


#===================================
# Analysis

#Rstudio T-test of the means

# Create difference between means of HK vs SF
mkcat <- mkcat %>% mutate(country_differences = hk_ppm - us_ppm)
mkcat %>% identify_outliers(country_differences)
mkcat %>% shapiro_test(country_differences) 
ggqqplot(mkcat, "country_differences")

# Calculate the T-test
t.test(mkcat$us_ppm, mkcat$hk_ppm, paried = TRUE)

# T-test Methalated vs Organic
totalmercury <- filter(mmdata, lab == "Eurofins HK" | lab == "Eurofins LA" && id %in% methylated$id)
methylated <- filter(mmdata, lab == "Wisconsin State")
t.test(totalmercury$mercury_ppm, methylated$mercury_ppm)

# Paird T-test FDA / dried fish
dried_fdacat <- filter(mkcat, market_category %in% fda$market_category)

t.test(dried_fdacat$mean_ppm, fda$mean_ppm, paired = TRUE, alternative = "greater")

#Two Way ANOVA
interaction <- aov(mercury_ppm ~ country * market_category, data = mmdata)
summary(interaction)

# Build the linear model
lmodel  <- lm(mercury_ppm ~ trophic_level,
             data = filter(mmdata, trophic_level != ""))
# Create a QQ plot of residuals
ggqqplot(residuals(lmodel))
shapiro_test(residuals(lmodel))

# Kruskal Wallace Mercury_ppm country
country.kruskal <- filter(mmdata, trophic_level != "") %>% kruskal_test(mercury_ppm ~ country)
country.kruskal
mmdata %>% kruskal_effsize(mercury_ppm ~ country)

# Kruskal Wallace Mercury_ppm country
trophic.kruskal <- filter(mmdata, trophic_level != "") %>% kruskal_test(mercury_ppm ~ trophic_level)
trophic.kruskal
mmdata %>% kruskal_effsize(mercury_ppm ~ trophic_level)

# Kruskal Wallace Mercury_ppm country
habitat.kruskal <- filter(mmdata, habitat != "") %>% kruskal_test(mercury_ppm ~ habitat)
habitat.kruskal
mmdata %>% kruskal_effsize(mercury_ppm ~ trophic_level)

#Two Way ANOVA
model <- aov(mercury_ppm ~ trophic_level * habitat, data = filter(mmdata, trophic_level != ""))
summary(model)

#===================================
# Graphing

# Box Plot with outliers labeled
ggbetweenstats(data = mmdata, 
               x = country,
               y = mercury_ppm,
               outlier.tagging = TRUE,
               outlier.label = market_category)

# Box Plot for trophic categorieswith outliers labeled
ggbetweenstats(data = filter(mmdata, trophic_level != ""), 
               x = trophic_level,
               y = mercury_ppm,
               outlier.tagging = TRUE,
               outlier.label = market_category)

# Box Plot for habitat categorieswith outliers labeled
ggbetweenstats(data = filter(mmdata, habitat != ""), 
               x = habitat,
               y = mercury_ppm,
               outlier.tagging = TRUE,
               outlier.label = market_category)

# Box Plot for Dried vs Wet
drywet <- paired(dried_fdacat$mean_ppm, fda$mean_ppm)
ggpaired(drywet, cond1 = "dried_fdacat$mean_ppm", cond2 = "fda$mean_ppm",
         fill = "condition", palette = "d")

# Graphed 95% confidence intervals for each market category. lul.
ggplot(mkcat, aes(x=mean_ppm, y=market_category)) +
  geom_bar(stat = 'identity') +
  geom_errorbarh(aes(xmin=mean_ppm-interval, xmax=mean_ppm+interval), width=.2, position=position_dodge(.9)) 


# Market Categories scaterplot with means
ggplot() +
  geom_point(data=mkcat, 
             aes(mean_ppm, fct_reorder(market_category,mean_ppm)), 
             col="forestgreen", 
             shape=21, 
             size=6, 
             fill = "forestgreen", 
             alpha=.5) +  
  geom_point(data=mmdata, 
             aes(mercury_ppm, market_category), 
             shape = 21, 
             alpha = 1, 
             size = 3, 
             fill = 'black', 
             colour = 'transparent') +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )

#  Trophic Levels scaterplot with medians
ggplot() +
  geom_point(data=trophic_level, 
             aes(median_ppm, fct_reorder(trophic_level,median_ppm)), 
             col="forestgreen", 
             shape=21, 
             size=6, 
             fill = "forestgreen", 
             alpha=.5) +
  geom_point(data=filter(mmdata, trophic_level != ""), 
             aes(mercury_ppm, trophic_level), 
             shape = 21, 
             alpha = 1, 
             size = 3, 
             fill = 'black', 
             colour = 'transparent') +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )

#  Trophic Levels scaterplot with medians
ggplot() +
  geom_point(data=habitat, 
             aes(median_ppm, fct_reorder(habitat,median_ppm)), 
             col="forestgreen", 
             shape=21, 
             size=6, 
             fill = "forestgreen", 
             alpha=.5) +
  geom_point(data=filter(mmdata, habitat != ""), 
             aes(mercury_ppm, habitat), 
             shape = 21, 
             alpha = 1, 
             size = 3, 
             fill = 'black', 
             colour = 'transparent') +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )


# Methyl vs Mercury
ggplot() +
  geom_point(data=mkcat, 
             aes(mean_ppm, fct_reorder(market_category,mean_ppm)), 
             col="forestgreen", 
             shape=21, 
             size=6, 
             fill = "forestgreen", 
             alpha=.1) +  
  geom_point(data=filter(mmdata, lab == "Eurofins HK" | lab == "Eurofins LA"), 
             aes(mercury_ppm, market_category), 
             shape = 21, 
             alpha = 1, 
             size = 3, 
             fill = 'darkred', 
             colour = 'transparent') +
  geom_point(data=filter(mmdata, lab == "Wisconsin State"), 
             aes(mercury_ppm, market_category), 
             shape = 21, 
             alpha = 1, 
             size = 3, 
             fill = 'orange', 
             colour = 'transparent') +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )

#Export Data Frame

#Add the median in each Market Category for the US Samples
mkcat$us_median_ppm <- mmdata %>%
  group_by(market_category, country) %>%
  summarise(country_median = median(mercury_ppm)) %>%
  filter(country == "US") %>%
  pull(country_median)

#Add the n in each Market Category for the US Samples
mkcat$us_n <- mmdata %>%
  group_by(market_category, country) %>%
  summarise(us_n = n()) %>%
  filter(country == "US") %>%
  pull(us_n)

#Add the median in each Market Category for the HK Samples
mkcat$hk_median_ppm <- mmdata %>%
  group_by(market_category, country) %>%
  summarise(country_median = median(mercury_ppm)) %>%
  filter(country == "HK") %>%
  pull(country_median)

#Add the n in each Market Category for the HK Samples
mkcat$hk_n <- mmdata %>%
  group_by(market_category, country) %>%
  summarise(us_n = n()) %>%
  filter(country == "HK") %>%
  pull(us_n)

write.csv(mkcat, "D:\\Dropbox\\Conservation\\Shark Stewards\\MethylMercury\\mkcat.csv")
write.csv(habitat, "D:\\Dropbox\\Conservation\\Shark Stewards\\MethylMercury\\habitat.csv")
write.csv(trophic_level, "D:\\Dropbox\\Conservation\\Shark Stewards\\MethylMercury\\trophic.csv")

