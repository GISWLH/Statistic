library(ggplot2)
library(ggrepel) 
library(maps)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
## whole map
world <- ne_countries(scale = "medium", returnclass = "sf")
site <- read.csv("D:/OneDrive/UCAS/stu/2023HanZY/2.csv")

ggplot() +
  geom_sf(data=world) +
#  geom_polygon(data = map_data('world'), aes(x = long, y = lat, group = group), fill = 'gray',) +
  theme_bw() +
  geom_point(data = site, aes(x = lon, y = lat, color = soil), size = 1.5) +
  scale_color_manual(values = c('#3CA6CC','#F99233')) +
  theme(panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = c(-150, -100, -50, 0, 50, 100, 150), expand = c(0, 0), 
                     labels = c('150°W', '100°W', '50°W', '0', '50°E', '100°E', '150°E')) +
  scale_y_continuous(breaks = c(-60, -30, 0, 30, 60), expand = c(0, 0), 
                     labels = c('60°S', '30°S', '0', '30°N', '60°N')) +
  labs(x = 'Longitude', y = 'Latitude', color = 'soil') +
  theme(axis.text = element_text(size=12,color="black", family="serif"),
        legend.text = element_text(size=12,color="black", family="serif"),
        legend.title = element_blank()) 

## sub map
ggplot() +
  geom_sf(data=world) +
  #  geom_polygon(data = map_data('world'), aes(x = long, y = lat, group = group), fill = 'gray',) +
  theme_bw() +
  geom_point(data = site, aes(x = lon, y = lat, color = soil), size = 6) +
  scale_color_manual(values = c('#3CA6CC','#F99233')) +
  theme(axis.ticks.x = element_blank())+
  theme(axis.ticks.y = element_blank())+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_blank())+
  coord_sf(xlim = c(100, 125), ylim = c(20, 40), expand = FALSE) 

## VPA
# 1.rda analysis
library(vegan)
library(tidyverse)
data(mite)
data(mite.env)
data(mite.pcnm)
rda_full <- rda(mite~., data = cbind(mite.pcnm, mite.env))
rda_null <- rda(mite~1, data = cbind(mite.pcnm, mite.env))

rda_back <- ordistep(rda_full, direction = 'backward',trace = 0)

# forward selection
rda_frwd <- ordistep(rda_null, formula(rda_full), direction = 'forward',trace = 0)

# bothward selection 
rda_both <- ordistep(rda_null, formula(rda_full), direction = 'both',trace = 0)


# 2.VPA analysis
# Divide the environmental data.frame by different categories.
df_env <- mite.env %>% select(WatrCont)
df_pcnm <- mite.pcnm %>% select(V2, V6)

# Perform variation partitioning analysis, the first variable is the community matrix
# the second and third variables are climate variable set and soil property variable set
vpt <- varpart(mite, df_env, df_pcnm)

vpt

plot(vpt, bg = 2:5, Xnames = c('env', 'pcnm'))

formula_env <- formula(mite ~ WatrCont + Condition(V2) + Condition(V6))
formula_pcnm <- formula(mite ~ Condition(WatrCont) + V2 + V6)

anova(rda(formula_env, data = cbind(mite.pcnm, mite.env)))
anova(rda(formula_pcnm, data = cbind(mite.pcnm, mite.env)))



mod <- varpart(mite, ~ SubsDens + WatrCont, ~ Substrate + Shrub + Topo,
               mite.pcnm, data=mite.env, transfo="hel")
plot(mod, bg=2:4, Xnames = c('SubsDens + WatrCont', 'Substrate + Shrub + Topo', 'pcnm'))
