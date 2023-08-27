

### Analysis associated with Weber et al. publication: “Isotopic niche overlap among foraging marine turtle species in the Gulf of Mexico”

# Study site map #

#load required packages
library(sf)
library(ggplot2)
library(ggspatial)
library(plotly)
library("MetBrewer")
library("patchwork")

###load data
Excel <- read.csv("MasterDB.csv")

#load shapefile of Florida
Florida<- st_read("FL.shp")


#First make the insert map of FL, make sure the rectangle encompasses the data points
ggplot(data=Florida)+
  geom_sf() +
  geom_point(data=Excel, aes(GPS_X, GPS_Y, color=Species))+
 scale_color_manual(values=c("#E31A1C", "#33A02C", "#0000FF")) +
  geom_rect(xmin=-82.9, xmax=-82.65, ymin=28.65, ymax=28.9, fill="transparent", color="red", size=1)+
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'none')

#same with the datapoints, make sure all are in frame
ggplot(data=Florida)+
  geom_sf() +
  geom_point(data=Excel, aes(GPS_X, GPS_Y, color=Species), size=3, alpha=1)+
  scale_color_manual(values=c("#E31A1C", "#33A02C", "#0000FF")) +
  lims(x=c(-82.9, -82.7), y=c(28.7, 28.87))+
  labs(x="Longitude", y="Latitude")+
  theme_bw()+
  theme(axis.title = element_text(size=16), 
        axis.text = element_text(size=12))+
  annotation_scale(location="bl", width_hint=0.4)+
  annotation_north_arrow(location="bl", which_north="true",
                         pad_x = unit(0.0, "in"), pad_y = unit(0.2, "in"),
                         style=north_arrow_fancy_orienteering)

# Calculate angle between river and x-axis
crystal_river <- st_coordinates(Florida[Florida$NAME == "Citrus",])
angle <- atan(diff(crystal_river[,2])/diff(crystal_river[,1])) * 180 / pi

#create inset map of Florida
p.inset <- ggplot(data=Florida)+
  geom_sf() +
  geom_sf_label(aes(label = "Florida"), size=4, fontface="bold") + #adds the "Florida" label on the inset map
   geom_point(data=Excel, aes(GPS_X, GPS_Y, color=Species))+
  scale_color_manual(values=c("#E31A1C", "#33A02C", "#0000FF")) +
  geom_rect(xmin=-83.1, xmax=-82.6, ymin=28.55, ymax=28.95, fill="red", size=1)+
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'none')

#create map showing turtle capture locations
p.main <- ggplot(data=Florida)+
  geom_sf() +
  geom_point(data=Excel, aes(GPS_X, GPS_Y, color=Species), size=3, alpha=1)+
  scale_color_manual(values=c("#E31A1C", "#33A02C", "#0000FF")) +
  lims(x=c(-82.93, -82.63), y=c(28.73, 28.94))+
  labs(x="Longitude", y="Latitude") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=16), 
        axis.text = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=14)) +
  annotation_scale(location="bl", width_hint=0.4)+
  annotation_north_arrow(location="bl", which_north="true",
                         pad_x = unit(0.0, "in"), pad_y = unit(0.2, "in"),
                         style=north_arrow_fancy_orienteering) +
annotate("text", x = -82.65, y = 28.92, label = "Crystal River", size = 4, fontface = "bold", angle = -20) #adding a label to Crystal River


#combine maps into one image
Fig1 <- p.main + inset_element(p.inset, left=0.01, top=0.98, right=0.4, bottom=0.55)


ggsave(filename = "Fig1.pdf", Fig1, width = 250, height = 234, units = "mm")

