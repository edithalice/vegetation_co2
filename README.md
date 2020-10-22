# Vegetation Cover and Its Effects on Local Air Quality
### By Edith Johnston

## Objective
The goal of this project was to try to model the effects of local surface vegetation cover on local air quality, beginning with carbon dioxide concentration, and potentially expanding to other pollutants if time permits. The idea was to have a model appicable to anywhere on the planet. I also wanted to build a autoregressive times series component into this model, for both the vegetation and CO2 levels - the idea being that, given a place which currently has x level of vegetation and y level of carbon dioxide, if we change x what happens to y.  

Other potential features I may try to add to the model include:
- Local geographic information, such as elevation, biome, or latitude (as a proxy for seasonality)
- Human related information, such as nearby population density, distance to closest urban or semi-urban area (and maybe a rank of said area by size), or perhaps country level-data such as clean air legislation or a ranking of countries by population's clean air habits (i.e. car driving vs public transport, etc)

Essentially the point of this project was to create a model that attempts to determine what, if any, is the effect of changing vegetation levels on CO2 levels at a local, rather than global level, with other factors held constant and accounting for the change of global CO2 levels due to climate change.

## Stretch Goal
Create an app where users can define vegetation levels in an area and view the resulting effect on air quality. Depending on what the model looks like, and what types of features I add, I may try to do this either by:
- providing a blank area for users to fill in with vegetation levels of their choice, then visualizing the resulting air quality
- providing a blank area based on an actual geographic location for users to fill with vegetation levels of their choice, then visualizing air quality
- providing the current vegetation levels of an actual geographic area for users to alter, then visualizing the resulting change in air quality

## Data Sources
- NASA's MODIS/Terra Vegetation Continuous Fields Yearly Dataset: "a global representation of surface vegetation cover as gradations of three ground cover components: percent tree cover, percent non-tree cover, and percent non-vegetated (bare)."
- NASA's TES Level 2 Dataset: concentration levels for a number of pollutant molecules, including carbon dioxide, carbon monoxide, methane, ozone, ammonia, and more.
- possibly more
  
# Work in Progress!
### **This repository will be updated with findings and code as soon as possible.**
