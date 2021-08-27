# Vegetation Cover and Its Effects on Local Air Quality
### By Edith Johnston

## Table of Contents
1. [Objective](#objective)
2. [Data](#data)
4. [Model](#model)
5. [Main Tools Used](#main-tools-used)
6. [Deliverables](#deliverables)

## Objective
The goal of this project was to try to model the effects of local surface vegetation cover on local air quality, beginning with carbon dioxide concentration, and potentially expanding to other pollutants if time permits. The idea was to have a model appicable to anywhere on the planet. I also wanted to build a autoregressive times series component into this model, for both the vegetation and CO2 levels - the idea being that, given a place which currently has x level of vegetation and y level of carbon dioxide, if we change x what happens to y.  

Other potential features I would like to add to the model include:
- Local geographic information, such as elevation, biome, or latitude (as a proxy for seasonality)
- Human related information, such as nearby population density, distance to closest urban or semi-urban area (and maybe a rank of said area by size), or perhaps country level-data such as clean air legislation or a ranking of countries by population's clean air habits (i.e. car driving vs public transport, etc)

Essentially the point of this project was to create a model that attempts to determine what, if any, is the effect of changing vegetation levels on CO2 levels at a local, rather than global level, with other factors held constant and accounting for the change of global CO2 levels due to climate change.

## Data 
### Data Sources
- NASA's MODIS/Terra Vegetation Continuous Fields Yearly Dataset: "a global representation of surface vegetation cover as gradations of three ground cover components: percent tree cover, percent non-tree cover, and percent non-vegetated (bare)."
- NASA's TES Level 2 Dataset: concentration levels for a number of pollutant molecules, including carbon dioxide, carbon monoxide, methane, ozone, ammonia, and more.
### Data Complications
In retrospect, the discrepancies in resolution between the two datasets I used placed some critical limitations on the initial goal of this project.
- Carbon dioxide data: Unfortunately, it turned out that this data is inaccurate at a scale any smaller than 5 degrees latitude by 5 degrees longitude, rendering the initial goal of this project somewhat moot. After averaging and smoothing the data, I was  able to create a usable model, but at a vastly different scale than initially planned.
- Vegetation data: This dataset contains 15 years worth of vegetation data for the entire planet, at the incredibly detailed scale of only 250m, meaning that a single time slice of data contains over 1 billion data points, while the full dataset was over 200GB. While this was good for my initial plans with this project, the limitations of the carbon dioxide data meant that I ended up having to average to a scale of 5 degrees latitude by 5 degrees longitude (>5km square at the equator). As a result, just processing the data necessitated a dask pipeline, let alone using it to build a model.

## Model
After cleaning and merging the datasets, I created both a simple linear regession model and an auto-regressive time series model. While the latter of which was significantly more accurate, this was somewhat inevitable given that global levels of CO2 have been trending upwards for years. With that in mind, any model would be severely limited unless it used the year at a feature. 

In both of these models, the vegetation levels had significantly less impact on the model than other features (altitude and year), but not a completely insignificant impact!

## Main Tools Used
 - Xarray
 - Dask
 - NetCDF
 - Jupyter Notebooks
 - numpy, pandas

## Deliverables:
### Data Pipelines
- [CO2 Data Pipeline](https://github.com/edithalice/vegetation_co2/blob/master/tes_pipeline.py)
- [Vegetation Data Pipeline](https://github.com/edithalice/vegetation_co2/blob/master/veg_pipeline_latitudes.py)
- [Pipeline Implementation](https://github.com/edithalice/vegetation_co2/blob/master/data_cleaning.ipynb) (not properly commented yet yikes)
### Presentation
- [Google Slides](https://docs.google.com/presentation/d/1-zSyVt9V-nQbift6YEkjtDlH88qqMbtuKbNOMJvBraM/edit?usp=sharing)
- [PDF](https://github.com/edithalice/vegetation_co2/blob/master/Vegetation%20Project.pdf)

