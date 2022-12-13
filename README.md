# Annual county-level lead (Pb) vehicle emissions data (1969 to 1981)

This page provides data and code for the paper "Effects of Early-Childhood Exposure to Ambient Lead and Particulate Matter on Adult Personality."

At this stage, only the metadata and R code are publicly available. The full dataset will be posted on this page upon publication of the paper.

## Methodology

![Methodology](methodology.png)

## Data Sources
- Quantity sold in gallons, by month and grade, of leaded and unleaded fuel, by state and for DC, 1970–1982. Ethyl Corporation (Petroleum Chemicals Division) provides these data in its Yearly Report of Gasoline Sales: By States. 
- EPA’s data on nitrogen oxides (NOx) emissions from motor vehicles by state and by Air Quality Control Regions covering all counties in all states, 1972–1980, from the annual National Emissions Report (of EPA’s National Emissions Data System), beginning in 1972. We use NOx emissions to apportion state-wide gasoline use to AQCRs.
- Pb content of gasoline (grams per gallon) by grade (premium, regular, unleaded), for each of 17 marketing districts, which together cover the United States, by season (summer = June, July and August; winter = December, January and February), for 1968–1982. The source is DOE Bartlesville Energy Research Center. 
- Census data on population by county and state, by year. Estimates of county populations by calendar year since 1969 are available [here]( https://www.census.gov/data/datasets/time-series/demo/popest/2010s-counties-detail.html). We use AQCR population to estimate per capita Pb emissions from motor vehicles within an AQCR. 
- We use census data to estimate population density in the core urban area within the county. Data for the 1980s are available at [this link]( https://www2.census.gov/prod2/decennial/documents/1980/1980censusofpopu8011uns_bw.pdf) and data for the 1970s are from [here]( https://www.census.gov/library/publications/1973/dec/population-volume-1.html#par_textimage_43). 
