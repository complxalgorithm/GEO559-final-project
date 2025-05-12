# import necessary libraries
library(tidyverse)
library(tidycensus)
library(tigris)
library(sf)
library(ipumsr)
library(osmdata)
library(terra)
library(fasterize)
library(distanceto)
library(rgee)
library(reticulate)
library(exactextractr)
library(mapview)
library(units)
library(biscale)
library(broom)
library(broom.mixed)
#library(lme4)
#library(gamm4)
library(glmmTMB)
library(DHARMa)
#library(mgcv)
library(spdep)
library(pROC)
library(RColorBrewer)
library(cowplot)
library(kableExtra)
library(gt)

# 
options(tigris_use_cache = TRUE)
census_api_key("1127cdcdf329b894af0f6fe4bad3a4a8412c7dcc")

# sf::sf_use_s2(FALSE)

###################################################################
###################################################################

# get buffalo city boundaries
buff_border <-
  read_csv('https://data.buffalony.gov/api/views/p4ak-r4fg/rows.csv?date=20250308&accessType=DOWNLOAD') %>%
  st_as_sf(wkt = 'Geometry', crs = 'EPSG:4269')

# get buffalo neighborhoods
buff_neighborhoods <- 
  read.csv('https://data.buffalony.gov/resource/ekfg-mtu8.csv') %>%
  st_as_sf(wkt = 'the_geom', crs = 'EPSG:4269') %>%
  rename(
    neighborhood_name = nbhdname,
    area_sqmi = sqmiles,
    geometry = the_geom
  ) %>%
  select(-c(nbhdnum, calcacres, objectid_1))

# import buffalo crimes
### 319,473 ROWS IN RAW FILE
###### 312,867 HAVE COORDINATES
######### 311,222 ARE WITHIN BUFFALO BOUNDARIES
crimes <- 
  read.csv('data/Crime_Incidents_20250303.csv', check.names = FALSE) %>%
  select(`Case Number`, `Incident Datetime`, `Incident Type Primary`, `Parent Incident Type`, 
         `Hour of Day`, `Day of Week`, Address, City, State, zip_code, 
         Latitude, Longitude, neighborhood, `Police District`) %>%
  mutate(
    Month = str_split(
      str_split(.$`Incident Datetime`, ' ', simplify = TRUE)[, 1], '/', simplify = TRUE)[, 1],
    Day = str_split(
      str_split(.$`Incident Datetime`, ' ', simplify = TRUE)[, 1], '/', simplify = TRUE)[, 2],
    Year = str_split(
      str_split(.$`Incident Datetime`, ' ', simplify = TRUE)[, 1], '/', simplify = TRUE)[, 3]
  ) %>%
  mutate_at('Year', as.numeric) %>%
  filter(Year >= 2006 & Year < 2025) %>%
  .[!((.$Latitude == 'UNKNOWN' | .$Longitude == 'UNKNOWN') | (.$Latitude == '' | .$Longitude == '')), ] %>%
  st_as_sf(coords = c('Longitude', 'Latitude'), crs = 'EPSG:4269')

buff_crimes <- crimes[data.frame(st_intersects(crimes, buff_border))$row.id,]

## 88,503 VIOLENT CRIMES
violent_crimes.buff <- buff_crimes %>% 
  filter(`Incident Type Primary` %in% 
           c('MURDER', 'CRIM NEGLIGENT HOMICIDE', 'MANSLAUGHTER', 'RAPE', 'AGGR ASSULT', 
           'AGG ASSULT ON P/OFFICER', 'ROBBERY', 'Robbery', 'ASSAULT', 'Assault', 
           'SEXUAL ABUSE', 'SODOMY', 'Other Sexual Offense'))

## 166,334 PROPERTY CRIMES
property_crimes.buff <- buff_crimes %>%
  filter(`Incident Type Primary` %in%
           c('UUV', 'THEFT OF SERVICES', 'LARCENY/THEFT', 'Theft', 'Theft of Vehicle', 'Breaking & Entering'))

# show bar graph of number of violent crimes per year between 2006 & 2023
ggplot(violent_crimes.buff, aes(Year)) + 
  geom_bar()

# show bar graph of number of property crimes per year between 2006 & 2023
ggplot(property_crimes.buff, aes(Year)) + 
  geom_bar()

# create line graph of homicides between 2006 and 2024 in Buffalo
homicide_types = c('MURDER', 'CRIM NEGLIGENT HOMICIDE', 'MANSLAUGHTER')

violent_crimes.buff %>% 
  st_drop_geometry %>% 
  filter(`Incident Type Primary` %in% homicide_types) %>% 
  group_by(Year) %>% 
  summarize(homicides = n()) %>%
  ungroup() %>%
  ggplot() +
    geom_line(aes(Year, homicides), color = 'blue') +
    geom_smooth(aes(Year, homicides), color = 'red', method = loess, se = FALSE) +
    ggtitle('Homicides in Buffalo, NY (2006 - 2024)') +
    ylab('Homicides') +
    scale_y_continuous(limits = c(25, 85), 
                       breaks = c(30, 40, 50, 60, 70, 80), 
                       labels = c(30, 40, 50, 60, 70, 80)) +
    scale_x_continuous(limits = c(2006, 2024), 
                       breaks = c(2006:2024),
                       labels = c(2006:2024)) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major.y = element_blank())

# import bpd cameras data
### 275 TOTAL CAMERAS
##### 83 CAMERAS WITH VALID INSTALLED DATE
bpd_cameras <- 
  read.csv('https://data.buffalony.gov/api/views/gpj7-v6rr/rows.csv?date=20250304&accessType=DOWNLOAD',
           check.names = FALSE) %>%
  mutate(year = as.numeric(str_split(.$`Installed Date`, '/', simplify = TRUE)[, 3])) %>%
  st_as_sf(coords = c('X', 'Y'), crs = 'EPSG:4269')

bpd_cameras[grepl('2010', bpd_cameras$`Installed Date`), 'Installed Date'] <- 2010
bpd_cameras[grepl('2010', bpd_cameras$`Installed Date`), 'year'] <- 2010

# map of valid and invalid police cameras
bpd_cameras %>% 
  mutate(ifValid = ifelse(is.na(year), 0, 1)) %>% 
  ggplot() + 
    geom_sf(data = buff_border, color = 'black', fill = 'transparent') + 
    geom_sf(aes(color = as.factor(ifValid))) +
    scale_color_manual(name = 'Camera Validity', labels = c('Invalid', 'Valid'), values = c('darkgreen', 'pink')) +
    theme(panel.background = element_rect(fill = '#FFFFFF'), 
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

# separate into dataframes for valid and invalid cameras
valid_cameras <- bpd_cameras %>% filter(!is.na(year))
invalid_cameras <- bpd_cameras %>% filter(is.na(year))

# create 0.25 mile buffer around all valid cameras
#cameras_buffer <- st_buffer(valid_cameras, dist = 402.336)

# import buffalo parks data
### 218 AREAS
##### 
####### 
buff_parks <- 
  read.csv('https://data.buffalony.gov/api/views/tmik-tgt9/rows.csv?date=20250304&accessType=DOWNLOAD',
                       check.names = FALSE) %>%
  st_as_sf(wkt = 'Geometry', crs = 'EPSG:4269') %>%
  filter(
    (`Park Class` %in% c('Major Park', 'Large Park', 'Midsize Park', 'Small Park', 'Parkway')) | (`Park Class` == 'Triangle' & Designated == 1)
  ) %>%
  mutate('Park Category' = 'City Park',
         Area_Sqmi = units::set_units(st_area(.), 'mi2'),
         Area_Acres = units::set_units(Area_Sqmi, 'acre')) %>%
  select('Park name', Address, Year, 'Park Class', 'Park Category', 'Maintained by',
         Neighborhood, 'Council District', Area_Sqmi, Area_Acres) %>%
  rename(Name = 'Park name',
         geometry = Geometry)

buff_parks

# import nys parks, then find Buffalo-based state park and fill in relevant information
nys_parks <- st_read('data/NYS_Park_Polygons.gpkg') %>% st_transform(st_crs(buff_border))

nys_parks_in_buff.idx <- 
  nys_parks %>% 
  st_intersects(., buff_border) %>% data.frame() %>% .$`row.id`

nys_parks_in_buff <- 
  nys_parks[nys_parks_in_buff.idx,] %>% 
  filter(Category == 'State Park') %>%
  select(Category, Label, Original) %>% 
  mutate(
    Address = '1111 Fuhrmann Blvd', 
    'Park Class' = 'Major Park', 
    'Maintained by' = 'NYSOPRHP', 
    Neighborhood = 'Central', 
    'Council District' = 'South', 
    Area_Sqmi = units::set_units(st_area(.), 'mi2'), 
    Area_Acres = units::set_units(Area_Sqmi, 'acre')
  ) %>% 
  rename(
    Name = Label, Year = Original, 'Park Category' = Category, geometry = SHAPE
  ) %>% 
  mutate_at(vars(Year), as.numeric) %>%
  relocate('Park Category', .after = 'Park Class') %>% 
  relocate(Year, .after = Address)

nys_parks_in_buff

# set relavent information for all 7 Buffalo parks under county jurisdiction
erie_park_addresses <- c('Aqua Ln', '3781 Main St', '20 Smith St', '1670 Seneca St', 
                          '11 Fuhrmann Blvd', 'Hertel Ave', '152 Bailey Ave')
erie_park_years <- c(2013, 1902, 2012, 1994, 1987, 1925, 1996)
erie_park_park_classes <- c('Small Park', 'Large Park', 'Small Park', 'Midsize Park',
                            'Large Park', 'Small Park', 'Small Park')
erie_park_neighborhoods <- c('Black Rock', 'University Heights', 'First Ward', 'Seneca-Cazenovia', 
                             'Central', 'Black Rock', 'Seneca-Cazenovia')
erie_park_districts <- c('North', 'University', 'Fillmore', 'Lovejoy', 
                         'South', 'North', 'Lovejoy')

# read in Buffalo-based county parks and add above info
erie_parks_in_buff <- 
  st_read('data/erie_county_parks_in_buffalo.gpkg') %>%
  st_make_valid() %>%
  select(NAME) %>%
  mutate(
    geom = st_simplify(geom, dTolerance = 0.001),
    Address = erie_park_addresses,
    Year = erie_park_years,
    'Park Class' = 'Major Park',
    'Park Category' = 'County Park',
    'Maintained by' = 'ECPD',
    Neighborhood = erie_park_neighborhoods,
    'Council District' = erie_park_districts,
    Area_Sqmi = units::set_units(st_area(.), 'mi2'),
    Area_Acres = units::set_units(Area_Sqmi, 'acre')
  ) %>%
  rename(Name = NAME, geometry = geom)

erie_parks_in_buff

# bind all parks data together to create df of all Buffalo-based city, county, and state parks
# filter out parks that are less than 0.25 acres and all monumental/memorial parks
all_buff_parks <- 
  bind_rows(buff_parks, nys_parks_in_buff, erie_parks_in_buff) %>%
  filter(
    Area_Acres >= units::set_units(0.25, 'acre') & (!((grepl('Monument', Name) | (grepl('Memorial', Name)))))
  )

# add years to certain parks that don't currently have them
all_buff_parks[all_buff_parks$Name == 'Outer Harbor Parkway', 'Year'] <- 2013
all_buff_parks[all_buff_parks$Name == 'McCarthy Park', 'Year'] <- 1973
all_buff_parks[all_buff_parks$Name == 'Remembrance Park', 'Year'] <- 2003
all_buff_parks[all_buff_parks$Name == 'Peter St.', 'Year'] <- 2018
all_buff_parks[all_buff_parks$Name == 'Arlington Park', 'Year'] <- 1866
all_buff_parks[all_buff_parks$Name == 'Bristol Emslie Playground', 'Year'] <- 1998
all_buff_parks[all_buff_parks$Name == 'Genesee Gateway Triangle', 'Year'] <- 2008
all_buff_parks[all_buff_parks$Name == 'Unity Island', 'Year'] <- 2000

all_buff_parks

# map of parks by category
ggplot() +
  geom_sf(data = buff_border, color = 'black') +
  geom_sf(data = all_buff_parks, aes(fill = `Park Category`))


# get Buffalo grocery stores
## 500 TOTAL WERE EXTRACTED
### 102 VALID STORES LOCATED WITHIN BUFFALO BOUNDARIES

# create overpass query
buffalo_stores_query <- 
  opq("Buffalo, New York") %>%
  add_osm_feature(key = 'shop',
                   value = c('supermarket', 'convenience', 'deli', 'general'))

# extract grocery store data as sf
all_stores_raw <- osmdata_sf(buffalo_stores_query)

# combine points and centroids of polygons into stores, filter for certain shop categories,
# then rename rows
stores <-
  bind_rows(
    all_stores_raw$osm_points,
    st_centroid(all_stores_raw$osm_polygons)
  ) %>%
  select(-c(`nysgissam:nysaddresspointid`, phone, website, url, store_ref, `brand:wikidata`,
            any_of(starts_with('contact')), email, image, level, official_name, opening_hours,
            rating, ref, source, check_date)) %>%
  filter(shop %in% c('convenience', 'deli', 'general', 'supermarket')) %>%
  st_transform(crs = st_crs(buff_border))

rownames(stores) <- 1:nrow(stores)

# get indices of stores within Buffalo's boundaries
stores_in_buff.idx <-
  stores %>% 
  st_intersects(., buff_border) %>% data.frame() %>% .$`row.id`

# filter for stores within Buffalo boundaries and store in buff_stores
buff_stores <- stores[stores_in_buff.idx,]

# plot Buffalo stores by shop category
ggplot() + 
  geom_sf(data = buff_border, fill = 'white', color = 'black') + 
  geom_sf(data = buff_stores, aes(color = shop))

########################################################################
### statistical analysis of relationship between violent crime rate and
### presence of police cameras at block group level
########################################################################

#####
# 
#####

# load all nhgis files into lists
nhgis_data_files <- list.files('data/nhgis/vars', pattern = 'zip', full.names = TRUE)
nhgis_sf <- list.files('data/nhgis/spatial', pattern = 'zip', full.names = TRUE)

# process data for each year and store in its own df in nhgis_datasets
nhgis_datasets <- list()
index = 1

for (year in 2010:2012) {
  #
  data_file <- nhgis_data_files[grepl(as.character(year), nhgis_data_files)]
  spatial_file <- nhgis_sf[grepl(as.character(year), nhgis_sf)]
  
  #
  data <- read_nhgis(data_file, verbose = FALSE)
  #codebook <- read_nhgis_codebook(data_file)
  #colnames(data) <- codebook$var_info$var_label
  
  #
  if (year == 2010) {
    data <- 
      data %>% 
      rename(
        tot_pop = 42, 
        white_pop = 43, 
        black_pop = 44, 
        native_american_pop = 45, 
        asian_pop = 46, 
        pacific_islander_pop = 47, 
        other_race_pop = 48, 
        multiracial_pop = 49, 
        non_hispanic_pop = 53, 
        hispanic_pop = 54, 
        tot_pop_over_25 = 55, 
        per_capita_income = 125, 
        ppl_in_poverty = 91, 
        male_no_schooling = 57,
        male_hs_dip_equiv = 65, 
        male_some_college_lt_1yr = 66, 
        male_some_college_gt_1yr_no_deg = 67, 
        male_associates_deg = 68, 
        male_bachelors_deg = 69, 
        male_masters_deg = 70, 
        male_professional_deg = 71, 
        male_doctoral_deg = 72, 
        female_no_schooling = 74,
        female_hs_dip_equiv = 82, 
        female_some_college_lt_1yr = 83, 
        female_some_college_gt_1yr_no_deg = 84, 
        female_associates_deg = 85, 
        female_bachelors_deg = 86, 
        female_masters_deg = 87, 
        female_professional_deg = 88, 
        female_doctoral_deg = 89
      ) %>% 
      select(GEOID, YEAR, COUNTY, 42:49, 53:55, 125, 91, 57, 65:72, 74, 82:89) %>%
      mutate(
        no_schooling = male_no_schooling + female_no_schooling,
        high_school_dip_equiv = male_hs_dip_equiv + female_hs_dip_equiv,
        some_college_lt_1yr = male_some_college_lt_1yr + female_some_college_lt_1yr,
        some_college_gt_1yr_no_deg = male_some_college_gt_1yr_no_deg + female_some_college_gt_1yr_no_deg,
        associates_deg = male_associates_deg + female_associates_deg,
        bachelors_deg = male_bachelors_deg + female_bachelors_deg,
        masters_deg = male_masters_deg + female_masters_deg,
        professional_deg = male_professional_deg + female_professional_deg,
        doctoral_deg = male_doctoral_deg + female_doctoral_deg
      ) %>%
      select(-c(starts_with('male'), starts_with('female')))
  } else if (year == 2011) {
    data <- 
      data %>% 
      rename(
        tot_pop = 43, 
        white_pop = 44, 
        black_pop = 45, 
        native_american_pop = 46, 
        asian_pop = 47, 
        pacific_islander_pop = 48, 
        other_race_pop = 49, 
        multiracial_pop = 50, 
        non_hispanic_pop = 54, 
        hispanic_pop = 55, 
        tot_pop_over_25 = 56, 
        per_capita_income = 126, 
        ppl_in_poverty = 92, 
        male_no_schooling = 58,
        male_hs_dip_equiv = 66, 
        male_some_college_lt_1yr = 67, 
        male_some_college_gt_1yr_no_deg = 68, 
        male_associates_deg = 69, 
        male_bachelors_deg = 70, 
        male_masters_deg = 71, 
        male_professional_deg = 72, 
        male_doctoral_deg = 73, 
        female_no_schooling = 75,
        female_hs_dip_equiv = 83, 
        female_some_college_lt_1yr = 84, 
        female_some_college_gt_1yr_no_deg = 85, 
        female_associates_deg = 86, 
        female_bachelors_deg = 87, 
        female_masters_deg = 88, 
        female_professional_deg = 89, 
        female_doctoral_deg = 90
      ) %>% 
      select(GEOID, YEAR, COUNTY, 43:50, 54:56, 126, 92, 58, 66:73, 75, 83:90) %>%
      mutate(
        no_schooling = male_no_schooling + female_no_schooling,
        high_school_dip_equiv = male_hs_dip_equiv + female_hs_dip_equiv,
        some_college_lt_1yr = male_some_college_lt_1yr + female_some_college_lt_1yr,
        some_college_gt_1yr_no_deg = male_some_college_gt_1yr_no_deg + female_some_college_gt_1yr_no_deg,
        associates_deg = male_associates_deg + female_associates_deg,
        bachelors_deg = male_bachelors_deg + female_bachelors_deg,
        masters_deg = male_masters_deg + female_masters_deg,
        professional_deg = male_professional_deg + female_professional_deg,
        doctoral_deg = male_doctoral_deg + female_doctoral_deg
      ) %>%
      select(-c(starts_with('male'), starts_with('female')))
  } else {
    data <- 
      data %>% 
      rename(
        tot_pop = 43, 
        white_pop = 44, 
        black_pop = 45, 
        native_american_pop = 46, 
        asian_pop = 47, 
        pacific_islander_pop = 48, 
        other_race_pop = 49, 
        multiracial_pop = 50, 
        non_hispanic_pop = 54, 
        hispanic_pop = 55, 
        tot_pop_over_25 = 56, 
        per_capita_income = 116, 
        ppl_in_poverty = 82, 
        no_schooling = 57,
        high_school_diploma = 72,
        ged = 73,
        some_college_lt_1yr = 74,
        some_college_gt_1yr_no_deg = 75,
        associates_deg = 76,
        bachelors_deg = 77,
        masters_deg = 78,
        professional_deg = 79,
        doctoral_deg = 80
      ) %>% 
      select(GEOID, YEAR, COUNTY, 43:50, 54:56, 116, 82, 57, 72:80) %>%
      mutate(high_school_dip_equiv = high_school_diploma + ged) %>%
      select(-c(high_school_diploma, ged))
  }
  
  # 
  data <- 
    data %>%
    filter(COUNTY == 'Erie County') %>%
    relocate(GEOID, .before = YEAR) %>%
    mutate(GEOID = str_split(GEOID, '15000US', simplify=T)[,2],
           YEAR = str_split(YEAR, '-', simplify = TRUE)[,2]) %>%
    select(-COUNTY)
  
  # 
  if (year == 2010) {
    sf <- read_ipums_sf(spatial_file, verbose=F) %>% 
      filter(COUNTYFP10 == '029') %>% 
      select(GEOID10) %>%
      rename(GEOID = GEOID10) %>%
      st_transform(crs = st_crs(buff_border))
    
    d <- inner_join(sf, data, by = 'GEOID') %>% 
      relocate(geometry, .after = last_col())
  } else {
    sf <- read_ipums_sf(spatial_file, verbose=F) %>% 
      filter(COUNTYFP == '029') %>% 
      select(GEOID) %>%
      st_transform(crs = st_crs(buff_border))
    
    d <- inner_join(sf, data, by = 'GEOID') %>% 
      relocate(geometry, .after = last_col())
  }
  
  # 
  nhgis_datasets[[index]] <- d
  
  index = index + 1
}

# 
erie_blocks.nhgis <-
  bind_rows(nhgis_datasets) %>%
  rename(year = YEAR) %>%
  mutate_at(vars(per_capita_income, year), as.numeric) %>%
  mutate(
    area_sqmi = as.numeric(st_area(.)) * 3.861E-7,
    tot_pop_density = tot_pop / area_sqmi,
    pop_25plus_density = tot_pop_over_25 / area_sqmi,
    total_income = per_capita_income * tot_pop,
    some_college_no_deg = some_college_lt_1yr + some_college_gt_1yr_no_deg,
    bachelors_deg_at_least = bachelors_deg + masters_deg + professional_deg + doctoral_deg,
    less_than_bachelors_deg = tot_pop_over_25 - bachelors_deg_at_least,
    less_than_high_school_grad = tot_pop_over_25 - (high_school_dip_equiv + some_college_no_deg + bachelors_deg_at_least)
  ) %>%
  relocate(geometry, .after = last_col())

# get ACS block group data for pre-covid and post-covid years and store in separate dfs
erie_blocks.pre_covid <- 
  map(2013:2019, get_acs_block_group_data) %>%
  bind_rows()

erie_blocks.post_covid <- 
  map(2020:2023, get_acs_block_group_data) %>%
  bind_rows()

# separate buffalo block groups into their own dataframes
block_buff_intersects.nhgis <- data.frame(st_intersects(st_centroid(erie_blocks.nhgis), buff_border))$row.id
missing_blocks.nhgis <- which(erie_blocks.nhgis$GEOID %in% c(360290011003))
buff_blocks.nhgis <-
  erie_blocks.nhgis[c(block_buff_intersects.nhgis, missing_blocks.nhgis),] %>%
  st_make_valid()

block_buff_intersects.pre_covid <- data.frame(st_intersects(st_centroid(erie_blocks.pre_covid), buff_border))$row.id
missing_blocks.pre_covid <- which(erie_blocks.pre_covid$GEOID %in% c(360290011003, 360290072021, 360290058022))
buff_blocks.pre_covid <- 
  erie_blocks.pre_covid[c(block_buff_intersects.pre_covid, missing_blocks.pre_covid),] %>%
  st_make_valid()

block_buff_intersects.post_covid <- data.frame(st_intersects(st_centroid(erie_blocks.post_covid), buff_border))$row.id
missing_blocks.post_covid <- which(erie_blocks.post_covid$GEOID %in% c(360290011003, 360290072021, 360290058022))
buff_blocks.post_covid <- 
  erie_blocks.post_covid[c(block_buff_intersects.post_covid, missing_blocks.post_covid),] %>%
  st_make_valid()

# combine nhgis (2010-2012) and pre-covid (2013-2019) dataframes into
# full Buffalo pre-covid dataframe, and remove empty rows
buff_blocks.pre_covid <- 
  bind_rows(buff_blocks.nhgis, buff_blocks.pre_covid) %>%
  filter(!is.na(year))

# list of all variables pulled
all_vars <- c(
  'tot_pop', 'white_pop', 'black_pop', 'native_american_pop', 'asian_pop', 'pacific_islander_pop',
  'other_race_pop', 'multiracial_pop', 'non_hispanic_pop', 'hispanic_pop', 'ppl_in_poverty',
  'tot_pop_over_25', 'no_schooling', 'some_college_lt_1yr', 'some_college_gt_1yr_no_deg',
  'associates_deg', 'bachelors_deg', 'masters_deg', 'professional_deg', 'doctoral_deg',
  'high_school_dip_equiv', 'some_college_no_deg', 'bachelors_deg_at_least',
  'less_than_bachelors_deg', 'less_than_high_school_grad', 'per_capita_income'
)

# get erie county block boundaries in 2020, which was when block group boundaries changed
# then get blocks in buffalo boundary
# used to estimate population in each segment of new boundaries
erie_blocks.2020 <- tigris::blocks(state = 'NY', county = 'Erie', year = 2020)
buff_blocks.2020 <- st_intersection(erie_blocks.2020, buff_border, model = 'open')

# block group boundaries in 2020
buff_block_groups.2020 <- buff_blocks.post_covid %>% filter(year == 2020)

# list of pre-covid years
pre_covid_yrs <- unique(buff_blocks.pre_covid$year)

# empty list of interpolated block group boundary dataframes for each pre-covid year
# will bind rows after interpolation
pre_covid_data <- list()

# estimate data for years between 2013-2019 using block group boundaries in 2020-2023
# applies dasymetric population-weighted areal interpolation technique (interpolation_pw)
index = 1
for (i in pre_covid_yrs) {
  # 
  yr_df <- buff_blocks.pre_covid[buff_blocks.pre_covid$year == i,] %>%
    select(-c(year, per_capita_income, tot_pop_density, pop_25plus_density, area_sqmi))
  
  # apply population-weighted areal interpolation technique
  # to conform boundaries and estimate the data
  interpolated <- 
    yr_df %>%
    interpolate_pw(
      ., buff_block_groups.2020, 
      to_id = 'GEOID', 
      extensive = TRUE, 
      weights = buff_blocks.2020, 
      weight_column = 'POP20', 
      crs = st_crs(buff_blocks.pre_covid)
    ) %>%
    mutate(
      year = as.numeric(i),
      per_capita_income = total_income / tot_pop,
      area_sqmi = as.numeric(st_area(.) * 3.861E-7),
      tot_pop_density = tot_pop / area_sqmi,
      pop_25plus_density = tot_pop_over_25 / area_sqmi
    )
  
  #
  pre_covid_data[[index]] <- interpolated
  
  index = index + 1
}

# combine interpolated 2010-2019 (which now follow 2020 block group boundaries)
# with post-COVID data, also round interpolated data and remove a couple columns
buff_blocks <-
  bind_rows(pre_covid_data) %>%
  mutate_at(vars(all_vars), ~ round(.)) %>%
  bind_rows(., buff_blocks.post_covid) %>%
  select(-total_income) %>%
  mutate(
    pct_below_poverty = ppl_in_poverty / tot_pop,
    pct_bachelors_deg_at_least = bachelors_deg_at_least / tot_pop_over_25,
    pct_less_than_high_school_grad = less_than_high_school_grad / tot_pop_over_25,
    pct_white_pop = white_pop / tot_pop,
    pct_black_pop = black_pop / tot_pop,
    pct_hispanic_pop = hispanic_pop / tot_pop
    #has_crime = NA_integer_,
    #has_camera = NA_integer_,
    #First_Camera_Year = NA_integer_,
    #Num_Cameras = NA_integer_,
    #Cameras_Per_1000 = NA_real_,
    #Post_1st_Cam_Install = NA_integer_,
    #time_to_install_first_cam = NA_integer_,
    #Camera_Group = NA_integer_,
    #num_violent_crimes = NA_integer_,
    #violent_crime_rate = NA_real_,
    #pct_park_proximity = NA_real_
  ) %>%
  relocate(year, .after = GEOID) %>%
  relocate(geometry, .after = last_col()) %>%
  relocate(area_sqmi, .before = geometry)

# 
interpolated_pre_covid_pop <- 
  bind_rows(pre_covid_data) %>%
  st_drop_geometry() %>%
  summarize(pop = sum(.$tot_pop, na.rm=T)) %>%
  .[, 'pop']

original_pre_covid_pop <-
  buff_blocks.pre_covid %>%
  st_drop_geometry() %>%
  summarize(pop = sum(.$tot_pop, na.rm=T)) %>%
  .[, 'pop']

cat(
  paste(
    paste('Original Pop', original_pre_covid_pop, sep=': '),
    paste('Interpolated Pop', round(interpolated_pre_covid_pop), sep=': '),
    paste('Difference: ', round(((interpolated_pre_covid_pop - original_pre_covid_pop) / original_pre_covid_pop) * 100, 3), '%', sep=''),
    sep = '\n'
  )
)

# initialize separate lists to hold avg_park_distance rasters and buffer rasters
avg_park_distance_rasters <- c()
buffer_rasters <- c()

# determine number of cameras in block group in year,
# as well as number of violent crimes and violent crime rate
# also process data in general in preparation for logistic regression analysis
buff_blocks <- 
  lapply(2010:2023, function(yr) {
    ##print(as.character(yr))
    
    # UTM zone 18 PCS
    proj_crs <- 'EPSG:32618'
    
    #
    yr_buff_blocks <- 
      buff_blocks %>% 
      filter(year == yr) %>% 
      st_transform(crs = proj_crs) %>%
      mutate(area_sqmi = as.numeric(units::set_units(st_area(.), 'mi2'))) %>%
      st_make_valid()
    
    #
    cams <- 
      valid_cameras[valid_cameras$year <= yr,] %>% 
      st_transform(crs = proj_crs) %>%
      st_make_valid()
    
    #
    if (year >= 2018) {
      cams <- bind_rows(cams, st_transform(bpd_cameras[is.na(bpd_cameras$year),], 
                                           crs = proj_crs))
      cams[is.na(cams)] <- 2018
    } else {
      cams <- cams
    }
    
    #
    crimes <- 
      violent_crimes.buff[violent_crimes.buff$Year == yr,] %>% 
      st_transform(crs = proj_crs) %>%
      st_make_valid()
    
    #
    cams_per_block <-
      st_join(
        yr_buff_blocks %>% select(GEOID),
        st_join(cams, yr_buff_blocks, join = st_within) %>%
          group_by(GEOID) %>%
          summarize(Num_Cameras = n(), .groups = 'drop'),
        suffix = c('', '_y')
      ) %>%
      select(-GEOID_y)
    
    #
    crimes_per_block <-
      st_join(
        yr_buff_blocks %>% select(GEOID),
        st_join(crimes, yr_buff_blocks, join = st_within) %>%
          group_by(GEOID) %>%
          summarize(num_violent_crimes = n(), .groups = 'drop'),
        suffix = c('', '_y')
      ) %>%
      select(-GEOID_y)
    
    #
    parks_yr <- 
      all_buff_parks[all_buff_parks$Year <= yr,] %>% 
      st_transform(crs = proj_crs) %>%
      st_make_valid()
    
    # dissolved buffer of parks
    parks_buffer <- 
      st_buffer(parks_yr, dist = 402.336) %>%
      st_union() %>%
      st_make_valid() %>%
      st_transform(crs = proj_crs)
    #st_collection_extract('POLYGON')
    
    # 
    ext <- st_bbox(yr_buff_blocks)
    
    #
    parks_temp_raster <- rast(
      extent = ext,
      resolution = 10,
      crs = proj_crs
    )
    
    # rasterize buffer (cells inside buffer = 1, outside = 0)
    buffer_raster <- rasterize(vect(parks_buffer), parks_temp_raster, 
                               field = 1, background = 0)
    
    # 
    buffer_rasters[[as.character(yr)]] <- buffer_raster
    
    # 
    pct_park_coverage <- exactextractr::exact_extract(buffer_raster, yr_buff_blocks, 'mean')
    
    #
    dist_raster <- distanceto::distance_raster(
      y = parks_yr %>% st_collection_extract('POLYGON'),
      cellsize = 30, 
      extent = ext, 
      check = FALSE
    )
    
    #
    avg_park_distances <- exactextractr::exact_extract(dist_raster, yr_buff_blocks, 'mean')
    
    avg_park_distance_rasters[[as.character(yr)]] <- avg_park_distances
    
    # 
    return(
      yr_buff_blocks %>%
        st_join(cams_per_block, st_nearest_feature, suffix = c('', '_y')) %>%
        select(-GEOID_y) %>%
        st_join(crimes_per_block, st_nearest_feature, suffix = c('', '_y')) %>%
        select(-GEOID_y) %>%
        mutate(
          Num_Cameras = replace_na(Num_Cameras, 0),
          num_violent_crimes = replace_na(num_violent_crimes, 0),
          pct_park_proximity = pct_park_coverage,
          avg_distance_to_park = avg_park_distances,
          violent_crime_rate = (num_violent_crimes / tot_pop) * 1000,
          has_crime = ifelse(num_violent_crimes > 0, 1, 0),
          has_camera = ifelse(Num_Cameras > 0, 1, 0),
          Cameras_Per_1000 = (Num_Cameras / tot_pop) * 1000
        ) %>%
        relocate(geometry, .after = last_col()) %>%
        st_transform(crs = 'EPSG:4269')
    )
  }) %>%
  bind_rows()

##########
# average June surface temperature for each block group
##########

# initialize GEE
py_interp = 'C:/Users/scsan/anaconda3/envs/geo559/python.exe'
conda_env = 'C:\\Users\\scsan\\anaconda3\\envs\\geo559'
google_user = 'scsandman6496@gmail.com'
GEE_API_KEY = 'AIzaSyAZ8l2z3JnchDwW1IKhxlmc43-qo9F_4wg'
OAUTH2_SECRET = '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'

Sys.setenv(HOME = 'C:/Users/scsan/Documents')

ee_install_set_pyenv(py_path = py_interp, py_env = 'geo559')
#ee_install()
#ee_install_upgrade()
ee_Authenticate(user = google_user, authorization_code = OAUTH2_SECRET)
ee_Initialize(user = google_user, cloud_api_key = GEE_API_KEY, 
              auth_mode = 'gcloud', credentials = OAUTH2_SECRET)

# 
month = 6
cloud_threshold = 15

###
# define get_lst_raster() function: get median June surface temp in area for
# each year between 2010-2023
###
get_lst_raster <- function(year, aoi) {
  #
  start_date <- sprintf('%d-%02d-01', year, month)
  end_date <- sprintf('%d-%02d-30', year, month)
  
  #
  if (year < 2013) {
    # 
    collection <- ee$ImageCollection('LANDSAT/LE07/C02/T1_L2')
    band <- 'ST_B6'
  } else {
    #
    collection <- ee$ImageCollection('LANDSAT/LC08/C02/T1_L2')
    band <- 'ST_B10'
  }
  
  # 
  img <-
    collection$
    filterBounds(aoi)$
    filterDate(start_date, end_date)$
    filter(ee$Filter$lt('CLOUD_COVER', cloud_threshold))$
    select(band)$
    median()$
    subtract(273.15)$
    multiply(9/5)$
    add(32.0)$
    rename('LST_C')$
    clip(aoi)
  
  #
  return(img)
}

###
# define get_mean_lst_by_bg() function - load year's raster then compute zonal mean
# surface temperature for each year
###
get_mean_lst_by_bg <- function(year) {
  cat('Processing Year: ', year, '\n')
  
  # get current year's block groups
  buff_bg_yr <- buff_blocks %>% filter(year == year)
  
  # convert year's block groups to EE and unified bounding box for image clipping
  buff_union <- st_union(buff_bg_yr)
  buff_ee <- sf_as_ee(buff_union)
  
  # 
  img <- get_lst_raster(year, buffalo_ee)
  img_rast <- ee_as_raster(image = img, region = buffalo_ee$geometry(), scale = 30)
  r <- rast(img_rast[[1]])
  names(r) <- 'LST_C'
  
  # 
  stats_df <- exact_extract(r, buff_bg_yr, 'mean') %>%
    bind_cols(buff_blocks) %>%
    select(GEOID, mean) %>%
    mutate(year = year)
  
  # 
  return(stats_df)
}

# get surface temperature for each block group for each year
map_dfr(2022:2023, get_mean_lst_by_bg) %>%
  rename(avg_lst_fahrenheit = mean)

#####################################

# plot per capita income by year
ggplot(buff_blocks) + 
  geom_sf(aes(fill = per_capita_income)) + 
  scale_fill_viridis_c(option = 'inferno') +
  facet_wrap(~year)

# plot pct below poverty by year
ggplot(buff_blocks) +
  geom_sf(aes(fill = pct_below_poverty)) +
  scale_fill_viridis_c(option = 'inferno') +
  facet_wrap(~year)

# plot pct less than high school grad by year
ggplot(buff_blocks) +
  geom_sf(aes(fill = pct_less_than_high_school_grad)) +
  scale_fill_viridis_c(option = 'inferno') +
  facet_wrap(~year)

# plot pct bachelors degree at least by year
ggplot(buff_blocks) +
  geom_sf(aes(fill = pct_bachelors_deg_at_least)) +
  scale_fill_viridis_c(option = 'inferno') +
  facet_wrap(~year)

#####
# bivariate choropleth map (violent crime rate and pct poverty)
#####

# determine classes
bivariate <- bi_class(buff_blocks %>% filter(year == 2023), x = violent_crime_rate, y = pct_below_poverty, 
                      style = 'quantile', dim = 3)

# create legend
bi_legend <- bi_legend(pal = 'GrPink2',
                       dim = 3,
                       xlab = 'Violent Crime Rate',
                       ylab = '% in Poverty',
                       size = 5)

# create map
bi_map <-
  ggplot() +
    geom_sf(data = bivariate, 
            aes(fill = bi_class), 
            color = 'lightgray', 
            show.legend = FALSE) +
    bi_scale_fill(pal = 'GrPink2', dim = 3) +
    bi_theme()

# combine map with legend
final_bi_map <-
  ggdraw() +
    draw_plot(bi_map, 0, 0, 1, 1) +
    draw_plot(bi_legend, 0.15, 0.65, 0.2, 0.2)

final_bi_map

#############################################################################

######################
# Regression Analysis
######################

# scale quantitative independent variables, then place each observation in a year group
buff_blocks.scaled <- 
  buff_blocks %>%
  mutate(
    per_capita_income = scale(per_capita_income)[,1],
    pct_below_poverty = scale(pct_below_poverty)[,1],
    pct_bachelors_deg_at_least = scale(pct_bachelors_deg_at_least)[,1],
    pct_black_pop = scale(pct_black_pop)[,1],
    pct_park_proximity = scale(pct_park_proximity)[,1],
    avg_distance_to_park = scale(avg_distance_to_park)[,1],
    year_group = case_when(
      year <= 2012 ~ '2010-2012',
      year <= 2015 ~ '2013-2015',
      year <= 2018 ~ '2016-2019',
      TRUE ~ '2020-2023'
    )
  )

# logistic regression with block group as random effect, 
# year_group as fixed effect, and scaled data
tmb_model <- glmmTMB(
  has_crime ~ has_camera + per_capita_income + pct_below_poverty + pct_bachelors_deg_at_least + 
    pct_black_pop + pct_park_proximity + factor(year_group) + (1 | GEOID),
  data = buff_blocks.scaled,
  family = binomial()
)

summary(tmb_model)

###
# model diagnostics using DHARMa
###

# simulate residuals to help assess model fit,
# then show diagnostic plots
tmb_model.simres <- simulateResiduals(tmb_model)
plot(tmb_model.simres)

###
# diagnose assumption of normally distributed estimated random effects
###

# get random effects, then extract random intercepts
ranef.d <- ranef(tmb_model)$cond$GEOID
random_intercepts <- ranef.d[['(Intercept)']]

# qq plot
qqnorm(random_intercepts)
#qqnorm(random_intercepts, col = 'red')

###
# Spatial Autocorrelation (Moran's I test)
###


### NEED TO CREATE buff_blocks.ranef BEFORE THIS STEP (DO BEFORE DIAGNOSTIC TESTING)

# centroids of each block group, then create neighbors based on distance
# (e.g., 10 nearest neighbors)
bg_centroids <- st_centroid(buff_blocks.ranef)
coords <- st_coordinates(bg_centroids)
neighbors <- knearneigh(coords, k = 10) %>% knn2nb()
weights <- nb2listw(neighbors, style = 'W')

# Moran's I test
moran_test <- moran.test(buff_blocks.ranef$random_intercept, weights)

moran_test

###
# accuracy evaluation
###

# get predicted probabilities
pred_probs <- predict(tmb_model, type = 'response')

# convert probabilities to binary predictions
pred_class <- ifelse(pred_probs > 0.5, 1, 0)

# create confusion matrix, then extract values
conf_matrix <- table(Predicted = pred_class, Actual = buff_blocks.scaled$has_crime)

TN <- conf_matrix[1,1]
FP <- conf_matrix[2,1]
FN <- conf_matrix[1,2]
TP <- conf_matrix[2,2]

# create confusion matrix dataframe
confusion_df <- tibble(
  Outcome = c('True Negatives', 'False Positives', 
              'False Negatives', 'True Positives'),
  Count = c(TN, FP, FN, TP)
)

# create gt table for confusion matrix
confusion_gt <-
  confusion_df %>%
  gt() %>%
  tab_header(title = 'Confusion Matrix: Crime Prediction Model') %>%
  fmt_number(columns = Count, decimals = 0)

# calculate metrics
accuracy <- (TP + TN) / sum(conf_matrix)
sensitivity <- TP / (TP + FN)
specificity <- TN / (TN + FP)

# create metrics dataframe
metrics_df <- tibble(
  Metric = c('Accuracy', 'Sensitivity (Recall)', 'Specificity'),
  Value = c(accuracy, sensitivity, specificity)
)

# create gt table for metrics
metrics_gt <- 
  metrics_df %>%
  mutate(Value = round(Value, 4)) %>%
  gt() %>%
  tab_header(title = 'Model Performance Metrics') %>%
  fmt_percent(columns = Value, decimals = 2)

# show gt tables of confusion matrix and metrics
confusion_gt
metrics_gt

###
# ROC Curve & AUC (evaluating binary classifiers)
###
# roc_obj <- roc(buff_blocks.scaled$has_crime, pred_probs)
# plot(roc_obj)
# auc(roc_obj)

###
# odds ratio interpretation
###

# extract fixed effects
# then exponentiate to get odds ratios
fixef_vals <- fixef(tmb_model)$cond
odds_ratios <- exp(fixef_vals)

# get standard errors
se_vals <- summary(tmb_model)$coefficients$cond[, 'Std. Error']

# compute 95% confidence interval
lower_ci <- exp(fixef_vals - 1.96 * se_vals)
upper_ci <- exp(fixef_vals + 1.96 * se_vals)

# get fixed effects data
fixed_effects <-
  data.frame(
    Estimate = fixef_vals,
    Odds_Ratio = odds_ratios,
    CI_Lower = lower_ci,
    CI_Upper = upper_ci,
    p_value = summary(tmb_model)$coefficients$cond[, 'Pr(>|z|)']
  )

# odds ratio table
odds_table <-
  fixed_effects %>%
  mutate(
    Predictor = rownames(.),
    Odds_Ratio = round(Odds_Ratio, 2),
    CI = paste0('[', round(CI_Lower, 2), ', ', round(CI_Upper, 2), ']'),
    p_value = formatC(p_value, format = 'e', digits = 2)
  ) %>%
  select(Predictor, Odds_Ratio, CI, p_value) %>%
  gt() %>%
  tab_header(
    title = 'Odds Ratios from Logit Mixed-Effects Model',
    subtitle = '95% Confidence Intervals and p-values'
  ) %>%
  fmt_markdown(columns = c(CI)) %>%
  cols_label(
    Predictor = 'Predictor',
    Odds_Ratio = 'Odds Ratio',
    CI = '95% CI',
    p_value = 'p-value'
  ) %>%
  tab_options(
    table.font.size = 'small',
    heading.align = 'left'
  )

odds_table

# add signifance to highlight predictors that are significant
fixed_effects_plot <-
  fixed_effects %>%
  mutate(
    Predictor = rownames(.),
    Significance = ifelse(p_value < 0.05, 'Significant', 'Not Significant')
  )

# plot fixed effects with confidence intervals and significance
ggplot(fixed_effects_plot, aes(x = reorder(Predictor, Odds_Ratio), y = Odds_Ratio)) +
  geom_point(aes(color = Significance), size = 3) +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper, color = Significance), width = 0.2) +
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'red') +
  scale_y_log10() + # log scale for better visualization
  coord_flip() +
  labs(
    title = 'Odds Ratios with 95% Confidence Intervals',
    subtitle = 'Logit Mixed-Effects Model (glmmTMB)',
    x = 'Predictor',
    y = 'Odds Ratio (log scale)',
    color = 'Significance'
  ) +
  theme_minimal(base_size = 12)

###
# random effects interpretation
###

# get the ordered unique GEOIDs
geoid_levels <- unique(tmb_model$frame$GEOID)

# extract random effects (each unique GEOID and its random intercept)
ranef_data <- 
  ranef(tmb_model)$cond$GEOID %>%
  as_tibble() %>%
  mutate(GEOID = geoid_levels) %>%
  rename(random_intercept = `(Intercept)`) %>%
  relocate(GEOID, .before = random_intercept)

# table of a sample of the random effects (GEOIDs and their random intercepts)
ranef_data %>%
  head(10) %>%
  kbl(
    digits = 3,
    caption = 'Sample of Random Intercepts by GEOID',
    col.names = c('GEOID', 'Random Intercept (log-odds)')
  ) %>%
  kable_styling(full_width = FALSE, bootstrap_options = c('striped', 'hover'))

# join random effects to block groups
buff_blocks.ranef <-
  buff_blocks %>%
  distinct(GEOID, .keep_all = TRUE) %>%
  left_join(ranef_data, by = 'GEOID')

# make map of random effects to see which areas have a higher risk of crime
# than expected possibly due to unmeasured factors, and which areas have lower
# than expected crime risk
ranef_map <-
  ggplot(buff_blocks.ranef) +
  geom_sf(aes(fill = random_intercept), color = NA) +
  scale_fill_distiller(
    palette = 'RdYlBu',
    name = 'Random Intercept\n(Log-Odds)',
    na.value = 'grey80',
    direction = -1
  ) +
  labs(
    title = 'Random Effects by Block Group',
    subtitle = 'Baseline Crime Risk After Controlling for Predictors'
  ) +
  theme_map()

ranef_map

###
# investigating persistent positive random effects by neighborhood
###

# see where block groups and neighborhoods intersect
# calculate area of each intersected piece
bg_neighborhoods_overlap <- 
  st_intersection(st_transform(buff_blocks.ranef, crs = 'EPSG:26918'), 
                  st_transform(buff_neighborhoods, crs = 'EPSG:26918')) %>%
  mutate(overlap_area = st_area(.)) %>%
  st_transform(crs = 'EPSG:4269')

# for each block group, keep the neighborhood where the overlap is largest
bg_neighborhoods <- 
  bg_neighborhoods_overlap %>%
  group_by(GEOID) %>%
  slice_max(overlap_area) %>%
  ungroup() %>%
  .[,-c(2:45)] %>%
  rename(area_sqmi = `area_sqmi.1`)

# compute average and maximum random intercept and number of block groups
# for each neighborhood
ranef_summary <-
  bg_neighborhoods %>%
  st_drop_geometry() %>%
  group_by(neighborhood_name) %>%
  summarize(
    avg_random_effect = mean(random_intercept, na.rm = TRUE),
    max_random_effect = max(random_intercept, na.rm = TRUE),
    num_blockgroups = n()
  )

# 
neighborhood_ranef <- left_join(buff_neighborhoods, 
                                ranef_summary, 
                                by = 'neighborhood_name')

# show table of avg random effects in descending order
neighborhood_ranef %>%
  st_drop_geometry %>%
  arrange(desc(avg_random_effect)) %>%
  gt()

# identify top 5 most at-risk neighborhoods, then get their centroids
# for labeling of top 5 most at-risk neighborhoods on the map
top5_neighborhoods <-
  neighborhood_ranef %>%
  arrange(desc(avg_random_effect)) %>%
  slice_head(n = 5)

top5_centroids <- st_centroid(st_transform(top5_neighborhoods, crs = 'EPSG:32618'))

# map neighborhood-level random effects
neighborhood_ranef_map <-
  ggplot(neighborhood_ranef) +
    geom_sf(aes(fill = avg_random_effect), color = 'white') +
    scale_fill_distiller(
      palette = 'RdYlBu',
      name = 'Random Intercept\n(Log-Odds)',
      na.value = 'grey80',
      direction = -1,
      limits = c(
        min(neighborhood_ranef$avg_random_effect),
        max(neighborhood_ranef$avg_random_effect)
      )
    ) +
  geom_sf_text(
    data = top5_centroids,
    aes(label = neighborhood_name),
    size = 2, color = 'black', fontface = 'bold'
  ) +
    labs(
      title = 'Average Unexplained Crime Risk by Neighborhood',
      subtitle = 'Based on Random Intercepts from Mixed-Effects Model'
    ) +
    theme_map()

neighborhood_ranef_map

#########################################################

# # event study approach
# event_study <- 
#   feols(
#     violent_crime_rate ~ i(time_to_install_first_cam, Camera_Group, ref = -1) + Num_Cameras + black_pop + pct_park_proximity | GEOID + year, 
#     data = buff_blocks %>% filter(First_Camera_Year != 0)
#   )
# 
# summary(event_study)
# broom::tidy(event_study)

# # difference-in-difference regression analysis
# did <- 
#   feols(
#     violent_crime_rate ~ Post_1st_Cam_Install:Camera_Group + Cameras_Per_1000 + black_pop + pct_below_poverty + pct_less_than_high_school_grad + pct_park_proximity | GEOID + year, 
#     data = buff_blocks %>% filter(First_Camera_Year != 0)
#   )
# 
# summary(did)
# broom::tidy(did)

##########
# SOURCES
##########
# IPUMS from NHGIS (2010-2012):
#   Steven Manson, Jonathan Schroeder, David Van Riper, Katherine Knowles, Tracy Kugler, Finn Roberts, and Steven Ruggles. IPUMS National Historical Geographic Information System: Version 19.0 [dataset]. Minneapolis, MN: IPUMS. 2024. http://doi.org/10.18128/D050.V19.0
# Dasymetric Mapping (Population-Weighted Areal Interpolation):
#   https://walker-data.com/census-r/spatial-analysis-with-us-census-data.html#population-weighted-areal-interpolation
#   https://compass.onlinelibrary.wiley.com/doi/pdf/10.1111/j.1749-8198.2009.00220.x?casa_token=xnpg8VtB2gsAAAAA%3AliXtuQU_NpOQfOsNo75ejIY-33YISTCP9j-2SO5fqGVssbgt620u_cC7u_oR-ek1F629nHNsUK4vIDYH
#
# Packages
#
# doParallel & foreach: 
#   https://cran.r-project.org/web/packages/doParallel/vignettes/gettingstartedParallel.pdf
# tidycensus: 
#   Walker K, Herman M (2025). tidycensus: Load US Census Boundary and Attribute Data as 'tidyverse' and 'sf'-Ready Data Frames. R package version 1.7.1, https://walker-data.com/tidycensus/. 
