###
# define get_block_group_data: retrieves census block group ACS data for Erie County tracts for a given year
###
get_acs_block_group_data <- function(year) {
  return(
    get_acs(
      geography = 'block group', 
      variables = c(
        tot_pop = 'B01003_001',
        white_pop = 'B02001_002',
        black_pop = 'B02001_003',
        native_american_pop = 'B02001_004',
        asian_pop = 'B02001_005',
        pacific_islander_pop = 'B02001_006',
        other_race_pop = 'B02001_007',
        multiracial_pop = 'B02001_009',
        non_hispanic_pop = 'B03003_002',
        hispanic_pop = 'B03003_003',
        per_capita_income = 'B19301_001',
        ppl_in_poverty = 'B17021_002',
        tot_pop_over_25 = 'B15003_001',
        no_schooling = 'B15003_002',
        high_school_diploma = 'B15003_017',
        ged = 'B15003_018',
        some_college_lt_1yr = 'B15003_019',
        some_college_gt_1yr_no_deg = 'B15003_020',
        associates_deg = 'B15003_021',
        bachelors_deg = 'B15003_022',
        masters_deg = 'B15003_023',
        professional_deg = 'B15003_024',
        doctoral_deg = 'B15003_025'
      ), 
      state = 'NY', 
      county = 'Erie', 
      geometry = TRUE, 
      output = 'wide', 
      year = year
    ) %>%
      select(-ends_with('M', ignore.case = FALSE)) %>%
      rename_at(vars(ends_with('E', ignore.case = FALSE)), ~ str_replace(., 'E$', '')) %>%
      mutate(year = year,
             area_sqmi = as.numeric(st_area(.)) * 3.861E-7,
             tot_pop_density = tot_pop / area_sqmi,
             pop_25plus_density = tot_pop_over_25 / area_sqmi,
             total_income = per_capita_income * tot_pop,
             high_school_dip_equiv = high_school_diploma + ged,
             some_college_no_deg = some_college_lt_1yr + some_college_gt_1yr_no_deg,
             bachelors_deg_at_least = bachelors_deg + masters_deg + professional_deg + doctoral_deg,
             less_than_bachelors_deg = tot_pop_over_25 - bachelors_deg_at_least,
             less_than_high_school_grad = tot_pop_over_25 - (high_school_dip_equiv + some_college_no_deg + bachelors_deg_at_least)
      ) %>%
      select(-c(NAM, high_school_diploma, ged))
  )
}
