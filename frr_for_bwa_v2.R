
rm(list=ls())
library(devtools)
library(eppasm)
library(rlang)
library(purrr)
library(dplyr)
library(stringr)
library(stringi)
library(dplyr)
setwd( "/Users/sak3699/Library/CloudStorage/OneDrive-HarvardUniversity/Postdoc Research Projects/Botswana/")
source('setpath.R')
bwa_path <- paste0(root, "spectrum/Botswana2023v4 WPP 02_03_2023 KOS.PJNZ")
df_area0_p <- readRDS(paste0(rds_dir, "/pool0_prev.rds"))  # Tsepamo


cd4_lab <-  c(">500", "350-499", "250-349", "200-249", "100-199", "50-99", "<50")
age_lab <- c( "15-16","17-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50+"  )
sex_lab <- c("male", "female")
art_lab <- c("art0mos", "art6mos", "art1yr")
hiv_lab <- c("hivn", "hivp")
agel <- 15:80

yr_start <- 1970
yr_end <- 2030
proj.years <- yr_start:yr_end
fp <- eppasm::prepare_directincid(bwa_path)
mod <- eppasm::simmod(fp)


dimnames(attr(mod, "hivpop")) <- list(cd4stage = cd4_lab , agegr = age_lab, sex = sex_lab, year =proj.years)
dimnames(attr(mod, "artpop")) <- list(artdur = art_lab, cd4stage = cd4_lab , agegr = age_lab, sex = sex_lab, year =proj.years)
dimnames(mod) <- list( age = agel, sex = sex_lab, c("hiv_neg", "hiv_pos"), year =proj.years)



## Optim fit ----
test_fun <- function (x,  mod, fp, ageseq, data, cd4, art ) {

  proj_years <- fp$ss$proj_start + seq_len(fp$ss$PROJ_YEARS) - 1L
  years <- proj_years
  df <- expand.grid(year = years, 
                    age =  ageseq,
                    sex = "female"
  )

  # take age groups: 15 to 49, both hIV negative and HIV positive people
  hivp <- mod [1:35, "female","hiv_pos",]
  hivn <- mod [1:35, "female","hiv_neg",]
  
  # sum by 5 year age group 
  hivp <- apply(hivp, 2, function (x) tapply(x, (seq_along(  x) - 1) %/% 5, sum))
  hivn <- apply(hivn, 2, function (x) tapply(x, (seq_along(  x) - 1) %/% 5, sum))
  
  # convert both from long to wide to get a dataset with year, age group, and number of HIV negative and positive women 
  hivpos <- as.data.frame(  hivp )  %>%
    tidyr::pivot_longer( cols = everything(), 
                         names_to = "year",      
                         values_to = "hivpos") %>%
    mutate (age = rep (  ageseq, length(years)))
  
  hivneg <- as.data.frame(  hivn )  %>%
    tidyr::pivot_longer( cols = everything(), 
                         names_to = "year",   
                         values_to = "hivneg") %>%
    mutate (age = rep (  ageseq, length(years)))
  
  hivs <- left_join(hivpos, hivneg)
  
  ## Calculate age-specific FRR given the CD4 and ART duration distribution
  
  ## sum the 15-17 and 15-19 groups for women living with HIV 
  hivpop_fert <- attr(mod, "hivpop")[ , fp$ss$h.fert.idx, fp$ss$f.idx, ]
  mat_hiv <- array (c (0, 0,0, 0),dim =c(length(cd4_lab),length(age_lab[1:8])-1,length(proj.years)), dimnames = list (cd4_lab, c("15-19", age_lab[3:8]), proj.years ))
  
  mat_hiv [,"15-19",] <- apply(  hivpop_fert [, 1:2,], c(1,3), sum) ## s
  mat_hiv [,2:7,] <-  hivpop_fert [, 3:8,]
  
  
  ## sum the 15-17 and 15-19 groups for women living with HIV by ART status 
  artpop_fert <- attr(mod, "artpop")[ , , fp$ss$h.fert.idx, fp$ss$f.idx, ]
  mat_art <- array (c (0, 0,0, 0),dim =c(3, length(cd4_lab),length(age_lab[1:8])-1,length(proj.years)), dimnames = list (c("art0mos", "art6mos", "art1yr" ), cd4_lab, c("15-19", age_lab[3:8]), proj.years ))
  
  mat_art [,,"15-19",] <- apply(    artpop_fert [,,1:2,], c(1,2,4), sum) 
  mat_art  [,,2:7,] <-    artpop_fert [,,3:8,]
  
  # for frr by CD4 count, remove the frr for 15-16 since its the same as 17-18 anyways 
  if (is.null(cd4)) {
    
  frr_cd4 <- fp$frr_cd4[, -1, , drop = FALSE]
  dimnames(   frr_cd4  ) <- list(cd4stage = cd4_lab, c("15-19", age_lab[3:8]), year =proj.years)
  
  } else {
  frr_cd4 <- array(x, dim = c(length(cd4_lab),length(age_lab[1:8])-1,length(proj.years))) 
  dimnames(  frr_cd4 ) <- list(cd4stage = cd4_lab, c("15-19", age_lab[3:8]), year =proj.years)
  }

  
  # for frr by ART count, remove the frr for 15-16 since its the same as 17-18 anyways 
  if (is.null(art)) { 
  
    frr_art <- fp$frr_art[,,-1,, drop = FALSE]
    dimnames(  frr_art) <-  list(art = art_lab, cd4stage = cd4_lab,  age = c("15-19", age_lab[3:8]), year =proj.years)
  
  } else { 
    frr_art <- array(x, dim = c(length (art_lab), length(cd4_lab),length(age_lab[1:8])-1,length(proj.years))) 
    dimnames(  frr_art) <- list(art = art_lab, cd4stage = cd4_lab,  age = c("15-19", age_lab[3:8]), year =proj.years)
    
    }
  

  # calculate the weighted FRR
  ha_frr <- (colSums(mat_hiv * frr_cd4) + colSums(mat_art * frr_art,,2)) / (colSums(mat_hiv) + colSums(mat_art ,,2))
  
  # ASFR is same by 5 yeag age group, so just select every 5th row 
  asfr_new <- fp$asfr[seq(1, nrow(fp$asfr), by = 5), ]
  
  # convert ASFR from long to wide 
  asfr <- as.data.frame(asfr_new)  %>%
    tidyr::pivot_longer( cols = everything(), # Select columns starting with 'Var'
                         names_to = "year",      # Name the new variable column
                         values_to = "asfr") %>%
    mutate (age = rep ( ageseq, length(years)))
  
  # covert wighted FRR from long to wide 
  ha_frr_new <- as.data.frame( ha_frr)  %>%
    tidyr::pivot_longer( cols = everything(), # Select columns starting with 'Var'
                         names_to = "year",      # Name the new variable column
                         values_to = "asfr") %>%
    mutate (age = rep ( ageseq, length(years)))
  
  # calculate the number of births among both HIV positve and hIV negtive women 
  births_a <-   asfr$asfr * (hivs$hivneg + hivs$hivpos)
  
  # HIV prevalence among pregnant women 
  pregprev_a <- 1 - hivs$hivneg / (hivs$hivneg + ha_frr_new$asfr * hivs$hivpos)
  df$pregprev_a <- pregprev_a 
  
  
  ## for now only age stratified 
  df <- df %>%
  mutate (age_group = case_when(age =="15" ~ "Y015_019", 
                                age == "20" ~ "Y020_024", 
                                age == "25" ~ "Y025_029", 
                                age == "30" ~ "Y030_034", 
                                age == "35" ~ "Y035_039", 
                                age == "40" ~ "Y040_044", 
                                age == "45" ~ "Y045_049", 
                                TRUE ~ NA)) 
  dt_out <- left_join( data, df) %>% filter (!is.na (pregprev_a))
  bin_out <- dbinom(dt_out$n_pos , dt_out$n_test, dt_out$pregprev_a, log = T )
  print(sum(bin_out))
  
  return(sum(bin_out))
  
}

##5 

 #test_cd4_opt <- optim(par = rep(1,2), fn = test_fun,  cd4=1, art =NULL,
  #                     mod = mod,
  #                     fp = fp,
  #                     ageseq = c(15, 20, 25, 30, 35, 40, 45), 
  #                     control = list(fnscale = -1, reltol = 1e-8),
  #                     data = df_area0_p)
 

## estimate 4 paramters ---
test_art_opt <- optim(par = rep(1,4) , fn = test_fun, cd4=NULL, art =1,
                      mod = mod,
                      fp = fp,
                      ageseq = c(15, 20, 25, 30, 35, 40, 45), 
                      control = list(fnscale = -1, reltol = 1e-8),
                      data = df_area0_p)

## Change FRR data with these parameters and re-run the simulation  ----
frr_art <- fp$frr_art[,,-1,, drop = FALSE]
dimnames(  frr_art) <-  list(art = art_lab, cd4stage = cd4_lab,  age = c("15-19", age_lab[3:8]), year =proj.years)

frr_art_new <- frr_art 
frr_art_new[,,c ("25-29"),] <- frr_art_new[,,c ("25-29"),] * 1.18 # this fits better but unable eto obtain the number from the optim()
frr_art_new[,,c ("30-34"),] <- frr_art_new[,,c ("30-34"),] * test_art_opt$par[2]
frr_art_new[,,c ("35-39"),] <- frr_art_new[,,c ("35-39"),] * test_art_opt$par[3]
frr_art_new[,,c ("40-44"),] <- frr_art_new[,,c ("40-44"),] * test_art_opt$par[3]





mod = mod;
fp = fp;
ageseq = c(15, 20, 25, 30, 35, 40, 45)

proj_years <- fp$ss$proj_start + seq_len(fp$ss$PROJ_YEARS) - 1L
years <- proj_years
df <- expand.grid(year = years, 
                  age =  ageseq,
                  sex = "female"
)

# take age groups: 15 to 49, both hIV negative and HIV positive people
hivp <- mod [1:35, "female","hiv_pos",]
hivn <- mod [1:35, "female","hiv_neg",]

# sum by 5 year age group 
hivp <- apply(hivp, 2, function (x) tapply(x, (seq_along(  x) - 1) %/% 5, sum))
hivn <- apply(hivn, 2, function (x) tapply(x, (seq_along(  x) - 1) %/% 5, sum))

# convert both from long to wide to get a dataset with year, age group, and number of HIV negative and positive women 
hivpos <- as.data.frame(  hivp )  %>%
  tidyr::pivot_longer( cols = everything(), 
                       names_to = "year",      
                       values_to = "hivpos") %>%
  mutate (age = rep (  ageseq, length(years)))

hivneg <- as.data.frame(  hivn )  %>%
  tidyr::pivot_longer( cols = everything(), 
                       names_to = "year",   
                       values_to = "hivneg") %>%
  mutate (age = rep (  ageseq, length(years)))

hivs <- left_join(hivpos, hivneg)
#1 (15), 6 (20), 11 (25), 16 (30),  21 (35), 26 (40), 31 ( 45), 35 (49 )


## Calculate age-specific FRR given the CD4 and ART duration distribution

## sum the 15-17 and 15-19 groups for women living with HIV 
hivpop_fert <- attr(mod, "hivpop")[ , fp$ss$h.fert.idx, fp$ss$f.idx, ]
mat_hiv <- array (c (0, 0,0, 0),dim =c(length(cd4_lab),length(age_lab[1:8])-1,length(proj.years)), dimnames = list (cd4_lab, c("15-19", age_lab[3:8]), proj.years ))

mat_hiv [,"15-19",] <- apply(  hivpop_fert [, 1:2,], c(1,3), sum) ## s
mat_hiv [,2:7,] <-  hivpop_fert [, 3:8,]


## sum the 15-17 and 15-19 groups for women living with HIV by ART status 
artpop_fert <- attr(mod, "artpop")[ , , fp$ss$h.fert.idx, fp$ss$f.idx, ]
mat_art <- array (c (0, 0,0, 0),dim =c(3, length(cd4_lab),length(age_lab[1:8])-1,length(proj.years)), dimnames = list (c("art0mos", "art6mos", "art1yr" ), cd4_lab, c("15-19", age_lab[3:8]), proj.years ))

mat_art [,,"15-19",] <- apply(    artpop_fert [,,1:2,], c(1,2,4), sum) 
mat_art  [,,2:7,] <-    artpop_fert [,,3:8,]


frr_art <- frr_art_new
frr_cd4 <- frr_cd4_new

# calculate the weighted FRR
ha_frr <- (colSums(mat_hiv * frr_cd4) + colSums(mat_art * frr_art,,2)) / (colSums(mat_hiv) + colSums(mat_art ,,2))

# ASFR is same by 5 yeag age group, so just select every 5th row 
asfr_new <- fp$asfr[seq(1, nrow(fp$asfr), by = 5), ]

# convert ASFR from long to wide 
asfr <- as.data.frame(asfr_new)  %>%
  tidyr::pivot_longer( cols = everything(), # Select columns starting with 'Var'
                       names_to = "year",      # Name the new variable column
                       values_to = "asfr") %>%
  mutate (age = rep ( ageseq, length(years)))

# covert wighted FRR from long to wide 
ha_frr_new <- as.data.frame( ha_frr)  %>%
  tidyr::pivot_longer( cols = everything(), # Select columns starting with 'Var'
                       names_to = "year",      # Name the new variable column
                       values_to = "asfr") %>%
  mutate (age = rep ( ageseq, length(years)))

# calclulate the number of births among both HIV positve and hIV negtive women 
births_a <-   asfr$asfr * (hivs$hivneg + hivs$hivpos)

# HIV prevalence among pregnant women 
pregprev_a <- 1 - hivs$hivneg / (hivs$hivneg + ha_frr_new$asfr * hivs$hivpos)
df$pregprev_a <- pregprev_a 

#print (sum(df$pregprev_a, na.rm = T))

## for now only age stratified 
df_n <- df %>%
  mutate (age_group = case_when(age =="15" ~ "Y015_019", 
                                age == "20" ~ "Y020_024", 
                                age == "25" ~ "Y025_029", 
                                age == "30" ~ "Y030_034", 
                                age == "35" ~ "Y035_039", 
                                age == "40" ~ "Y040_044", 
                                age == "45" ~ "Y045_049", 
                                TRUE ~ NA))  %>%
  dplyr::rename (prevalence = pregprev_a) %>%
  filter (year >=2014 ) %>%
  select (age_group, year, prevalence ) %>%
  mutate (scenario = "Spectrum")




## bind to tsepamo
pregprev_tse <- df_area0_p %>%
  dplyr:: select(year, age_group, tse_prevalence, scenario) %>%
  rename (prevalence  = tse_prevalence)

pregprev <- rbind(df_n, pregprev_tse)
pregprev <- pregprev  %>%
  filter (year < 2023 & age_group != "Y015_049")


## plot it ----

ggplot( data = pregprev, aes(x = year, y = prevalence, color = scenario)  ) +
  geom_point() +
  facet_wrap (~ age_group) +
  stat_smooth(se = F)+ 
  scale_color_manual(values = c('#46ACC8', '#B40F20'), name = "Data source") +
  scale_x_continuous(   breaks = seq(2014, 2025, by = 2)) +
  labs(
    x = "Year", 
    y = "HIV prevalence among pregnant women")  +
  theme_bw() +
  theme (axis.title.x = element_text(size = 12, face = "bold"),
         axis.title.y = element_text(size = 12, face = "bold" ))



ggsave(file = paste0(figs_dir, "/", "pregprev_hiv_changedfrr.png"), device = "png", 
       width = 24, height =18, units = c("cm"), dpi = 300, limitsize = TRUE )








## Manual fit (ignore)  -----


extract_eppasm_pregprev_change <- function(mod, fp, ageseq,
                                           cd4_1519,
                                           cd4_2024,
                                           cd4_2529,
                                           cd4_3034,
                                           cd4_3539,
                                           cd4_4044,
                                           cd4_4549,
                                           art_1519,
                                           art_2024,
                                           art_2529,
                                           art_3034,
                                           art_3539,
                                           art_4044,
                                           art_4549 ) {
  
  proj_years <- fp$ss$proj_start + seq_len(fp$ss$PROJ_YEARS) - 1L
  years <- proj_years
  df <- expand.grid(year = years, 
                    age =  ageseq,
                    sex = "female"
  )
  # View(mod [1:35, "female","hiv_pos",])
  
  # take age groups: 15 to 49, both hIV negative and HIV positive people
  hivp <- mod [1:35, "female","hiv_pos",]
  hivn <- mod [1:35, "female","hiv_neg",]
  
  # sum by 5 year age group 
  hivp <- apply(hivp, 2, function (x) tapply(x, (seq_along(  x) - 1) %/% 5, sum))
  hivn <- apply(hivn, 2, function (x) tapply(x, (seq_along(  x) - 1) %/% 5, sum))
  
  # convert both from long to wide to get a dataset with year, age group, and number of HIV negative and positive women 
  hivpos <- as.data.frame(  hivp )  %>%
    tidyr::pivot_longer( cols = everything(), 
                         names_to = "year",      
                         values_to = "hivpos") %>%
    mutate (age = rep (  ageseq, length(years)))
  
  hivneg <- as.data.frame(  hivn )  %>%
    tidyr::pivot_longer( cols = everything(), 
                         names_to = "year",   
                         values_to = "hivneg") %>%
    mutate (age = rep (  ageseq, length(years)))
  
  hivs <- left_join(hivpos, hivneg)
  #1 (15), 6 (20), 11 (25), 16 (30),  21 (35), 26 (40), 31 ( 45), 35 (49 )
  
  
  ## Calculate age-specific FRR given the CD4 and ART duration distribution
  
  ## sum the 15-17 and 15-19 groups for women living with HIV 
  hivpop_fert <- attr(mod, "hivpop")[ , fp$ss$h.fert.idx, fp$ss$f.idx, ]
  mat_hiv <- array (c (0, 0,0, 0),dim =c(length(cd4_lab),length(age_lab[1:8])-1,length(proj.years)), dimnames = list (cd4_lab, c("15-19", age_lab[3:8]), proj.years ))
  
  mat_hiv [,"15-19",] <- apply(  hivpop_fert [, 1:2,], c(1,3), sum) ## s
  mat_hiv [,2:7,] <-  hivpop_fert [, 3:8,]
  
  
  ## sum the 15-17 and 15-19 groups for women living with HIV by ART status 
  artpop_fert <- attr(mod, "artpop")[ , , fp$ss$h.fert.idx, fp$ss$f.idx, ]
  mat_art <- array (c (0, 0,0, 0),dim =c(3, length(cd4_lab),length(age_lab[1:8])-1,length(proj.years)), dimnames = list (c("art0mos", "art6mos", "art1yr" ), cd4_lab, c("15-19", age_lab[3:8]), proj.years ))
  
  mat_art [,,"15-19",] <- apply(    artpop_fert [,,1:2,], c(1,2,4), sum) 
  mat_art  [,,2:7,] <-    artpop_fert [,,3:8,]
  
  # for frr by CD4 count, remove the frr for 15-16 since its the same as 17-18 anyways 
  
  frr_cd4 <- fp$frr_cd4[, -1, , drop = FALSE]
  dimnames(   frr_cd4  ) <- list(cd4stage = cd4_lab, c("15-19", age_lab[3:8]), year =proj.years)
  
  
  # Reduce ages 25-29 and higher by a certain % 
  
  redcd4 <- array(NA, dim = c(length(cd4_lab), 7, length(proj.years)))
  dimnames(   redcd4  ) <- list(cd4stage = cd4_lab, c("15-19", age_lab[3:8]), year =proj.years)
  redcd4[,"15-19",] <- cd4_1519
  redcd4[,"20-24",] <- cd4_2024
  redcd4[,"25-29",] <- cd4_2529
  redcd4[,"30-34",] <- cd4_3034
  redcd4[,"35-39",] <- cd4_3539
  redcd4[,"40-44",] <- cd4_4044
  redcd4[,"45-49",] <- cd4_4549
  
  
  frr_cd4 <-  frr_cd4*redcd4
  
  
  # for frr by ART count, remove the frr for 15-16 since its the same as 17-18 anyways 
  frr_art <- fp$frr_art[,,-1,, drop = FALSE]
  dimnames(  frr_art) <-  list(art = art_lab, cd4stage = cd4_lab,  age = c("15-19", age_lab[3:8]), year =proj.years)
  
  # reduce art: 
  
  redart <- array(NA, dim = c(3, length(cd4_lab), 7, length(proj.years)))
  dimnames(   redart ) <- list( art = art_lab, cd4stage = cd4_lab, c("15-19", age_lab[3:8]), year =proj.years)
  redart[,,"15-19",] <- art_1519
  redart[,,"20-24",] <- art_2024
  redart[,,"25-29",] <- art_2529
  redart[,,"30-34",] <- art_3034
  redart[,,"35-39",] <- art_3539
  redart[,,"40-44",] <- art_4044
  redart[,,"45-49",] <- art_4549
  
  
  frr_art <-   frr_art*redart
  
  
  
  # calculate the weighted FRR
  ha_frr <- (colSums(mat_hiv * frr_cd4) + colSums(mat_art * frr_art,,2)) / (colSums(mat_hiv) + colSums(mat_art ,,2))
  
  # ASFR is same by 5 yeag age group, so just select every 5th row 
  asfr_new <- fp$asfr[seq(1, nrow(fp$asfr), by = 5), ]
  
  # convert ASFR from long to wide 
  asfr <- as.data.frame(asfr_new)  %>%
    tidyr::pivot_longer( cols = everything(), # Select columns starting with 'Var'
                         names_to = "year",      # Name the new variable column
                         values_to = "asfr") %>%
    mutate (age = rep ( ageseq, length(years)))
  
  # covert wighted FRR from long to wide 
  ha_frr_new <- as.data.frame( ha_frr)  %>%
    tidyr::pivot_longer( cols = everything(), # Select columns starting with 'Var'
                         names_to = "year",      # Name the new variable column
                         values_to = "asfr") %>%
    mutate (age = rep ( ageseq, length(years)))
  
  # calclulate the number of births among both HIV positve and hIV negtive women 
  births_a <-   asfr$asfr * (hivs$hivneg + hivs$hivpos)
  
  # HIV prevalence among pregnant women 
  pregprev_a <- 1 - hivs$hivneg / (hivs$hivneg + ha_frr_new$asfr * hivs$hivpos)
  df$pregprev_a <- pregprev_a 
  
  return(df)
}



pregprev_5age <- extract_eppasm_pregprev_change (mod = mod, fp = fp, ageseq = c(15, 20, 25, 30, 35, 40, 45), 
                                                 cd4_1519 = 1.1,
                                                 cd4_2024 = 1,
                                                 cd4_2529 = 1.15,
                                                 cd4_3034 = 1,
                                                 cd4_3539 = 1.3,
                                                 cd4_4044 = 1,
                                                 cd4_4549 = 1,
                                                 art_1519 = 1.1,
                                                 art_2024 = 1,
                                                 art_2529 = 1.15,
                                                 art_3034 = 1.3,
                                                 art_3539 = 1.65,
                                                 art_4044 = 1.6,
                                                 art_4549 = 1)








