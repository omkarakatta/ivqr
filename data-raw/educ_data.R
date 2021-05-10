### Meta -------------------------
###
### Title: Prepare data from Autor et al. (2017) Replication
###
### Description: Clean and merge the data sets of Autor et al. (2017) for
### education-related replication
### See ~/BFI/5_IVQR_GP/notes/explore_autor2017_data/explore.rmd and
### ~/BFI/5_IVQR_GP/notes/apply_ivqr/apply_ivqr4.R for the initial drafts
###
### There are two main data sets in the Autor et al. (2017) replication files.
### The `detroit` data set is what they use to run the QR and IVQR regressions
### in Matlab. Second, the `ddata` data set is the raw Stata file containing the
### education variables of interest. Only the `ddata` data set has the education
### variables, so we would like to use the columns of the `ddata` data set that
### correspond to the QR and IVQR analysis performed on the `detroit` data.
### The `detroit` data set has unnamed column variables, so the Matlab code
### refers to their column numbers. Hence, our first job will be to assign
### column names to `detroit`. Then, we retrieve the column names based on the
### column numbers. Using these names, we select the appropriate columns in
### `ddata`. We also need to create the relevant education variables, i.e.,
### interactions between education and the endogeneous variables/instruments.
###
### Author: Omkar A. Katta
###

### Preliminaries -------------------------

library(dplyr)
library(readr) # read fixed-width files
library(haven) # read dta files from Stata
library(data.table) # read CSVs

### Prepare `detroit` -------------------------

detroit_path <- here::here("data-raw/detroit_full.txt")
detroit_raw <- read_fwf(detroit_path, fwf_empty(detroit_path))

detroit <- detroit_raw %>%
  rename(female = X1) %>%
  rename(wh_race = X2) %>%
  rename(oth_race = X3) %>%
  rename(age = X4) %>%
  rename(age2 = X5) %>%         # age2 does not appear to be the same age^2
  rename(wageadjh1_8 = X6) %>%
  rename(qtremph1_8 = X7) %>%
  rename(yq19994 = X8) %>%
  rename(yq20001 = X9) %>%
  rename(yq20002 = X10) %>%
  rename(yq20003 = X11) %>%
  rename(yq20004 = X12) %>%
  rename(yq20011 = X13) %>%
  rename(yq20012 = X14) %>%
  rename(yq20013 = X15) %>%
  rename(yq20014 = X16) %>%
  rename(yq20021 = X17) %>%
  rename(yq20022 = X18) %>%
  rename(yq20023 = X19) %>%
  rename(yq20024 = X20) %>%
  rename(yq20031 = X21) %>%
  rename(wageadj2_4 = X22) %>%
  rename(wageadj5_8 = X23) %>%
  rename(wageadj2_8 = X24) %>%
  rename(tempemp = X25) %>%
  rename(regemp = X26) %>%
  rename(emp = X27) %>%
  rename(js7_tern = X28) %>%
  rename(new_ear2_8 = X29) %>%
  rename(regern2_8 = X30) %>%
  rename(tmpern2_8 = X31) %>%
  rename(Idistyr_11 = X32) %>%
  rename(Idistyr_12 = X33) %>%
  rename(Idistyr_13 = X34) %>%
  rename(Idistyr_20 = X35) %>%
  rename(Idistyr_21 = X36) %>%
  rename(Idistyr_22 = X37) %>%
  rename(Idistyr_23 = X38) %>%
  rename(Idistyr_30 = X39) %>%
  rename(Idistyr_31 = X40) %>%
  rename(Idistyr_32 = X41) %>%
  rename(Idistyr_33 = X42) %>%
  rename(Idistyr_41 = X43) %>%
  rename(Idistyr_42 = X44) %>%
  rename(Idistyr_43 = X45) %>%
  rename(Idistyr_52 = X46) %>%
  rename(Idistyr_60 = X47) %>%
  rename(Idistyr_61 = X48) %>%
  rename(Idistyr_62 = X49) %>%
  rename(Idistyr_63 = X50) %>%
  rename(Idistyr_70 = X51) %>%
  rename(Idistyr_71 = X52) %>%
  rename(Idistyr_73 = X53) %>%
  rename(Idistyr_82 = X54) %>%
  rename(Idistyr_91 = X55) %>%
  rename(Idistyr_100 = X56) %>%
  rename(Idistyr_101 = X57) %>%
  rename(Idistyr_102 = X58) %>%
  rename(Idistyr_103 = X59) %>%
  rename(Idistyr_110 = X60) %>%
  rename(Idistyr_112 = X61) %>%
  rename(Idistyr_121 = X62) %>%
  rename(Idistyr_122 = X63) %>%
  rename(Idistyr_123 = X64) %>%
  rename(temp_pct = X65) %>%
  rename(reg_pct = X66) %>%
  rename(emp_pct = X67) %>%
  rename(regern2_4 = X68) %>%
  rename(tmpern2_4 = X69) %>%
  rename(regern5_8 = X70) %>%
  rename(tmpern5_8 = X71)

### Get names of variables from `detroit` -------------------------

# TODO: I *think* QR and IVQR code from Autor et al. missed some covariates that
# the OLS and 2SLS regressions used. How should we handle this?
# For now, I will not uses these covariates since we are most interested in QR
# and IVQR.

covs <- c(
  1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
  21, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 43, 44, 45, 46, 47, 48,
  49, 50, 51, 52, 53, 54, 56, 57, 58, 59, 60, 62, 63, 64
  # 42, 55, 61  # missed this in the matlab code
)

# get names of variables from detroit data set
y_name <- names(detroit)[24] # we are going to use wageadj2_8, not new_ear2_8
d_name <- names(detroit)[27]
z_name <- names(detroit)[67]
covs_name <- names(detroit)[covs]

educfe_name <- c("hsgrad", "hsplus") # hsless is going to be reference
dint_name <- paste0("emp_", educfe_name)
zint_name <- paste0("emp_pct_", educfe_name)

### Prepare `ddata` -------------------------

# The raw file is the result of applying the Stata code from Autor et al. (2017)
# on ddata_2021.dta which was given to us by the coauthors themselves.
# This contains education variables.

ddata_path <- here::here("data-raw/ddata_2021_cleaned.csv")
ddata_raw <- as.data.frame(fread(ddata_path))
ddata <- ddata_raw %>%
  mutate(emp_hsless = emp * hsless) %>%
  mutate(emp_hsgrad = emp * hsgrad) %>%
  mutate(emp_hsplus = emp * hsplus) %>%
  mutate(emp_educunkn = emp * educunkn) %>%
  mutate(emp_pct_hsless = emp_pct * hsless) %>%
  mutate(emp_pct_hsgrad = emp_pct * hsgrad) %>%
  mutate(emp_pct_hsplus = emp_pct * hsplus) %>%
  filter(educunkn == 0) # remove any observations whose education is unknown
colnames(ddata) <- gsub("^_", "", colnames(ddata)) # remove leading underscores

### Split `ddata` into matrices of interest -------------------------

# Outcome matrix
Y_educ <- as.matrix(ddata %>% select(all_of(y_name)))
colnames(Y_educ) <- y_name

# Endogeneous matrix
D_educ <- as.matrix(ddata %>% select(all_of(d_name)))
colnames(D_educ) <- d_name

# Instrument matrix
Z_educ <- as.matrix(ddata %>% select(all_of(z_name)))
colnames(Z_educ) <- z_name

# Covariate matrix
X_educ <- as.matrix(ddata %>% select(all_of(covs_name)))
colnames(X_educ) <- covs_name

# Education FE matrix
EducFE_educ <- as.matrix(ddata %>% select(all_of(educfe_name)))
colnames(EducFE_educ) <- educfe_name

# Education Interaction matrix
DInt_educ <- as.matrix(ddata %>% select(all_of(dint_name)))
colnames(DInt_educ) <- dint_name

# Instrument Interaction matrix
ZInt_educ <- as.matrix(ddata %>% select(all_of(zint_name)))
colnames(ZInt_educ) <- zint_name

### Save education data -------------------------

usethis::use_data(y_name, overwrite = TRUE)
usethis::use_data(d_name, overwrite = TRUE)
usethis::use_data(z_name, overwrite = TRUE)
usethis::use_data(covs_name, overwrite = TRUE)
usethis::use_data(educfe_name, overwrite = TRUE)
usethis::use_data(dint_name, overwrite = TRUE)
usethis::use_data(zint_name, overwrite = TRUE)

usethis::use_data(Y_educ, overwrite = TRUE)
usethis::use_data(D_educ, overwrite = TRUE)
usethis::use_data(Z_educ, overwrite = TRUE)
usethis::use_data(X_educ, overwrite = TRUE)
usethis::use_data(EducFE_educ, overwrite = TRUE)
usethis::use_data(DInt_educ, overwrite = TRUE)
usethis::use_data(ZInt_educ, overwrite = TRUE)
