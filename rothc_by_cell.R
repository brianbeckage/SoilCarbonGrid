
rothc_by_cell <- function(cell_id, gcm, gis_dat_path, gcm_dat_path, lm_dat_path, out_dat_path) {

  ## Load required libraries
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(readxl))
  # devtools::install_github('MPIBGC-TEE/SoilR-exp/pkg')
  suppressPackageStartupMessages(library(SoilR))
  
  
  ## Set global parameters in list
  params_ls <- list(cell_id = cell_id,
                    gcm = gcm,
                    gis_dat_path = gis_dat_path,
                    gcm_dat_path = gcm_dat_path,
                    lm_dat_path = lm_dat_path ,
                    lu_labs = c("crops", "hay", "pasture"),
                    t0_yr = 2022,
                    n_yrs_init = 750,
                    n_yrs = 77, # sim from 2022 - 2099 (end of epscor data)
                    soil_thick = 30,
                    dpm_rpm_ratio = "default")
  
  # create vectors of monthly timesteps for spinup and simulation
  params_ls$tsteps_spinup <- seq(1/12, params_ls$n_yrs_init, by=1/12)
  params_ls$tsteps <- seq(1/12, params_ls$n_yrs, by=1/12)
  
  # initialize lists to store input and output data
  gis_dat_ls <- gcm_dat_ls <- lm_dat_ls <- scn_dat_ls <- res_ls <- list()
  
  
  ## Calculate and store key input data for RothC from gis_dat in the "gis_dat_ls" list object
  pop_gis_dat_ls <- function(cell_id, gis_dat_path) {
    
    # load gis data file row of interest
    gis_dat_row <- read_xlsx(gis_dat_path) %>% filter(.data[["cell_id"]] == .env[["cell_id"]])
    
    # pull relevant data into gis_dat_ls
    gis_dat_ls$cell_id <- gis_dat_row$cell_id
    gis_dat_ls$clim$evap <- gis_dat_row %>% select(contains("gldas")) %>% unlist(use.names = F)  
    gis_dat_ls$soil <- list(clay = gis_dat_row$clay,
                            bulk_dens = gis_dat_row$bulk_dens,
                            som = list(crops = gis_dat_row$som_crops,
                                       hay = gis_dat_row$som_hay,
                                       pasture = gis_dat_row$som_pasture)
    ) 
    
    vb_fact <- 0.5 # modern revision of Van Bemmelen factor
    
    # calculate soil organic carbon in mg/ha for each land use in cell
    for (lu_lab in params_ls$lu_labs) {
      som <- gis_dat_ls$soil$som[[lu_lab]] # percent SOM
      soc <- som * vb_fact # times modified V.B. factor
      gis_dat_ls$soil$soc_mg_ha[[lu_lab]] <- gis_dat_ls$soil$bulk_dens * 3000 * soc / 100 # use bulk density to find soc mg/ha
    }
    
    gis_dat_ls
  }
  
  gis_dat_ls <- pop_gis_dat_ls(params_ls$cell_id, params_ls$gis_dat_path)
  
  
  ## Import EPSCoR GCM climate data
  pop_gcm_dat_ls <- function(cell_id, gcm) {
    
    gcm_dat <- read_csv(paste(params_ls$gcm_dat_path, cell_id, ".csv", sep = ""), show_col_types = F)
    gcm_dat_ls$P <- gcm_dat$P * 30.437 # rectify precip values from mm/day to mm/month
    gcm_dat_ls$Tavg <- colMeans(rbind(gcm_dat$Tmin, gcm_dat$Tmax))
    
    gcm_dat_ls
  }
  
  gcm_dat_ls <- pop_gcm_dat_ls(params_ls$cell_id, params_ls$gcm)
  
  if (sum(gcm_dat_ls[["P"]]) == 0 || sum(gcm_dat_ls[["Tavg"]]) == 0) {
    print(paste("Error: No GCM data available for cell #", params_ls$cell_id, ", ", params_ls$gcm, sep = ""))
    break
  }
  
  ## Set land management parameters from file
  pop_lm_dat_ls <- function(lm_dat_file) {
    
    # load data from excel file into list object
    lm_excel_ls <- list(crops = read_excel(lm_dat_file, sheet = "Crops") %>% suppressMessages(),
                        hay = read_excel(lm_dat_file, sheet = "Hay") %>% suppressMessages(),
                        pasture = read_excel(lm_dat_file, sheet = "Pasture") %>% suppressMessages()
    )
    
    for (lu_lab in names(lm_excel_ls)) {
      lm_dat <- lm_excel_ls[[lu_lab]]
      lm_dat_ls[[lu_lab]] <- list(plt_res_c = as.numeric(lm_dat[2,2:13]),
                                  manure_c = as.numeric(lm_dat[3,2:13]),
                                  soil_cov = as.numeric(lm_dat[4,2:13]))
    }
    
    lm_dat_ls
  }
  
  lm_dat_ls <- pop_lm_dat_ls(params_ls$lm_dat_path)
  
  
  ## Calculate effects of temp. and precip. on decomposition rates
  pop_scn_dat_ls_xi <- function(gis_dat_ls, lm_dat_ls, scn_dat_ls, params_ls) {
    
    # populate scn_dat_ls with tstep, gis, fT, fW, and xi data
    for (lu_lab in params_ls$lu_labs) {
      
      # soil covered or bare each month? (0,1) - same for all treatments
      s_cov <- lm_dat_ls[[lu_lab]]$soil_cov
      
      ## SPINUP
      # spinup temp and prcp vectors are monthly average over 30 year period prior to t0, repeated for all years of spinup period
      temp_vec <- tibble(mo = rep(seq(1,12), 30),
                         t = gcm_dat_ls$Tavg[((41*12)+1):(71*12)]) %>% # 1992 - 2022 (30 yrs prior to sim t0)) 
        group_by(mo) %>%
        summarize(t_mo_avg = mean(t)) %>%
        pull(t_mo_avg) %>% unlist() %>%
        rep(params_ls$n_yrs_init)
      
      prcp_vec <- tibble(mo = rep(seq(1,12), 30),
                         p = gcm_dat_ls$P[((41*12)+1):(71*12)]) %>% # 1992 - 2022 (30 yrs prior to sim t0)) 
        group_by(mo) %>%
        summarize(p_mo_avg = mean(p)) %>%
        pull(p_mo_avg) %>% unlist() %>%
        rep(params_ls$n_yrs_init)
      
      # evap is from gis_dat
      evap_vec <- rep(gis_dat_ls$clim$evap,
                      params_ls$n_yrs_init)
      
      # fT vector; temperature effects per month
      fT <- fT.RothC(temp_vec)
      
      # fW vector; moisture effects per month
      # under vegetated conditions
      fW_veg <- fW.RothC(P = prcp_vec,
                         E = evap_vec,
                         S.Thick = params_ls$soil_thick,
                         pClay = gis_dat_ls$soil$clay,
                         pE = 1, bare = F)$b
      
      # under bare soil conditions
      fW_bare <- fW.RothC(P = prcp_vec,
                          E = evap_vec,
                          S.Thick = params_ls$soil_thick,
                          pClay = gis_dat_ls$soil$clay,
                          pE = 1, bare = T)$b
      
      # combine vegetated and bare months
      fW <- numeric(params_ls$n_yrs_init * 12)
      fW[s_cov == T] <- fW_veg[s_cov == T]
      fW[s_cov == F] <- fW_bare[s_cov == F]
      
      # compute xi vector by multiplying fT*fW
      scn_dat_ls[[lu_lab]]$xi_spinup <- fT * fW
      
      
      ## SIMULATION
      # temp and prcp projections going forward from 2022
      temp_vec <- gcm_dat_ls$Tavg %>% tail(params_ls$n_yrs * 12) # sim values = last n_yrs * 12
      
      prcp_vec <- gcm_dat_ls$P %>% tail(params_ls$n_yrs * 12) # sim values = last n_yrs * 12
      
      # evap is from gis_dat
      evap_vec <- rep(gis_dat_ls$clim$evap,
                      params_ls$n_yrs)
      
      # fT vector; temperature effects per month
      fT <- fT.RothC(temp_vec)
      
      # fW vector; moisture effects per month
      # under vegetated conditions
      fW_veg <- fW.RothC(P=prcp_vec,
                         E=evap_vec,
                         S.Thick = params_ls$soil_thick,
                         pClay = gis_dat_ls$soil$clay,
                         pE = 1, bare = F)$b
      
      # under bare soil conditions
      fW_bare <- fW.RothC(P=prcp_vec,
                          E=evap_vec,
                          S.Thick = params_ls$soil_thick,
                          pClay = gis_dat_ls$soil$clay,
                          pE = 1, bare = T)$b
      
      # combine vegetated and bare months
      fW <- numeric(params_ls$n_yrs * 12)
      fW[s_cov == T] <- fW_veg[s_cov == T]
      fW[s_cov == F] <- fW_bare[s_cov == F]
      
      # compute xi vector by multiplying fT*fW
      scn_dat_ls[[lu_lab]]$xi_sim <- fT * fW
    }
    
    scn_dat_ls
  }
  
  scn_dat_ls <- pop_scn_dat_ls_xi(gis_dat_ls, lm_dat_ls, scn_dat_ls, params_ls)
  
  
  ## Calculate IOM proportion
  pop_c_stocks_init <- function(scn_dat_ls, gis_dat_ls, params_ls) {
    
    for (lu_lab in params_ls$lu_labs) {
      
      soc_mg_ha <- gis_dat_ls$soil$soc_mg_ha[[lu_lab]] # pull soc (in Mg/Ha) from gis_dat_ls
      
      FallIOM <- 0.049 * soc_mg_ha^(1.139) # Falloon IOM function - gives estimate of prop. soc that's inert
      
      # DPM, RPM, BIO, HUM, IOM
      scn_dat_ls[[lu_lab]]$c_stocks_init <- c(0, 0, 0, 0, FallIOM)
    }
    
    scn_dat_ls
  }
  
  scn_dat_ls <- pop_c_stocks_init(scn_dat_ls, gis_dat_ls, params_ls)
  
  
  ## Define spinup runs
  run_model_spinup <- function(lu_lab, gis_dat_ls, lm_dat_ls, scn_dat_ls, params_ls, plt_res_c = "default") {
    
    # if simply spinning up the model, use value in scn_dat_ls; when used in fitness fxn can specify plt_res_c
    if (!is.numeric(plt_res_c) && plt_res_c == "default") {
      plt_res_c <- scn_dat_ls[[lu_lab]]$plt_res_c_calib
    } 
    
    # format input data to df; normalize to 12 months
    plt_res_c_df <- data.frame(years = params_ls$tsteps_spinup,
                               plt_res_c = rep(plt_res_c * 12,
                                               params_ls$n_yrs_init))
    
    manure_c_df <- data.frame(years = params_ls$tsteps_spinup,
                              manure_c = rep(lm_dat_ls[[lu_lab]]$manure_c * 12,
                                             params_ls$n_yrs_init))
    
    xi_df <- data.frame(years = params_ls$tsteps_spinup,
                        xi = scn_dat_ls[[lu_lab]]$xi_spinup)
    
    # load the model
    model <- RothCModel(t = params_ls$tsteps_spinup,
                        C0 = scn_dat_ls[[lu_lab]]$c_stocks_init,
                        In = plt_res_c_df,
                        FYM = manure_c_df,
                        clay = gis_dat_ls$soil$clay,
                        xi = xi_df)
    
    c_stocks <- getC(model) # calc stocks for each pool per month
    c_stocks
  }
  
  
  ## Adjust below-ground plant matter input levels so baseline C stock matches empirical observations
  # define fitness function for optimization
  fit_fxn <- function (bgc, lu_lab, gis_dat_ls, lm_dat_ls, scn_dat_ls, params_ls) {

    # assume below-ground c is only returned to soil when ground vegetated
    grw_seas <- lm_dat_ls[[lu_lab]]$soil_cov # e.g. c(0,0,0,1,1,1,1,1,1,1,1,0)
    bgc_vec <- grw_seas * bgc / sum(grw_seas)
    
    # sum above and below ground plant carbon inputs
    plt_res_c <- lm_dat_ls[[lu_lab]]$plt_res_c + bgc_vec
    
    c_stocks <- run_model_spinup(lu_lab, gis_dat_ls, lm_dat_ls, scn_dat_ls, params_ls, plt_res_c = plt_res_c)
    
    soc_mg_ha_spinup <- tail(c_stocks, n = 12) %>% # SOC (Mg/Ha) after spinup, avg of last yr.
      rowSums() %>% mean()
    
    soc_mg_ha_emp <- gis_dat_ls$soil$soc_mg_ha[[lu_lab]] # empirical SOC (Mg/Ha)
    
    abs(soc_mg_ha_spinup - soc_mg_ha_emp) # value to minimize is absolute difference between empirical and model C stocks
  }
  
  # run spinup optimizations for each land use
  for (lu_lab in params_ls$lu_labs) {

    # run optimization
    bgc_opt <- optimize(fit_fxn, c(0, 13), tol = .1, # can change search interval and tolerance for optimization here
                        lu_lab, gis_dat_ls, lm_dat_ls, scn_dat_ls, params_ls) # params to pass to fit fxn
    
    if (bgc_opt$objective > 1) {
      print(paste("Warning: Cell #", params_ls$cell_id, ", ", lu_lab, ", did not converge during spinup: diff = ", round(bgc_opt$objective, digits = 2), sep = ""))
    }
    
    # store bgc/yr in scn_dat_ls for later analysis
    scn_dat_ls[[lu_lab]]$bgc_yr <- bgc_opt$minimum
    
    # store resultant calibrated plt_res_c in scn_dat_ls
    # assume below-ground c is only returned to soil when ground vegetated
    grw_seas <- lm_dat_ls[[lu_lab]]$soil_cov # e.g. c(0,0,0,1,1,1,1,1,1,1,1,0)
    bgc_vec <- grw_seas * bgc_opt$minimum / sum(grw_seas)
    
    # re-calc plt_res_c in scn_dat_ls to include below-ground c
    scn_dat_ls[[lu_lab]]$plt_res_c_calib <- lm_dat_ls[[lu_lab]]$plt_res_c + bgc_vec
  }
  
  # run final spinups and save results
  for (lu_lab in params_ls$lu_labs) {
    res_ls[[lu_lab]]$spinup$c_stocks <- run_model_spinup(lu_lab, gis_dat_ls, lm_dat_ls, scn_dat_ls, params_ls)
    
    scn_dat_ls[[lu_lab]]$c_stocks_init <- as.numeric(tail(res_ls[[lu_lab]]$spinup$c_stocks, n = 1))
  }

  ## Execute simulation runs
  # define simulation function
  run_model_sim <- function(lu_lab_old, lu_lab_new, scn_dat_ls, gis_dat_ls, res_ls, params_ls) {
    
    # initial c stocks are based on spinup of original land use / management
    c_stocks_init <- scn_dat_ls[[lu_lab_old]]$c_stocks_init
    
    # rest of data is based on new land use / management
    plt_res_c <- rep(scn_dat_ls[[lu_lab_new]]$plt_res_c_calib * 12, 
                     params_ls$n_yrs)
    
    # format input data to df  
    plt_res_c_df <- data.frame(years = params_ls$tsteps,
                               plt_res_c = plt_res_c)
    
    manure_c_df <- data.frame(years = params_ls$tsteps,
                              manure_c = rep(lm_dat_ls[[lu_lab_new]]$manure_c * 12,
                                             params_ls$n_yrs))
    
    xi_df <- data.frame(years = params_ls$tsteps,
                        xi = scn_dat_ls[[lu_lab_new]]$xi_sim)
    
    # load the model
    model <- RothCModel(t = params_ls$tsteps,
                        C0 = c_stocks_init,
                        In = plt_res_c_df,
                        FYM = manure_c_df,
                        clay = gis_dat_ls$soil$clay,
                        xi = xi_df)
    
    c_stocks <- getC(model) # calc stocks for each pool per month
  }
  
  # run the simulations
  for (lu_lab in params_ls$lu_labs) {

    # business as usual only for now
    res_ls[[lu_lab]]$bau$c_stocks <- 
      run_model_sim(lu_lab, lu_lab, scn_dat_ls, gis_dat_ls, res_ls, params_ls) 
  }
  
  
  ## Package data for export
  c_pool_labs <- c("DPM", "RPM", "BIO", "HUM", "IOM")
  out_dat <- tibble(year = seq(2022, 2098+11/12, by = 1/12))
  for (lu_lab in params_ls$lu_labs) {
    suppressWarnings(sim_dat_tib <- as_tibble(res_ls[[lu_lab]]$bau$c_stocks))
    names(sim_dat_tib) <- paste(lu_lab, c_pool_labs, sep = "_")
    # quickly plot results if desired
    # sim_dat_tib %>%
    #   pivot_longer(cols = everything(), names_to = "c_pool", values_to = "c_stock") %>%
    #   mutate(year = rep(seq(2022, 2098+11/12, by = 1/12), 5)) %>%
    #   ggplot(aes(x = year, y = c_stock, color = c_pool)) +
    #   geom_line()
    out_dat <- cbind(out_dat, sim_dat_tib)
  }

  write_csv(out_dat, file = paste(out_dat_path, "sim_data_", params_ls$gcm, "_", params_ls$cell_id, ".csv", sep = ""))
  
  print(paste("Cell #", params_ls$cell_id, " simulation completed", sep = ""))
}

# Sample function call
# rothc_by_cell(cell_id = 15782,
#               gcm = "mri-cgcm3",
#               gis_dat_path = "input_data/cell_gis_dat_filter_ag.xlsx",
#               gcm_dat_path = "input_data/mri-cgcm3_1950-2099_by_cell/",
#               lm_dat_path = "input_data/rothc_lm_in_dat.xlsx",
#               out_dat_path = "output_data/")


args <- commandArgs(trailingOnly=TRUE)
# First argument is cell_id
# Second argument is gcm
# Third argument is gis_dat_path
# Fourth argument is gcm_dat_path
# Fifth argument is lm_dat_path
# Sixth argument is out_dat_path
rothc_by_cell(args[1], args[2], args[3], args[4], args[5], args[6])


