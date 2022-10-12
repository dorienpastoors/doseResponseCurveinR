

# get the fit for existing data with good curves
#plateReaderData_forSimFit <- plateReaderData %>% filter( timepoint == "t=24h" & compound == "compoundA")

temp_table <- read.table("plateReaderData_forSimFit.txt", stringsAsFactors = TRUE)

fit_list <- list()



for(k_cellType in levels(plateReaderData_forSimFit$cellType)){ 
  print(k_cellType)
  plateReaderData_temp <- plateReaderData_forSimFit %>% 
    filter(cellType == k_cellType & concentration_uM > 0)%>%
    mutate(y = value, x = log10(concentration_uM))%>%
    dplyr::select(x, y)
  #print(plateReaderData_temp)
  
  
  print(temp_result <- data.frame(cellType = k_cellType,
                                  IC50 = ED(drm_fit, 50)[[1]], error = ED(drm_fit, 50)[[2]]))
  
  drm_fit <- suppressMessages(drm( y ~x , fct = L.4(), data = plateReaderData_temp))
  
  
  
  fit_list <- c(fit_list, list(drm_fit))
  
}


# extract coefficients and play around
temp_simdata <- rdrm(3, fct = L.4(), coef(fit_list[[1]]), xerror = sort(unique(plateReaderData$concentration_uM)) , xpar = 1, yerror = "rnorm", ypar = c(0, 0.01), 
                     onlyY = TRUE)

simdata_df <- as.data.frame(t(temp_simdata[["y"]]))

simdata_df$conc <- sort(unique(plateReaderData$concentration_uM))

colnames(simdata_df)[1:3] <- c("rep1", "rep2", "rep3")

simdata_long <- simdata_df %>% pivot_longer(cols = c("rep1", "rep2", "rep3") )


ggplot(data = simdata_long, aes(x = conc, y = value))+
  ggtitle(paste0(round(coef(fit_list[[1]]), digits = 2), collapse = ";"))+
  geom_point()+
  scale_x_log10()


temp_coef1 <- coef(fit_list[[1]])

temp_coef2 <- coef(fit_list[[2]])

# get fewer concentrations
concentrations <- sort(unique(plateReaderData$concentration_uM))

concentrations <- concentrations[-2]

# function to get different coefficients and errors and write files

try_coeffs <- function(coef_mod, concentrations, ypar){
  
  temp_simdata <- rdrm(3, fct = L.4(), coef_mod, xerror = concentrations , xpar = 1, yerror = "rnorm", ypar = ypar, 
                       onlyY = TRUE)
  
  simdata_df <- as.data.frame(t(temp_simdata[["y"]]))
  
  simdata_df$conc <- concentrations
  
  colnames(simdata_df)[1:3] <- c("rep1", "rep2", "rep3")
  
  simdata_long <- simdata_df %>% pivot_longer(cols = c("rep1", "rep2", "rep3") )
  
  
  ggplot(data = simdata_long, aes(x = conc, y = value))+
    ggtitle(paste0(round(coef_mod, digits = 2), collapse = ";"))+
    geom_point()+
    scale_x_log10()
  
  ggsave(paste0("simdata/simDoseResp", paste0(round(coef_mod, digits = 2), collapse = "_"), ".png"))
  write.table(x = simdata_df %>%arrange(desc(conc))%>%dplyr::select(-conc), 
              file = paste0("simdata/simDoseResp", paste0(round(coef_mod, digits = 2), collapse = "_"), ".txt"), 
              quote = FALSE, row.names = FALSE, sep = "\t")
  
}


# generate sim data
temp_coef1a <- temp_coef1 

try_coeffs(temp_coef1a, concentrations = concentrations, ypar = c(0,0.02))


temp_coef1a[4] <- 0.1

try_coeffs(temp_coef1a, concentrations = concentrations, ypar = c(0,0.01))





temp_coef1a[4] <- 9

try_coeffs(temp_coef1a, concentrations = concentrations, ypar = c(0,0.03))


temp_coef1a[4] <- 20

try_coeffs(temp_coef1a, concentrations = concentrations, ypar = c(0,0.04))




temp_coef2a <- temp_coef2

temp_coef2a[4] <- 0.1

try_coeffs(temp_coef2a, concentrations = concentrations, ypar = c(0,0.01))

temp_coef2a[4] <- 0.2

try_coeffs(temp_coef2a, concentrations = concentrations, ypar = c(0,0.05))



temp_coef2a[4] <- 20

try_coeffs(temp_coef2a, concentrations = concentrations, ypar = c(0,0.02))

temp_coef2a[4] <- 2

try_coeffs(temp_coef2a, concentrations = concentrations, ypar = c(0,0.01))



# generate background
bckgr1 <- rnorm(n = 12, mean = 0.3, sd = 0.01)

write.table(bckgr1, "simdata/background_24hrs_rnorm_0.3.txt", quote = FALSE, row.names = FALSE, sep = "/t")

bckgr1 <- rnorm(n = 12, mean = 0.27, sd = 0.01)

write.table(bckgr1, "simdata/background_24hrs_rnorm_0.27.txt", quote = FALSE, row.names = FALSE, sep = "/t")
