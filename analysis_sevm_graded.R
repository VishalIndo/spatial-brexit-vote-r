## ----setup, include=TRUE, message=FALSE----------------------------------------------------
require(knitr)
require(sf)
require(spdep)
require(tmap)
require(ggplot2)
require(ggpubr)
require(dplyr)
require(corrplot)
require(MASS)
require(GGally)
require(ggeffects)
require(sjPlot)
require(spatialreg)
require(spgwr)




## ------------------------------------------------------------------------------------------

setwd("C:/Users/visha/Downloads")
brexit <- st_read(dsn="brexit.shp", layer = "brexit")
st_crs(brexit) <- "EPSG:27700"




## ------------------------------------------------------------------------------------------

# Calculate the surface point of the Brexit geometries
posurfBrexit <- st_point_on_surface(brexit)

# Create a nearest neighbor object for k=1
brexitK1nb <- st_coordinates(posurfBrexit) %>%
  knearneigh(k = 1) %>%
  knn2nb()

# Create a polygon neighbor object
brexitPolynb <- poly2nb(brexit)

# Union the polygon neighbor and nearest neighbor objects
brexitUnionedNb <- union.nb(brexitPolynb, brexitK1nb)

# Check if the neighbor object is symmetric
if (!is.symmetric.nb(brexitUnionedNb)) {
  brexitUnionedNb <- make.sym.nb(brexitUnionedNb)
}

# Calculate distances and create a weighted neighbor list
dlist <- nbdists(brexitUnionedNb, st_coordinates(posurfBrexit)) %>%
  lapply(function(x) 1/x)

brexitUnionedListw_d <- nb2listw(brexitUnionedNb, glist = dlist, style = "W")





## ------------------------------------------------------------------------------------------
# Fit a Moran's Eigenvector model using a quasipoisson family
meQp <- spatialreg::ME(
  Leave ~ prop_65_ov + prop_no_qu + prop_non_w, 
  family = quasipoisson, 
  data = brexit, 
  offset = log(Valid_Vote), 
  listw = brexitUnionedListw_d
)

# Display summary of the fitted Moran's Eigenvector model
summary(meQp)

# Fit a GLM model using quasipoisson family and including the fitted values from the Moran's Eigenvector model
model_SevmQp <- glm(
  Leave ~ prop_65_ov + prop_no_qu + prop_non_w + fitted(meQp), 
  family = quasipoisson, 
  data = brexit, 
  offset = log(Valid_Vote)
)

# Display summary of the fitted GLM model
summary(model_SevmQp)


## ------------------------------------------------------------------------------------------
# Fit a generalized linear model (GLM) using quasipoisson family
modelGlm <- glm(
  Leave ~ prop_65_ov + prop_no_qu + prop_non_w + offset(log(Valid_Vote)), 
  data = brexit, 
  family = quasipoisson
)

# Display the summary of the fitted GLM model
summary(modelGlm)





## ------------------------------------------------------------------------------------------
# Calculate the pseudo R-squared for the GLM model
pseudo_r2_glm <- 1 - modelGlm$deviance / modelGlm$null.deviance
# Print the pseudo R-squared for the GLM model
cat("Pseudo R-squared for the GLM model:", pseudo_r2_glm, "\n")

# Calculate the pseudo R-squared for the spatial error model (SEVM) using quasipoisson
pseudo_r2_sevmqp <- 1 - model_SevmQp$deviance / model_SevmQp$null.deviance
# Print the pseudo R-squared for the SEVM model
cat("Pseudo R-squared for the SEVM model:", pseudo_r2_sevmqp, "\n")




## ------------------------------------------------------------------------------------------
# Extract and print the coefficients of the GLM model
glm_coefficients <- coef(modelGlm)
cat("Coefficients of the GLM model:\n")
print(glm_coefficients)


## ------------------------------------------------------------------------------------------
# Extract and print the coefficients of the SEVM model
sevmqp_coefficients <- coef(model_SevmQp)
cat("Coefficients of the SEVM model:\n")
print(sevmqp_coefficients)



## ------------------------------------------------------------------------------------------
# Perform Moran's I test for spatial autocorrelation on the residuals of the SEVM model
moranSevmQp <- lm.morantest(model_SevmQp, listw = brexitUnionedListw_d)
cat("Moran's I test for SEVM model residuals:\n")
print(moranSevmQp)

# Perform Moran's I test for spatial autocorrelation on the residuals of the GLM model
moranGlmQp <- lm.morantest(modelGlm, listw = brexitUnionedListw_d)
cat("Moran's I test for GLM model residuals:\n")
print(moranGlmQp)

# Interpretation

cat("The Moran's I test assesses the presence of spatial autocorrelation in the model residuals.\n")
cat("A significant Moran's I in the GLM model residuals would suggest that spatial autocorrelation is present, indicating that the model has not fully accounted for spatial dependencies.\n")
cat("By comparing the Moran's I results from the SEVM and GLM models, we can determine if the SEVM model has successfully reduced spatial autocorrelation in the residuals.\n")
cat("If the Moran's I statistic is lower and less significant in the SEVM model compared to the GLM model, it indicates that the spatial error model has effectively accounted for spatial dependencies, reducing autocorrelation in the residuals.\n")





## ------------------------------------------------------------------------------------------
# Example function definition (replace with your actual function)
some_function <- function() {
  # Your function code here
  data.frame(matrix(rnorm(100), ncol=10))
}

# Generate sevQp using the function
sevQp <- some_function()
print(str(sevQp)) # Verify the structure



## ------------------------------------------------------------------------------------------
BrexitSvem <- st_read("brexit.shp")

# Verify that BrexitSvem is loaded correctly
print(st_geometry_type(BrexitSvem))
print(summary(BrexitSvem))


## ----1-------------------------------------------------------------------------------------
# Loop through each spatial eigenvector and plot separately

for (i in seq_along(colnames(sevQp))) {
  eigenvector_name <- colnames(sevQp)[i]
  
  # Create a map for each eigenvector
  map <- tm_shape(BrexitSvem) + 
    tm_polygons(col = eigenvector_name, palette = "-RdBu", lwd = 0.5,   
                n = 7, midpoint = 0, legend.show = FALSE) + 
    
    tm_layout(main.title = paste("Selected spatial eigenvector:", eigenvector_name), 
              legend.outside = TRUE, 
              panel.labels = eigenvector_name,
              panel.show = TRUE ) +
    tm_scale_bar()
  
  # Print the map
  print(map)
}



## ------------------------------------------------------------------------------------------
# Extract the coefficients, standard errors, and p-values from the model summary
theCoefficientDf <- summary(model_SevmQp)$coefficients %>% as.data.frame()

# Identify the rows corresponding to the fitted spatial eigenvectors
idx <- which(substr(row.names(theCoefficientDf), 1, 6) == "fitted")

# Extract the signs and coefficients of the spatial eigenvectors
theSign <- sign(theCoefficientDf[idx, 1])
theCoef <- theCoefficientDf[idx, 1]

# Identify the column indices of the spatial eigenvectors in the BrexitSvem data
idx <- which(substr(names(BrexitSvem), 1, 3) == "vec")
sevNames <- names(BrexitSvem)[idx]

# Initialize variables to store the new names for sign-corrected eigenvectors and their effects
sevNamesSignCor <- NULL
sevNamesEffect <- NULL
i <- 1

# Loop through each spatial eigenvector name
for(aName in sevNames) {
  newName <- paste0("sign_corr", aName)        # Create a new name for the sign-corrected eigenvector
  newNameEffect <- paste0("effect", aName)     # Create a new name for the eigenvector's effect
  
  # Apply sign correction to the eigenvector
  BrexitSvem[, newName] <- st_drop_geometry(BrexitSvem)[, aName] * theSign[i]
  
  # Calculate the effect of the coefficient at the response scale
  BrexitSvem[, newNameEffect] <- exp(st_drop_geometry(BrexitSvem)[, aName] * theCoef[i])
  
  # Update lists with the new names
  sevNamesSignCor <- c(sevNamesSignCor, newName)
  sevNamesEffect <- c(sevNamesEffect, newNameEffect)
  
  i <- i + 1
}

# The BrexitSvem object now contains sign-corrected eigenvectors and their effects.



## ------------------------------------------------------------------------------------------

sevNamesEffect <- sevNamesEffect[sevNamesEffect != "effectvec73"]
sevNamesSignCor <- sevNamesSignCor[sevNamesSignCor != "sign_corrvec73"]

# Now plot the maps for the remaining spatial eigenvectors
for (i in seq_along(sevNamesEffect)) {
  effect_name <- sevNamesEffect[i]
  
  # Create a map for each spatial eigenvector effect
  map <- tm_shape(BrexitSvem) + 
    tm_polygons(col = effect_name, palette = "-RdBu", lwd = 0.5,   
                n = 7, midpoint = 1, legend.show = TRUE, title = "Effect at\nresponse scale") + 
    tm_layout(main.title = paste("Spatial Eigenvector Effect:", effect_name), 
              legend.outside = TRUE,
              attr.outside = TRUE, panel.show = TRUE, 
              panel.labels = colnames(sevQp)[i], title.size = 0.7) +
    tm_scale_bar()
  
  # Print the map
  print(map)
}




## ------------------------------------------------------------------------------------------
# Calculate the product of the spatial eigenvector effects for each row
BrexitSvem$sevEffectResponseProd <- apply(
  st_drop_geometry(BrexitSvem[, sevNamesEffect]), 
  MARGIN = 1, 
  FUN = prod
)

# Create a map displaying the combined effect at the response scale
map <- tm_shape(BrexitSvem) + 
  tm_polygons(
    col = "sevEffectResponseProd", 
    palette = "-RdBu", 
    lwd = 0.5,   
    n = 7, 
    midpoint = 1,  
    legend.show = TRUE, 
    title = "Combined effect at\nresponse scale"
  ) + 
  tm_layout(
    main.title = "Selected spatial eigenvectors", 
    legend.outside = TRUE,
    attr.outside = TRUE, 
    panel.show = TRUE, 
    title.size = 0.7
  ) +
  tm_scale_bar()

# Print the map
print(map)





## ------------------------------------------------------------------------------------------
# Extract the coordinates from the spatial surface points of the Brexit data
brexit_coords <- st_coordinates(posurfBrexit)



## ------------------------------------------------------------------------------------------
brexit_sp <- as(brexit, "Spatial")


## ------------------------------------------------------------------------------------------
# Determine the optimal bandwidth for Geographically Weighted Regression (GWR) using AICc
optimal_bw <- bw.gwr(
  formula = Leave ~ prop_65_ov + prop_no_qu + prop_non_w, 
  data = brexit_sp, 
  approach = "AICc", 
  adapt = TRUE
)



## ------------------------------------------------------------------------------------------
# Perform Geographically Weighted Regression (GWR) using the optimal bandwidth
gwr_model <- gwr.basic(
  formula = Leave ~ prop_65_ov + prop_no_qu + prop_non_w, 
  data = brexit_sp,
  bw = optimal_bw,  # Use the previously determined optimal bandwidth
  adapt = TRUE
)

# Extract the regression coefficients for prop_no_qu from the GWR model
brexit$gwr_coeff_prop_no_qu <- gwr_model$SDF$prop_no_qu

# Display a summary of the GWR model
summary(gwr_model)


## ------------------------------------------------------------------------------------------
# Select the optimal bandwidth for Geographically Weighted Regression (GWR)
optimal_bandwidth <- gwr.sel(
  Leave ~ prop_65_ov + prop_no_qu + prop_non_w + offset(log(Valid_Vote)), 
  data = brexit, 
  coords = st_coordinates(posurfBrexit)
)

# Perform GWR using the selected bandwidth
gwr_model <- gwr(
  formula = Leave ~ prop_65_ov + prop_no_qu + prop_non_w + offset(log(Valid_Vote)), 
  data = brexit, 
  coords = st_coordinates(posurfBrexit), 
  bandwidth = optimal_bandwidth
)

# Add the spatially varying coefficient for prop_no_qu to the brexit data frame
brexit$gwr_coeff_prop_no_qu <- gwr_model$SDF@data$prop_no_qu

# Plot the spatially varying regression coefficient for prop_no_qu
tm_shape(brexit) + 
  tm_polygons(
    col = "gwr_coeff_prop_no_qu", 
    palette = "RdBu", 
    n = 7, 
    title = "Spatially Varying Coefficient\nfor prop_no_qu"
  ) + 
  tm_layout(
    main.title = "Spatially Varying Regression Coefficient", 
    legend.outside = TRUE
  ) +
  tm_scale_bar()



## ------------------------------------------------------------------------------------------
# Extract the GWR coefficients for the selected variables from the GWR model
gwr_coefficients <- as.data.frame(gwr_model$SDF@data[, c("prop_65_ov", "prop_no_qu", "prop_non_w")])

# Combine SEVM coefficients and GWR coefficients into a dataframe for comparison
coefficients_comparison <- data.frame(
  Variable = c("prop_65_ov", "prop_no_qu", "prop_non_w"),
  SEVM_Coefficients = coef(model_SevmQp)[2:4],  # Extracting SEVM coefficients for the variables
  GWR_Coefficients = sapply(gwr_coefficients, mean)  # Calculating the mean of GWR coefficients for each variable
)

# Display the comparison dataframe
print(coefficients_comparison)



## ------------------------------------------------------------------------------------------
# Deviance for SEVM model
deviance_sevm <- 1 - model_SevmQp$deviance / model_SevmQp$null.deviance

# Deviance for GWR model (requires manual calculation based on residuals)
gwr_residuals <- gwr_model$SDF@data$residuals
gwr_deviance <- 1 - (sum(gwr_residuals^2) / model_SevmQp$null.deviance)

# Compare the explained deviance
print(paste("Explained Deviance - SEVM Model:", deviance_sevm))
print(paste("Explained Deviance - GWR Model:", gwr_deviance))



## ------------------------------------------------------------------------------------------
# Add spatially varying coefficients to the spatial data frame
brexit$gwr_coeff_prop_no_qu <- gwr_model$SDF@data$prop_no_qu
brexit$gwr_coeff_prop_65_ov <- gwr_model$SDF@data$prop_65_ov
brexit$gwr_coeff_prop_non_w <- gwr_model$SDF@data$prop_non_w

# Map for prop_no_qu
map_prop_no_qu <- tm_shape(brexit) + 
  tm_polygons(
    col = "gwr_coeff_prop_no_qu", 
    palette = "RdBu", 
    n = 7, 
    midpoint = NA,  # Show the full spectrum of the color palette
    title = "Spatially Varying Coefficient\nfor prop_no_qu"
  ) + 
  tm_layout(
    main.title = "Spatially Varying Regression Coefficient for prop_no_qu", 
    legend.outside = TRUE,
    legend.outside.position = "right",
    legend.outside.size = 0.5
  ) +
  tm_scale_bar()

# Map for prop_65_ov
map_prop_65_ov <- tm_shape(brexit) + 
  tm_polygons(
    col = "gwr_coeff_prop_65_ov", 
    palette = "RdBu", 
    n = 7, 
    midpoint = NA,  # Show the full spectrum of the color palette
    title = "Spatially Varying Coefficient\nfor prop_65_ov"
  ) + 
  tm_layout(
    main.title = "Spatially Varying Regression Coefficient for prop_65_ov", 
    legend.outside = TRUE,
    legend.outside.position = "right",
    legend.outside.size = 0.5
  ) +
  tm_scale_bar()

# Map for prop_non_w
map_prop_non_w <- tm_shape(brexit) + 
  tm_polygons(
    col = "gwr_coeff_prop_non_w", 
    palette = "RdBu", 
    n = 7, 
    midpoint = NA,  # Show the full spectrum of the color palette
    title = "Spatially Varying Coefficient\nfor prop_non_w"
  ) + 
  tm_layout(
    main.title = "Spatially Varying Regression Coefficient for prop_non_w", 
    legend.outside = TRUE,
    legend.outside.position = "right",
    legend.outside.size = 0.5
  ) +
  tm_scale_bar()

# Print maps separately
print(map_prop_no_qu)
print(map_prop_65_ov)
print(map_prop_non_w)




## ------------------------------------------------------------------------------------------
# Map the residuals from SEVM and GLM models
tm_shape(brexit) + 
  tm_polygons(col = "gwr_coeff_prop_65_ov", palette = "-RdBu", lwd = 0.5, 
              n = 7, midpoint = 0, legend.show = TRUE, title = "GLM Residuals") + 
  tm_layout(main.title = "Residuals for GLM Model", 
            legend.outside = TRUE, 
            panel.show = TRUE) + 
  tm_scale_bar()

tm_shape(brexit) + 
  tm_polygons(col = "gwr_coeff_prop_non_w", palette = "-RdBu", lwd = 0.5, 
              n = 7, midpoint = 0, legend.show = TRUE, title = "SEVM Residuals") + 
  tm_layout(main.title = "Residuals for SEVM Model", 
            legend.outside = TRUE, 
            panel.show = TRUE) + 
  tm_scale_bar()

