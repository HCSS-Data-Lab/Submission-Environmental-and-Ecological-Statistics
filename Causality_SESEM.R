#================= Information ===========================>

# Implementation of SESEM (Spatially Explicit Structural Equation Modeling) package.
# SESEM use paper: https://eprints.qut.edu.au/88787/1/88787.pdf
# SESEM package: https://cran.r-project.org/web/packages/sesem/sesem.pdf

# Date: 12/06/24
# Author: Rens van Dam

#================= Import libraries ======================>

library(sesem)
library(lavaan)
library(gplots)
library(ggplot2)

############################# Data preperation #####################################

#======================================

regions <- read.csv("iraq_regions_coordinatesadded.csv")
distances <- read.csv("network_distances.csv")
dist <- distances$distance
governorates <- regions$ADM1_EN
districts <- regions$ADM2_EN
x <- regions$ADM3_X
y <- regions$ADM3_Y

#=========================================

data <- read.csv("iraqmunicipalities2024.csv", stringsAsFactors = FALSE)

conflict_events <- data$ACLED_events_total
civilian_events <- data$ACLED_events_Violence_against_civilians
battle_events <- data$ACLED_events_Battles

#scaling the variables that are continuous
population_density <- scale(data$CIESIN_GPWv411_GPW_Population_Density_population_density_max)
skin_temperature <- scale(data$ECMWF_ERA5_LAND_MONTHLY_skin_temperature_max_std)
wheat_production <- scale(data$spam2020v2r0_global_prod_P_WHEA_A_mean)
latent_heatflux <- scale(data$NASA_FLDAS_NOAH01_C_GL_M_V001_Qle_tavg_max_std)
soil_water <- scale(data$ECMWF_ERA5_LAND_MONTHLY_volumetric_soil_water_layer_3_min_std)

datafile <- data.frame(x, y,
					conflict_events,
					civilian_events,
					battle_events,
					population_density,
					skin_temperature,
					wheat_production,
					latent_heatflux,
					soil_water)


############################# Bin creation #####################################

# SESEM takes spatially dependent variables into account by 
# performing the fit at different lag distances. Fundamentally,
# this means that you regress on variables at different distances
# from the dependent variable. We put the other datapoints in
# bins.
# Note: we only use datapoints of the lowest 20% of distances.
# Note: We have about 200 datapoints, so about 40000 distances.

districtdist <- dist

automatic_bins <- make.bin(districtdist, type="s.size", s.size = 500, p.dist = 20) #binsize 2000 ; can also try 
automatic_bins_size <- automatic_bins[1][[1]]
automatic_bins_name <- automatic_bins[2][[1]]

plotbin(districtdist, automatic_bins_size)

district_covar <- make.covar(datafile,districtdist,
                             automatic_bins_size,
                             automatic_bins_name)

############################# Spatial model specification #####################################

causal_effects <- 'total effect level' #every variable is given as a parameter.
causal_effects <- 'path level' #every path is given as a parameter to calculate its coefficient

if (causal_effects == 'total effect level') {

spatial_model <-'
	wheat_production ~ a*latent_heatflux + b*soil_water
	skin_temperature ~ c*latent_heatflux + d*wheat_production
	population_density ~ e*skin_temperature + f*soil_water
	civilian_events ~ g*population_density
	battle_events ~ h*civilian_events + i*skin_temperature
	conflict_events ~ j*battle_events + k*civilian_events

	tot_wheat_production := d*e*g*h*j + d*i*j + d*e*g*k
	tot_latent_heatflux := c*e*g*h*j + a*d*e*g*h*j + c*i*j + a*d*i*j + c*e*g*k + a*d*e*g*k
	tot_soil_water := b*d*e*g*h*j + f*g*h*j + b*d*i*j + b*d*e*g*k + f*g*k
	tot_skin_temperature := e*g*h*j + i*j + e*g*k
	tot_population_density := g*h*j + g*k
	tot_civilian_events := h*j + k
	tot_battle_events := j
	'
}

if (causal_effects == 'path level') {

spatial_model <-'
	wheat_production ~ a*latent_heatflux + b*soil_water
	skin_temperature ~ c*latent_heatflux + d*wheat_production
	population_density ~ e*skin_temperature + f*soil_water
	civilian_events ~ g*population_density
	battle_events ~ h*civilian_events + i*skin_temperature
	conflict_events ~ j*battle_events + k*civilian_events

	tot_wheat_production_1 := d * e * g * h * j
	tot_wheat_production_2 := d * i * j
	tot_wheat_production_3 := d * e * g * k
	tot_latent_heatflux_1 := c * e * g * h * j
	tot_latent_heatflux_2 := a * d * e * g * h * j
	tot_latent_heatflux_3 := c * i * j
	tot_latent_heatflux_4 := a * d * i * j
	tot_latent_heatflux_5 := c * e * g * k
	tot_latent_heatflux_6 := a * d * e * g * k
	tot_soil_water_1 := b * d * e * g * h * j
	tot_soil_water_2 := f * g * h * j
	tot_soil_water_3 := b * d * i * j
	tot_soil_water_4 := b * d * e * g * k
	tot_soil_water_5 := f * g * k
	tot_skin_temperature_1 := e * g * h * j
	tot_skin_temperature_2 := i * j
	tot_skin_temperature_3 := e * g * k
	tot_population_density_1 := g * h * j
	tot_population_density_2 := g * k
	tot_civilian_events_1 := h * j
	tot_civilian_events_2 := k
	tot_battle_events_1 := j
	'
}

############################# runModels function #####################################

runModels<-function(spatial_model,covdata){
	bin.summary<-covdata[[1]]
	model.fit<-matrix(numeric(1*length(bin.summary[,1])),ncol=length(bin.summary[,1]))
	names(model.fit)<-bin.summary[,1]
	par.name<-matrix(numeric(1*length(bin.summary[,1])),ncol=1)
	unst.est<-matrix(numeric(1*length(bin.summary[,1])),ncol=length(bin.summary[,1]))
	names(unst.est)<-bin.summary[,1]
	est.se<-matrix(numeric(1*length(bin.summary[,1])),ncol=length(bin.summary[,1]))
	names(est.se)<-bin.summary[,1]
	p.val<-matrix(numeric(1*length(bin.summary[,1])),ncol=length(bin.summary[,1]))
	names(p.val)<-bin.summary[,1]
	std.est<-matrix(numeric(1*length(bin.summary[,1])),ncol=length(bin.summary[,1]))
	names(std.est)<-bin.summary[,1]
	varname.r.square<-matrix(numeric(1*length(bin.summary[,1])),ncol=1)
	r.square<-matrix(numeric(1*length(bin.summary[,1])),ncol=length(bin.summary[,1]))
	names(r.square)<-bin.summary[,1]
	mod.parname<-matrix(numeric(1*length(bin.summary[,1])),ncol=1)
	mod.indices<-matrix(numeric(1*length(bin.summary[,1])),ncol=length(bin.summary[,1]))
	names(mod.indices)<-bin.summary[,1]
	for (i in 2:length(bin.summary[,1])-1) {
		bin_name<-bin.summary[i,1]
		print(factor(bin_name))### need to make this cleaner
		covmatrix<-data.frame(covdata[[4]][,,i])
		names(covmatrix)<-covdata[[2]]
		covmatrix<-as.matrix(round(covmatrix,8))#rounding needed to make matrix symmetrical 
												#difference in number of digits saved between 
												#upper and lower; as.matrix needed to convert 
												#dataframe to matrix for sem
		semmodel<-sem(spatial_model,sample.cov=covmatrix,sample.nobs=bin.summary[i,4])	
		if(inspect(semmodel,"converged")==FALSE) {print("Warning: convergence failure")}### figure out how to get to print w/o quotes
		options(warn=1)# print warnings when they occur
		if(inspect(semmodel,"converged")==TRUE) {
		model.fit[i]<-data.frame(fitMeasures(semmodel))
			unst.est[i]<-data.frame(parameterEstimates(semmodel)$est)
			est.se[i]<-data.frame(parameterEstimates(semmodel)$se)
			p.val[i]<-data.frame(parameterEstimates(semmodel)$pvalue)
			std.est[i]<-data.frame(parameterEstimates(semmodel,standardized=T)$std.all)
			if(i==1) {varname.r.square<-(rownames(data.frame((inspect(semmodel,"rsquare")))))}
			r.square[i]<-data.frame(inspect(semmodel,"rsquare"))
			#mod.indices[i]<-data.frame(modificationIndices(semmodel)$mi) 
			
			######### mod indices may need to be disabled for convergence issues
}}
		covmatrix<-data.frame(covdata[[3]])#nonspatial (flat) matrix
		names(covmatrix)<-covdata[[2]]
		covmatrix<-as.matrix(round(covmatrix,8)) 
		semmodel<-sem(spatial_model,sample.cov=covmatrix,sample.nobs=bin.summary[length(bin.summary[,1]),4])
		par.name<-data.frame(seq(1:length(parameterEstimates(semmodel)[,1])),paste(parameterEstimates(semmodel)[,"lhs"],parameterEstimates(semmodel)[,"op"],parameterEstimates(semmodel)[,"rhs"]))
		colnames(par.name)<-c("parameter.number","parameter.name")
		model.fit[length(bin.summary[,1])]<-data.frame(fitMeasures(semmodel))
		unst.est[length(bin.summary[,1])]<-data.frame(parameterEstimates(semmodel)[,"est"])
		#print(unst.est)
		est.se[length(bin.summary[,1])]<-data.frame(parameterEstimates(semmodel)[,"se"])
		p.val[length(bin.summary[,1])]<-data.frame(parameterEstimates(semmodel)[,"pvalue"])
		std.est[length(bin.summary[,1])]<-data.frame(parameterEstimates(semmodel,standardized=T)["std.all"])
		r.square[length(bin.summary[,1])]<-data.frame(inspect(semmodel,"rsquare"))
		mod.parname<-data.frame(paste(modificationIndices(semmodel)[,"lhs"],modificationIndices(semmodel)[,"op"],modificationIndices(semmodel)[,"rhs"]))
		colnames(mod.parname)<-c("parameter name")
		mod.indices[length(bin.summary[,1])]<-data.frame(modificationIndices(semmodel)[,"mi"])
	return(list(
	model_fit<-data.frame(model.fit,
	row.names=names(fitMeasures(semmodel))), # dataframe containing model fit indices for all bins
	par_name<-data.frame(par.name), # parameter names in rows
	unst_est<-data.frame(unst.est), # unstandardized parameter estimates in rows, column for each lag distance bin
	est_se<-data.frame(est.se), # std error of unstandardized parameter estimates in rows, column for each lag distance bin
	p_val<-data.frame(p.val), # p-value for each unstandardized parameter estimates in rows, column for each lag distance bin
	std_est<-data.frame(std.est),# standardized parameter estimates in rows, column for each lag distance bin
	varname_r_square<-data.frame(varname.r.square),#names of dependent variables
	r_square<-data.frame(r.square),#returns r2 values for dependent observed variables
	mod_parname<-data.frame(mod.parname),# list of all parameters for which there is a modification index in rows, column for each lag distance bin
	print(mod.indices),
	## MOD INDICES GIVES ERROR BECAUSE IT CAN NOT MAKE DATAFRAME WITH COLUMNS OF DIFFERENT ROW LENGTH ##
	#mod_indices<-data.frame(mod.indices),	# modification indices in rows, column for each lag distance bin
	bin.summary
	))
	options(warn=0) # return to default warning
}	
#print(colnames(district_covar))

spatial_model_results_districts <- runModels(spatial_model, district_covar)

################################### Results #########################################

print(bin.results(spatial_model_results_districts,bin="binflat") ) # prints results for flat (nonspatial) bin
print(bin.rsquare(spatial_model_results_districts,bin="binflat") ) # print Rsquared values for flat (nonspatial) bin
print(bin.results(spatial_model_results_districts,bin="Bin1") )# prints results for bin1 
print(bin.results(spatial_model_results_districts,bin="Bin2") )

#======================== Plotting ===========================>

# Plot of the model fit
plotmodelfit(spatial_model_results_districts,add.line="smooth",
             rmsea_err=F,
             pch=16,lty=1)

# Plot of the paths
if (causal_effects == 'total effect level') {
	plotpath(spatial_model_results_districts,path.type = 'user', selectpath = c(1,2,3,4,5,6,7,8,9,10,11), pch=11)
	plotpath(spatial_model_results_districts,path.type = 'user', selectpath = c(21,22,23,24,25,26,27), pch=11)
}
if (causal_effects == 'path level') {
	plotpath(spatial_model_results_districts,path.type = 'user', selectpath = c(1,2,3,4,5,6,7,8,9,10,11), pch=11)
	plotpath(spatial_model_results_districts,path.type = 'user', selectpath = c(21,22,23,24,25,26,27,28,29), pch=11)
	plotpath(spatial_model_results_districts,path.type = 'user', selectpath = c(30,31,32,33,34,35,36,37), pch=11)
	plotpath(spatial_model_results_districts,path.type = 'user', selectpath = c(38,39,40,41,42), pch=11)
}

modelfit = spatial_model_results_districts[[1]]
pathcoeffs_names = spatial_model_results_districts[[2]]
pathcoeffs = spatial_model_results_districts[[3]]
errorcoeffs = spatial_model_results_districts[[4]]
pvalcoeffs = spatial_model_results_districts[[5]]
rsquares = spatial_model_results_districts[[8]]
bins_info = spatial_model_results_districts[[11]]

print(summary)
write.csv(modelfit, "modelfit.csv", row.names = TRUE)
write.csv(pathcoeffs_names, "pathcoeffs_names.csv", row.names = TRUE)
write.csv(pathcoeffs, "pathcoeffs.csv", row.names = TRUE)
write.csv(errorcoeffs, "errorcoeffs.csv", row.names = TRUE)
write.csv(pvalcoeffs, "pvalcoeffs.csv", row.names = TRUE)
write.csv(rsquares, "rsquares.csv", row.names = TRUE)
write.csv(bins_info, "bins_info.csv", row.names = TRUE)