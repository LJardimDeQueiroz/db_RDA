#Libraries:
library("vegan")
library("PCNM") 
library("vegan") 
library("ggplot2")
library("ggrepel")


#############################
#############################
#############################

#List of Files:

# (I)   "gen_trip.txt": pairwise-Fst matrix.
# (II)  "geo_trip.txt": geographical distance matrix
# (III) "veg_trip6_rawData.txt": floodplain vegetation composition
# (IV)  "Water color.txt": table including both water transparency and water pH
# (V)   "watefall.txt": geographical position of each sampling site #according to the Teotônio Falls, if up- or downstream
# (VI)  "floodplain.txt":floodplain size, expressed by a size index.


#Calling the variables:

# (I) Genetic distance
	
	gen <- read.table("gen_trip.txt", header=TRUE, fill=TRUE)
	gen_dist <- as.dist(gen)
	gen_matrix <- as.matrix(gen_dist)
	gen_dist <- as.dist(gen_matrix[order(rownames(gen_matrix)),
			order(colnames(gen_matrix))])
	

# (II) Geographical distance 
	
	geo <- read.table("geo_trip.txt", header=TRUE, fill=TRUE)
	geo_dist <- as.dist(geo)
	geo_dist <- as.matrix(geo_dist)
	geo_dist <- as.dist(geo_dist[order(rownames(geo_dist)),
			order(colnames(geo_dist))])
	geo_matrix <- as.matrix(geo_dist)

	# Applying a Principal Coordinates of Neighbourhood Matrix (PCNM)
	# following Borcard & Legendre 2002):
	 
		geo_pcnm <- PCNM(geo_dist)
		
		min <- geo_pcnm$thresh #truncation distance
		nb_ev_geo <- length(geo_pcnm$values) #number of eigenvalues

		geo_pcnm$expected_Moran
		geo_pcnm$Moran_I #Moran'I for each PCNM variables

		#Selecting only the significantly positive axis
		select_geo <- which(geo_pcnm$Moran_I$Positive == TRUE) 
		geo_pcnm_pos <-  as.data.frame(geo_pcnm$vectors) [, select_geo]

		#Renaming the colunms:
		colnames(geo_pcnm_pos) <- c("geo1", "geo2", "geo3", "geo4")
		

		#Plotting the output of PCNM in a two-dimensional space:
		
		geo_pcnm_pos$site  <- c("a1", "a2", "a3", "a4", "a5", "a6", "ara", "aru", 
		"b1", "cau", "ctl", "jac", "m1", "n1", "pur", "sam", 
		"slo", "sot", "t1")
			plot(geo_pcnm_pos$geo1, geo_pcnm_pos$geo4)
			ggplot(data=geo_pcnm_pos, aes(x = geo1, y=geo4)) +theme_bw() +
			geom_text_repel(aes(label = site),
			box.padding =unit(1, "lines")) +
			geom_point(shape =21, colour="black",
			fill ="green", size = 5, stroke =.1)


# (III) Floodplain vegetation composition
	
	veg <- read.table("veg_trip6_rawData.txt", header=TRUE)
	
	#Applying a Factor Analysis:
		veg_fact2 <- factanal(scale(veg), factors = 2, rotation="varimax", 
				scores=c("regression"))
		veg_facta <- veg_fact2$scores
	
		colnames(veg_facta) <- c("veg1","veg2")


# (IV) Water color: water transparency (mm) and pH

	color <- read.table("water_quality.txt", header=TRUE)
	
	
# (V) Waterfall: if 0, upstream of the TeotÙnio Falls; if 1 = downstream
	waterfall <- read.table("waterfall.txt", header=TRUE)

# (VI) Floddplain size
	floodplain <- read.table("floodplain.txt", header=TRUE)
	
	
#Merging the variables in a single file
		
	# geographical distance + floodplain vegetation composition
	trip <- merge (geo_pcnm_pos,veg_facta, by="row.names", all=TRUE)
	row.names(trip) <-  trip$Row.names
	trip <- trip[,-1]
	
	# + water color (water transparency and pH)
	trip <- merge(trip,color, by="row.names", all=TRUE)
	row.names(trip) <- trip$Row.names
	trip <- trip[,-1]
	
	# + waterfall 
	trip <- merge(trip, waterfall, by="row.names", all=TRUE)	
	row.names(trip) <- trip$Row.names
	trip <- trip[,-1]	
	
	# + floodplain size
	trip <- merge(trip, floodplain,  by="row.names", all=TRUE)
	row.names(trip) <- trip$Row.names
	trip <- trip[,-1]	
	
	fac_waterfall <- which(colnames(trip)=="waterfall")
	trip[fac_waterfall] <- factor(trip[[fac_waterfall]]) #waterfall as factor
	is.factor(trip$waterfall)


#################################

#NULL MODEL:

rdaNullModel <- capscale(gen_dist ~ 1, data=trip)


# COMPLETE MODEL:

rdaFullModel <- capscale(gen_dist ~ scale(geo1) + scale(geo2) + scale(geo3) + scale(geo4)
            + scale(transp) + scale(ph) + scale(floodplain) + scale(veg1) + scale(veg2)
            + waterfall,
            data=trip)


extractAIC(rdaFullModel)
vif.cca(rdaFullModel)
anova(rdaFullModel)

	#Calculating variance partitioning and Adjuste-R2 of the model


		VarpartFullModel <- varpart(gen_dist,
					 ~ scale(geo1)+ scale(geo2)+scale(geo3)+ scale(geo4),
                    ~ waterfall,
 					 ~ scale(transp) + scale(ph) + scale(floodplain) + scale(veg1) + scale(veg2)
					 , data=trip)

		plot(VarpartFullModel, digits = 2, cutoff = -Inf)


#FINDING THE BEST MODEL:

rda1 <- step(rdaNullModel, scope=formula(rdaFullModel),
	  test="perm", direction="both", steps=100)
anova.cca(rda1, by="mar")


#RUNNING BEST MODEL:
	
rdaBestModel <- capscale(gen_dist ~ scale(geo1)+scale(geo4)+waterfall + scale(transp) 
	    + scale(floodplain) + scale(veg1), 
	    data=trip, scale(TRUE))

summary(rdaBestModel)
extractAIC(rdaBestModel)
vif.cca(rdaBestModel)
anova(rdaBestModel)

anova(rdaBestModel, by="mar")
plot(rdaBestModel)

BestModelVarpart <- varpart(gen_dist, ~scale(geo1) + scale(geo4),
                        ~ waterfall,
                        ~ scale(transp) + scale(floodplain) + scale(veg1), data=trip)
plot(BestModelVarpart, digits = 2, cutoff = -Inf)

BestModelVarpartIBE <- varpart(gen_dist, ~ scale(floodplain),~ scale(veg1),~ scale(transp), data=trip)
plot(BestModelVarpartIBE, digits = 2, cutoff = -Inf)


#RUNNING COMPETITIVE MODELS ACCORDING TO THE
#OUTPUT OF "rda1"

#Competitive model 1
#Excluding "geo4" from the "best model":

rdaCompModel1<- capscale(gen_dist ~ scale(geo1) + waterfall + scale(transp) 
	    	+ scale(floodplain) + scale(veg1), 
	   	   data=trip)

anova(rdaCompModel1, by="mar", permu=10000)

extractAIC(rdaCompModel1)
vif.cca(rdaCompModel1)
plot(rdaCompModel1)

CompModel1Varpart = varpart(gen_dist, ~ scale(geo1), ~waterfall, ~ scale(transp)
			+ scale(floodplain) + scale(veg1),
			data=trip)

plot(CompModel1Varpart , digits = 2, cutoff = -Inf)


#Competitive model 2
#Excluding "geo4" and "waterfall" from the "best model":

rdaCompModel2<- capscale(gen_dist ~ scale(geo1)  + scale(transp) 
	    	+ scale(floodplain) + scale(veg1), 
	   	    data=trip)

anova(rdaCompModel2, by="mar", permu=10000)

extractAIC(rdaCompModel2)
vif.cca(rdaCompModel2)
plot(rdaCompModel2)

CompModel2Varpart = varpart(gen_dist, ~ scale(geo1), ~ scale(transp)
			+ scale(floodplain) + scale(veg1),
			data=trip)

plot(CompModel2Varpart , digits = 2, cutoff = -Inf)


#Competitive model 3
#Excluding "geo4", "waterfall" and "floodplain" from the best model:

rdaCompModel3<- capscale(gen_dist ~ scale(geo1)+ scale(transp) 
	    + scale(veg1), 
	    data=trip)

anova(rdaCompModel3, by="mar", permu=10000)

extractAIC(rdaCompModel3)
vif.cca(rdaCompModel3)


CompModel3Varpart = varpart(gen_dist, ~ scale(geo1),~ scale(transp) 
	    + scale(veg1), data=trip)

plot(CompModel3Varpart , digits = 2, cutoff = -Inf)


#Competitive model 4
#geo1+transp:

rdaCompModel4<- capscale(gen_dist ~ scale(geo1) + scale(transp), 
	    data=trip)

anova(rdaCompModel4, by="mar", permu=10000)

extractAIC(rdaCompModel4)
vif.cca(rdaCompModel4)

CompModel4Varpart = varpart(gen_dist, ~ scale(geo1), ~ scale(transp)
			, data=trip)

plot(CompModel4Varpart , digits = 2, cutoff = -Inf)


#Competitive model 5
#waterfall+transp:

rdaCompModel5<- capscale(gen_dist ~ scale(geo1),
	    data=trip)

anova(rdaCompModel5, by="mar", permu=10000)

extractAIC(rdaCompModel5)
vif.cca(rdaCompModel5)

CompModel5Varpart = varpart(gen_dist, ~ scale(geo1),
			, data=trip)

plot(CompModel5Varpart , digits = 2, cutoff = -Inf)


#Competitive model 6
#no geo1+geo4

rdaCompModel6<- capscale(gen_dist ~ scale(transp)+scale(floodplain)+scale(veg1)+waterfall,
	    data=trip)

anova(rdaCompModel6, by="mar", permu=10000)

extractAIC(rdaCompModel6)
vif.cca(rdaCompModel6)

CompModel6Varpart = varpart(gen_dist, ~ scale(transp)+scale(floodplain)+scale(veg1), ~waterfall
			, data=trip)

plot(CompModel6Varpart , digits = 2, cutoff = -Inf)
