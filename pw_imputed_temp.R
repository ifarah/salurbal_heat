setwd('~/Downloads')

library(sp)
library(raster)
library(rgdal)
library(velox)

options(scipen=999)

pop <- read_sf("worldpoppers_bueno_convec.shp")
#Read in files
temp <- raster("2001_Q/2001_Q1.nc", varname='t2m')

library(rgdal)
teow<-readOGR(".","worldpoppers_bueno_convec")

ext <-  extent(temp)
xy <- abs(apply(as.matrix(bbox(ext)), 1, diff))
n <- 10
r <- raster(ext, ncol=xy[1]*n, nrow=xy[2]*n)
xxx <- teow[,-2]


rre <-rasterize(xxx, r)
plot(rre)

grain_poly = rasterToPolygons(t) %>% 
  st_as_sf()

plot(grain_poly)

levels(rre)

?raster::levels

writeRaster(rre, "worldpop9x9.tif", format="GTIFF")

writeRaster(rre, filename="worldpop9x9.tif", format="GTiff", overwrite=TRUE)
wp9=raster("worldpop9x9.tif")



cmask <- raster("L1AD_mask_SALID.tif")






#BUENAZO

library(sf)
#Read in files
pop <- read_sf("worldpoppers_bueno_convec.shp")
temp <- raster::brick("2001_Q/2001_Q1.nc", varname='t2m')
BB <- read_sf("level1_gcs_modify3.shp")
AD <- readOGR("cdmxAD.shp")

boundary=st_transform(st_as_sf(AD),4326)
t=temp-273.15
grain_poly = rasterToPolygons(t) %>% 
  st_as_sf()

#Loop empezaria aqui
t_cr <- st_overlay(grain_poly, boundary)

pop2=st_transform(pop,4326)
p<- st_intersection(pop2, boundary)

p$w=p$sum/sum(p$sum)

t_cr$ID=1:nrow(t_cr)
p$ID=1:nrow(p)

pp=as.data.frame(cbind(p$ID, p$w))
colnames(pp)=c("ID","w")

m=merge(pp,t_cr, on='ID')

x_temp=mean(m$X2.metre.temperature)
pw_temp=sum(m$X2.metre.temperature*m$w)

x

setwd('~/Downloads')
#import required libraries
library(maptools)
library(raster)
library(rgdal)
library(data.table)
library(sf)

# 1. Create raster of worldpop:
#Take temperature raster as a template
temp <- raster("2001_Q/2001_Q1.nc", varname='t2m')
#Read in shapefile created in GEE and clean it
teow<-readOGR(".","worldpoppers_bueno_convec")
xxx <- teow[,-2]

#Establish extent and create raster frame
ext <-  extent(temp)

xy <- abs(apply(as.matrix(bbox(ext)), 1, diff))
n <- 10
r <- raster(ext, ncol=xy[1]*n, nrow=xy[2]*n)
#Rasterize the shp
rre <-rasterize(xxx, r)

# 2. Clean temperature data
#list files (in this case raster nc) and stack quarters together
#grids <- list.files("./2001_Q/", pattern = "*.nc$")
grids=c("2001_Q1.nc", "2001_Q2.nc","2001_Q3.nc", "2001_Q4.nc")


s <- stack(paste0("./2001_Q/", grids))

# 3. read-in the polygon shapefile
#poly <- readShapePoly("cdmxAD.shp")
#BB <- readShapePoly("level1_gcs_modify3.shp")
#ok=BB@data$SALID1[1:]

#st=brick("2001_Q/2001_Q1.nc", varname='t2m')
ss=s-273.15
BB <- read_sf("level1_gcs_modify3.shp")

rm(list= ls()[!(ls() %in% c('rre', 'ss', 'BB'))])

save.image(file='pw.RData')
rm(list = ls())
.rs.restartR()
load('pw.RData')
#Remove items from RAM
rm(list= ls()[!(ls() %in% c('rre', 's', 'BB', 'ok2'))])
rm(u)
d=list()

ok=c(101101, 101103)
ok=tail(BB$SALID1,363)

ok2=head(ok,3)
ok2=split(ok, ceiling(seq_along(ok)/3))
ok2[1]
ok2[[4]]
p

rre
rr <- disaggregate(rre, 10)

rr

e <- extract(rre, BB, fun=sum, na.rm=T, df=T)
whut=lapply(rre, function(x) { x/e$layer } )









# 1. Create raster of worldpop:
#Take temperature raster as a template
temp <- raster("2001_Q/2001_Q1.nc", varname='t2m')
#Read in shapefile created in GEE and clean it
teow<-readOGR(".","worldpoppers_bueno_convec")
xxx <- teow[,-2]

#Establish extent and create raster frame
ext <-  extent(temp)

xy <- abs(apply(as.matrix(bbox(ext)), 1, diff))
n <- 10
r <- raster(ext, ncol=xy[1]*n, nrow=xy[2]*n)
#Rasterize the shp
rre <-rasterize(xxx, r)

# 2. Clean temperature data
#list files (in this case raster nc) and stack quarters together
#grids <- list.files("./2001_Q/", pattern = "*.nc$")
grids=c("2001_Q1.nc", "2001_Q2.nc","2001_Q3.nc", "2001_Q4.nc")
s <- stack(paste0("./2001_Q/", grids))

# 3. Create population weights 
BB <- read_sf("level1_gcs_modify3.shp")
out <- extract(rre, BB)
df <- data.frame(SALID1=BB$SALID1, SUM=unlist(lapply(out, sum)))
result = mapply(FUN = `/`, lapply(out, function(x) { x } ), lapply(out, function(x) { sum(x) } ), SIMPLIFY = FALSE)

# 4. Loop to multiple by daily temperature
d=list()
start_time <- Sys.time()
for(j in 1:2){
  #Extract daily raster to multiply it by the population weight
  out2 <- extract(s[[j]], BB)
  result_f = mapply(FUN = `*`, result, out2, SIMPLIFY = FALSE)
  
  df2 <- data.frame(SALID1=BB$SALID1, SUM=unlist(lapply(result_f, sum))-273.15, j)
  
  
  d[[j]] <- df2
}

end_time <- Sys.time()
end_time - start_time


df_total <-  (do.call(rbind,d))


#for(j in BB$SALID1){
for(i in 8:10){
  load('pw.RData')
for(j in ok2[[i]]){
  #load('pw.RData')
  poly = BB[BB$SALID1 == j,]
  poly=as_Spatial(poly)
  
  # 4. Create population weights
  #Extract sum of population by city
  p <- extract(rre, poly, fun=sum, na.rm=T, df=T)
  #Create population weight by pixel
  pw<- calc(rre, function(x) x/p$layer) 
  
  #Check if sum=1 for each city
  #checking <- extract(pw, poly, fun=sum, na.rm=T, df=T)
  
  
  #Takes a while
  #Apply population weight to each temperature pixel
  tw=pw*s
  
  rm(pw)
  rm(p)
  
  #provides all hours for a year for ONE CITY
  #extract raster cell count of weighted temperatures (sum) within each polygon area (poly)
  #for (i in 1:length(grids)){
  #  final_tempw <- extract(tw, poly, fun=sum, na.rm=TRUE, df=TRUE)
  #}
  
  final_tempw <- extract(tw, poly, fun=sum, na.rm=TRUE, df=TRUE)
  rm(tw)
  u=transpose(final_tempw)
  rm(final_tempw)
  values = seq(from = as.Date("2001-01-01"), to = as.Date("2001-12-31"), by = 'day')
  df=cbind(as.data.frame((u[-1,])),values)
  df$ID=poly@data$SALID1
  colnames(df)=c("temp_pw","date","ID")
  df$temp_pw=df$temp_pw-273.15
  write.csv(df,paste0("2001/df_",j,".csv"))
  #d[[j]] <- df
  #rm(list = ls())
  #.rs.restartR()
}
  rm(list = ls())
  #.rs.restartR()
}
  
end_time <- Sys.time()
end_time - start_time
df_total <-  do.call(rbind,d)

write.csv(df_total, file = "./pw_temp_2001_Q1.csv")


#extract raster cell count (sum) within each polygon area (poly)
for (i in 1:length(rre)){
  p <- extract(rre, poly, fun=sum, na.rm=TRUE, df=TRUE)
}

#write to a data frame
df <- data.frame(ex)


#write to a CSV file
write.csv(df, file = "./path/to/data/CSV.csv")

x




















getValues(r)
getValues(s)


library(tidync)
rec_txa =
  # get the netcdf into a data frame
  tidync(
    here('analysis', 'data', 'record-ts', '2001_Q/2001_Q1.nc')) %>%
  activate('t2m') %>% .
hyper_tibble() 


r <- raster(ncol=36, nrow=18)
r[] <- runif(ncell(r))
cellStats(r, mean)
s <- r
s[] <- round(runif(ncell(r)) * 5)
zonal(r, s, 'sum')
## [1] 0.5179682


x
















wow=getValues(rre)
wow

t=temp-273.15
library(stars)

(nct = st_rasterize(pop["sum"], dx = 951, dy = 951))

library(fasterize)
install.packages("fasterize")
p1 <- rbind(c(-180,-20), c(-140,55), c(10, 0), c(-140,-60), c(-180,-20))
hole <- rbind(c(-150,-20), c(-100,-10), c(-110,20), c(-150,-20))
p1 <- list(p1, hole)
p2 <- list(rbind(c(-10,0), c(140,60), c(160,0), c(140,-55), c(-10,0)))
p3 <- list(rbind(c(-125,0), c(0,60), c(40,5), c(15,-45), c(-125,0)))
pols <- st_sf(value = rep(1,3),
              geometry = st_sfc(lapply(list(p1, p2, p3), st_polygon)))
r <- raster(pop, res = 100)
r <- fasterize(pop, r, field = "value", fun="sum")
plot(pols)


r <- raster(pop, res = 1)
r <- fasterize(pop, r, field = "value", fun="sum")


r <- fasterize(pols, r, field = "value", fun="sum")

r=st_rasterize(
  pop,
  template = st_as_stars(st_bbox(pop)),
  file = tempfile(),
  driver = "GTiff"
)


r

microfinance and democracy



plot(nct)

AD <- readOGR("cdmxAD.shp")

951, 951, 904401 
temp 
pop
#1.1 crop the temperature raster
dem_crop <- crop(temp, AD)
plot(dem_crop)
plot(AD,add=T)

as.matrix(pop$sum)

(x = st_as_binary(pop$sum))
rp <- rasterize(pop, r, 'AREA')

m=pop %>% 
  as.data.frame %>% 
  sf::st_as_sf(coords = c(1,2))

st_linestring(matrix(ncol=2,byrow=TRUE))

#1.2 crop the WorldPop raster
boundary=st_transform(st_as_sf(AD),4326)
pop2=st_transform(pop,4326)
ppl_crop <- st_intersection(pop2, boundary)
ppl_crop$percy=ppl_crop$sum/sum(ppl_crop$sum)

C= dem_crop- 273.15

vectorized=rasterToPolygons(C)
t=st_transform(st_as_sf(vectorized),4326)
C_bound <- st_intersection(t, boundary)
plot(C_bound)


temp_weighted=C_bound*ppl_crop$percy

temp_AD_hour=sum(as.vector(temp_weighted$X2.metre.temperature), na.rm=T)


ok=as.data.frame(colSums(Filter(is.numeric, temp_weighted)))

n<-dim(ok)[1]
df<-as.data.frame(ok[1:(n-4),])


plot(temp)







pop <- read_sf("worldpoppers_bueno_convec.shp")
#Read in files
temp <- brick("2001_Q/2001_Q1.nc", varname='t2m')
AD <- readOGR("level1_gcs_modify3.shp")

#1.1 crop the temperature raster
dem_crop <- crop(temp, AD)

#1.2 crop the WorldPop raster
boundary=st_transform(st_as_sf(AD),4326)
pop2=st_transform(pop,4326)
ppl_crop <- st_intersection(pop2, boundary)
ppl_crop$percy=ppl_crop$sum/sum(ppl_crop$sum)

C= dem_crop- 273.15

vectorized=rasterToPolygons(C)
t=st_transform(st_as_sf(vectorized),4326)
C_bound <- st_intersection(t, boundary)

plot(ppl_crop)
temp_weighted=C_bound*ppl_crop$percy

temp_AD_hour=sum(as.vector(temp_weighted$X2.metre.temperature), na.rm=T)

#foo = readOGR("worldpoppers_vec1.kml")


ok=as.data.frame(colSums(Filter(is.numeric, temp_weighted)))

n<-dim(ok)[1]
df<-as.data.frame(ok[1:(n-4),])
