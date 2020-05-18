#load('pw.RData')
library(maptools)
library(raster)
library(rgdal)
library(data.table)
library(sf)

################################################################
# 1. Create raster of worldpop:
#Take temperature raster as a template
temp <- raster("2001_Q/2001_Q1.nc", varname='t2m')
#Read in shapefile created in GEE and clean it (L1AD)
teow<-readOGR(".","worldpoppers_bueno_convec")
xxx <- teow[,-2]
#Establish extent and create raster frame
ext <-  extent(temp)

xy <- abs(apply(as.matrix(bbox(ext)), 1, diff))
n <- 10
r <- raster(ext, ncol=xy[1]*n, nrow=xy[2]*n)
#Rasterize the shp
rre <-rasterize(xxx, r)

################################################################
# 2. Clean temperature data
#list files (in this case raster nc) and stack quarters together
#grids <- list.files("./2001_Q/", pattern = "*.nc$")
grids=c("2001_Q1.nc", "2001_Q2.nc","2001_Q3.nc", "2001_Q4.nc")
s <- stack(paste0("./2001_Q/", grids))

################################################################
# 3. Create population weights for L1AD
BB <- read_sf("level1_gcs_modify3.shp")
out <- extract(rre, as_Spatial(BB))
df <- data.frame(SALID1=BB$SALID1, SUM=unlist(lapply(out, sum)))
result = mapply(FUN = `/`, lapply(out, function(x) { x } ), lapply(out, function(x) { sum(x) } ), SIMPLIFY = FALSE)

################################################################
# 4. Loop to multiple by daily temperature
k=list()
start_time <- Sys.time()
for(j in 1:365){
  #Extract daily raster to multiply it by the population weight
  out2 <- raster::extract(s[[j]], BB)
  result_f = mapply(FUN = `*`, result, out2, SIMPLIFY = FALSE)
  
  df2 <- data.frame(SALID1=BB$SALID1, SUM=unlist(lapply(result_f, sum))-273.15,
                    mean=unlist(lapply(out2, mean))-273.15,j)
  
  
  k[[j]] <- df2
}

end_time <- Sys.time()
end_time - start_time
df_total <-  (do.call(rbind,k))

values = seq(from = as.Date("2001-01-01"), to = as.Date("2001-12-31"), by = 'day')
vv=as.data.frame(rep(values, each=371))
colnames(vv)="date"
dff=cbind(df_total,vv)

colnames(dff)=c("SALID1","temp_pw","temp_x","time","date")

write.csv(dff,paste0("2001_complete/df_","2001",".csv"))
