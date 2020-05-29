#.rs.restartR()
#rm(list = ls())
setwd('~/Downloads')
load('pw_L1_guf.RData')

library(raster)
library(sf)
library(fasterize)
library(sp)
library(maptools)
library(rgdal)
library(data.table)
library(velox)

BB <- read_sf("level1_gcs_modify3.shp")
UX <- read_sf("level1_gcs_ux_modify4_20191127.shp")

# Find "results Yang" in Google Drive
#GUF
# BB_PANPERU <- BB%>%dplyr::filter(Country=="Peru"| Country=="Panama")
# UX_PANPERU <- UX%>%dplyr::filter(Country=="Peru"| Country=="Panama")
# gufy=sf::read_sf("guf_panperu_vectorized_complete.shp")
# temp <- raster("results_yang/2010_Q1.nc", varname='t2m')
# rt_guf <- fasterize(gufy, temp, field = "sum", fun="sum")

# out_guf <- (velox(rt_guf))$extract(sp=as_Spatial(BB_PANPERU), small=T)
# sum_out <- (velox(rt_guf))$extract(sp= as_Spatial(BB_PANPERU), small=T, fun = function(x) sum(x, na.rm = TRUE))
# result_guf = mapply(FUN = `/`, lapply(out_guf, function(x) { x } ),sum_out, SIMPLIFY = FALSE)

# out_guf_UX <- (velox(rt_guf))$extract(sp=as_Spatial(UX_PANPERU), small=T)
# sum_out_UX <- (velox(rt_guf))$extract(sp= as_Spatial(UX_PANPERU), small=T, fun = function(x) sum(x, na.rm = TRUE))
# result_guf_UX = mapply(FUN = `/`, lapply(out_guf_UX, function(x) { x } ),sum_out_UX, SIMPLIFY = FALSE)
# save.image(file='pw_L1_guf.RData') 

load('pw_L1_guf.RData')
each_year_L1_leap_GUF=function(year){
  grids=c((paste0(year,"_Q1.nc")), (paste0(year,"_Q2.nc")),(paste0(year,"_Q3.nc")), (paste0(year,"_Q4.nc")))
  s <- stack(paste0("./results_yang/", grids))
  
  #save.image(file='pw.RData')
  
  # 4. Loop to multiple by daily temperature
  k=list()
  #Change to 366 in leap years: 2004, 2008
  for(j in 1:366){
    #Extract daily raster to multiply it by the population weight
    #GUF
    out2_guf<- (velox(s[[j]]))$extract(sp = as_Spatial(BB_PANPERU), small=T)
    result_f_guf = mapply(FUN = `*`, result_guf, out2_guf, SIMPLIFY = FALSE)
    
    out2_guf_UX <- (velox(s[[j]]))$extract(sp = as_Spatial(UX_PANPERU), small=T)
    result_f_guf_UX = mapply(FUN = `*`, result_guf_UX, out2_guf_UX, SIMPLIFY = FALSE)
    
    df2_guf <- data.frame(SALID1=BB_PANPERU$SALID1, SUM=unlist(lapply(result_f_guf, sum, na.rm=TRUE))-273.15,
                          mean=unlist(lapply(out2_guf, mean, na.rm=TRUE))-273.15,
                          SUM_UX=unlist(lapply(result_f_guf_UX, sum, na.rm=T))-273.15,
                          mean_UX=unlist(lapply(out2_guf_UX, mean, na.rm=TRUE))-273.15, j)
    
    k[[j]] <- df2_guf
    
  }
  
  df_total_guf <-  (do.call(rbind,k))

  values = seq(from = as.Date(paste0(year,"-01-01")), to = as.Date(paste0(year,"-12-31")), by = 'day')
  vv_guf=as.data.frame(rep(values, each=26))
  colnames(vv_guf)="date"
  dff_guf=cbind(df_total_guf,vv_guf)
  
  colnames(dff_guf)=c("SALID1","ADtemp_gufw","ADtemp_x",
                      "UXtemp_gufw","UXtemp_x","time","date")
  
  dff_guf=dff_guf[-6]
  print(year)
  write.csv(dff_guf,paste0("final_results/L1/L1_",year,"_GUF.csv"))

}

start_time <- Sys.time()
each_year_L1_leap_GUF("2004")
each_year_L1_leap_GUF("2008")
each_year_L1_leap_GUF("2012")
end_time <- Sys.time()
print(end_time - start_time)

