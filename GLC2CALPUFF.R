#GLC2CALPUFF by Chengwei Lu lcw@cdaes.cn
#20210211
library(raster)

TerrDir="E:\\обть\\calpuff\\CALPUFF6\\TERREL\\"
CTGPDir="E:\\обть\\calpuff\\CALPUFF6\\CTGPROC\\"
DataDIR="E:\\обть\\calpuff\\CALPUFF6\\CALPUFFLU\\"
AimCenter=c(101.814486,26.503327)
Pharse_Terrel=T
if(Pharse_Terrel){
  TerrelInfo=read.table(paste0(TerrDir,"TERR1KM.LST"),header=F,skip=94,nrows = 4)[,4:5]
  TerrelInfo$V5=-1*TerrelInfo$V5
  lats=TerrelInfo$V4
  lons=TerrelInfo$V5
}else{
  lats=c(26.3492928,26.3500004,26.6695518,26.6688404)
  lons=c(102.004395,101.650002,101.650002,102.005966)
}

AimCenter[2]-mean(lats)
AimCenter[1]-mean(lons)

use_pzh=F
use_bld=T

glc2017=raster(paste0(DataDIR,"re_sc_glc_1.tif"),as.is=T)
if(!use_pzh){
  hgt250=raster(paste0(DataDIR,"SRTM_NE_250m.tif"),as.is=T)
}else{
  hgt250=raster(paste0(DataDIR,"PZHHRZ.tif"),as.is=T)
  hgt250=projectRaster(hgt250,crs=crs(glc2017))
}

if(use_bld){
  hgt250=raster(paste0(DataDIR,"CD_BLD.tif"),as.is=T)
}

luvalid=crop(glc2017,extent(min(lons)-0.1,max(lons)+0.1,min(lats)-0.1,max(lats)+0.1))
htvalid=crop(hgt250,extent(min(lons)-0.1,max(lons)+0.1,min(lats)-0.1,max(lats)+0.1))
#remapping to calpuff
values(luvalid)[values(luvalid)==1]=20
values(luvalid)[values(luvalid)==2]=30
values(luvalid)[values(luvalid)==3]=40
values(luvalid)[values(luvalid)==4]=40
values(luvalid)[values(luvalid)==5]=60
values(luvalid)[values(luvalid)==6]=50
values(luvalid)[values(luvalid)==7]=40
values(luvalid)[values(luvalid)==8]=10
values(luvalid)[values(luvalid)==9]=70

spplot(luvalid,col.regions=rainbow(10))
spplot(htvalid,col.regions = terrain.colors(9))


lulist=as.data.frame(as(luvalid, 'SpatialPointsDataFrame'))
names(lulist)=c("VALUE","LON","LAT")
callu=data.frame(ID=1:nrow(lulist),CAT=lulist$VALUE,LAT=lulist$LAT,LON=lulist$LON)

write.table(file=paste0(CTGPDir,"CALPUFF.LU"),callu,sep=",",col.names = F,row.names = F)
writeRaster(htvalid,filename=paste0(TerrDir,"CALPUFFHT.tif"),overwrite=T,options="COMPRESS=NA")
