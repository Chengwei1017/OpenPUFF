#PlotCALPUFF by Chengwei Lu lcw@cdaes.cn
#20210214
library(ggplot2)
library(metR)
library(geosphere)
library(dplyr)
library(foreach)
library(doParallel)
library(ggmap)

GRID_ORIGIN=c(101.68,26.38)
gdis=0.3

StartDate="2021-01-01"
EndDate="2021-01-31"

pol="pm10"

picsize=110
picscale=2
filetype="png"
picdpi=96

mcores=10

xrg=c(6,24)
yrg=c(6,24)
trimpic=T

withrec=c(101.814486,26.503327)
witharem="TESTEM"

TerrDir="E:\\下载\\calpuff\\CALPUFF6\\TERREL\\"
TerrelInfo=read.table(paste0(TerrDir,"TERR1KM.LST"),header=F,skip=94,nrows = 4)[,4:5]
TerrelInfo$V5=-1*TerrelInfo$V5
lats=TerrelInfo$V4
lons=TerrelInfo$V5

windowsFonts(HT=windowsFont("黑体"),
             RMN=windowsFont("Times New Roman"),
             ST=windowsFont("宋体"))

plottheme=theme(
  axis.text.x = element_text(family="RMN",size=9,angle=90),
  axis.text.y = element_text(family="RMN",size=9,angle=0),
  axis.title=element_text(family="RMN",size=11, colour="black",face="bold"),
  strip.text =element_text(family="RMN",size=9, colour="black"),
  strip.background = element_rect(fill = "gray90"),
  legend.title=element_blank(), #element_text(family="HT", size=size_t,angle=0,face="plain"),
  legend.text=element_text(family="RMN",size=11, colour="black",face="plain"),
  legend.position="top", 
  legend.direction="horizontal",
  legend.key.size=unit(0.5,"cm")
)

toRasterGrid=function(iGRID){
  tgrid=t(iGRID)
  ogrid=tgrid[nrow(tgrid):1,]
  return(ogrid)
}

toNCGrid=function(iGRID){
  tgrid=t(iGRID)
  ogrid=tgrid[,ncol(tgrid):1]
  return(ogrid)
}

gen_mapid=function(matr){
  matr=t(matr)
  x=ncol(matr)
  y=nrow(matr)
  df=data.frame()
  pb=txtProgressBar(max=x,style=3)
  for(xi in 1:x){
    for(yi in 1:y)
      df=rbind(df,data.frame(xdis=(xi-1)*gdis+gdis/2,ydis=(yi-1)*gdis+gdis/2))
    setTxtProgressBar(pb,xi)
  }
  close(pb)
  return(df)
}

format_calpuff=function(mapi,matr){
  res_for=data.frame(mapi,z=as.vector(as.matrix(matr)))
  return(res_for)
}

maxtrix2df_withdis=function(matr){
  matr=t(matr)
  x=ncol(matr)
  y=nrow(matr)
  df=data.frame()
  pb=txtProgressBar(max=x,style=3)
  for(xi in 1:x){
    for(yi in 1:y)
      df=rbind(df,data.frame(xdis=(xi-1)*gdis+gdis/2,ydis=(yi-1)*gdis+gdis/2,z=matr[xi,yi]))
    setTxtProgressBar(pb,xi)
  }
  close(pb)
  return(df)
}

uv2df_withdis=function(matu,matv){
  matu=t(matu)
  matv=t(matv)
  x=ncol(matu)
  y=nrow(matu)
  df=data.frame()
  pb=txtProgressBar(max=x,style=3)
  for(xi in 1:x){
    for(yi in 1:y)
      df=rbind(df,data.frame(xdis=(xi-1)*gdis+gdis/2,ydis=(yi-1)*gdis+gdis/2,u=matu[xi,yi],v=matv[xi,yi]))
    setTxtProgressBar(pb,xi)
  }
  close(pb)
  return(df)
}

metpath="E:\\下载\\calpuff\\CALPUFF6\\PRTMET\\"
conpath="E:\\下载\\calpuff\\CALPUFF6\\CALPOST\\"
pltpath="E:\\下载\\calpuff\\CALPUFF6\\PLOTS\\"
emspath="E:\\下载\\calpuff\\CALPUFF6\\CALPUFFEM\\"


dir.create(pltpath,recursive = T)

mapinfo=gen_mapid(read.table(file=paste0(metpath,"qaterr.grd"),header=F,skip=5))

grddata=format_calpuff(mapinfo,read.table(file=paste0(metpath,"qaterr.grd"),header=F,skip=5))
lucdata=format_calpuff(mapinfo,read.table(file=paste0(metpath,"qaluse.grd"),header=F,skip=5))

pos_ll=c(101.713385,26.750737)

get_pos_info=function(pos_ll){
  rec=data.frame(
    rec_x=distGeo(p1=GRID_ORIGIN,p2=c(pos_ll[1],GRID_ORIGIN[2]))/1000,
    rec_y=distGeo(p1=GRID_ORIGIN,p2=c(GRID_ORIGIN[1],pos_ll[2]))/1000
  )
  delmet=abs(rec$rec_x-grddata$xdis)+abs(rec$rec_y-grddata$ydis)
  extpmt=which(delmet==min(delmet))
  if(min(delmet)>2*gdis) print("Warning: receptor seems to be outside of the domain.")
  
  pos_info=data.frame(
    lat=pos_ll[2],
    lon=pos_ll[1],
    xds=rec$rec_x,
    yds=rec$rec_y,
    hgt=grddata[extpmt,]$z,
    luc=lucdata[extpmt,]$z
  )
  return(pos_info)
}

plot_lu=function(lucdata){
  p=ggplot()+geom_tile(aes(x=xdis,y=ydis,fill=as.character(z)),data=lucdata,color="gray30",size=0.1)+plottheme+theme_bw()
  #p=p+geom_text(aes(x=xdis,y=ydis,label=z),lucdata,size=2,color="white")
  p=p+scale_fill_manual("Land use",
                        labels = c("10" = "Urban",
                                   "20" = "Irrigated agriculture",
                                   "30" = "Range land",
                                   "40" = "Forest land",
                                   "50" = "Water",
                                   "51" = "Water",
                                   "60" = "Wetlands",
                                   "70" = "Barren land"),
                        values=c("red","orange","lightgreen","darkgreen","blue","blue","cyan","brown")
  )
  p=p+xlab("X Dis. to origin(Km)")+ylab("Y Dis. to origin(Km)")
  p=p+coord_fixed()
  if(!is.na(witharem)){
    mapar=read.csv(file=paste0(emspath,witharem,"_AR_CLEM.MAP"),header=T)
    p=p+geom_polygon(aes(x=DISX,y=DISY,group=SNAM),mapar,size=0.8,color="red",fill=NA)
  }
  
  if(!is.na(withrec[1])){
    rec=data.frame(
      rec_x=distGeo(p1=GRID_ORIGIN,p2=c(withrec[1],GRID_ORIGIN[2]))/1000,
      rec_y=distGeo(p1=GRID_ORIGIN,p2=c(GRID_ORIGIN[1],withrec[2]))/1000
    )
    p=p+geom_point(aes(x=rec_x,y=rec_y),rec,size=2,shape=17,color="purple")
  }
  p
}

plot_ht=function(grddata){
  p=ggplot()
  p=p+geom_contour_fill(aes(x=xdis,y=ydis,z=z),grddata,binwidth=50)+scale_fill_gradientn("Terrain",colours = terrain.colors(20))+plottheme+theme_bw()+coord_fixed()
  p=p+xlab("X Dis. to origin(Km)")+ylab("Y Dis. to origin(Km)")
  
  if(!is.na(witharem)){
    mapar=read.csv(file=paste0(emspath,witharem,"_AR_CLEM.MAP"),header=T)
    p=p+geom_polygon(aes(x=DISX,y=DISY,group=SNAM),mapar,size=0.8,color="red",fill=NA)
  }
  
  if(!is.na(withrec[1])){
    rec=data.frame(
      rec_x=distGeo(p1=GRID_ORIGIN,p2=c(withrec[1],GRID_ORIGIN[2]))/1000,
      rec_y=distGeo(p1=GRID_ORIGIN,p2=c(GRID_ORIGIN[1],withrec[2]))/1000
    )
    p=p+geom_point(aes(x=rec_x,y=rec_y),rec,size=2,shape=17,color="purple")
  }
  p
}

plot_map=function(){
  bbox=c(left=min(lons),bottom=min(lats),right=max(lons),top=max(lats))
  bmap=get_stamenmap(bbox = bbox,zoom = 12, maptype = "terrain")
  p=ggmap(bmap)+plottheme+theme_bw()
  if(!is.na(witharem)){
    mapar=read.csv(file=paste0(emspath,witharem,"_AR_CLEM.MAP"),header=T)
    p=p+geom_polygon(aes(x=LONS,y=LATS,group=SNAM),mapar,size=0.8,color="red",fill=NA)
  }
  
  if(!is.na(withrec[1])){
    p=p+geom_point(aes(x=withrec[1],y=withrec[2]),size=2,shape=17,color="purple")
  }
  
  p=p+xlab("Longitude")+ylab("Latitude")
  return(p)
}

plot_conc=function(myear,mmonth,mday,mhour,dhgt=50,maxconc=100,metonly=F,vec_skip=1){
  dconc=maxconc/10
  conc_breaks=MakeBreaks(10)
  conname=paste0(myear,"_m",mmonth,"_d",mday,"_",mhour,"00(utc+0000)_l00_",pol,"_1hr_conc.dat")
  metname=paste0(myear,"_M",mmonth,"_D",mday,"_",mhour,"00(UTC-0000)_L01_1HR")
  
  u_temp=format_calpuff(mapinfo,read.table(file=paste0(metpath,metname,".usp"),header=F,skip=5))
  v_temp=format_calpuff(mapinfo,read.table(file=paste0(metpath,metname,".vsp"),header=F,skip=5))
  uv=cbind(u_temp,v_temp$z)
  names(uv)=c("xdis","ydis","u","v")
  if(!file.exists(paste0(conpath,conname))){
    metonly=T
    print("No conc file found for this timestep, plot wind only.")
  }
  p=ggplot()
  if(!metonly){
    conc=read.table(file=paste0(conpath,conname),header=F,skip=6)
    names(conc)=c("xdis","ydis","conc")
    p=p+geom_contour_fill(aes(x=xdis,y=ydis,z=conc),data=conc,binwidth = dconc)
    p=p+scico::scale_fill_scico("Conc.",palette = "bilbao",breaks=seq(0,maxconc,dconc),limits=c(0,maxconc))+guides(fill = guide_colorsteps())
  }
  p=p+plottheme+theme_bw()
  p=p+geom_contour2(aes(x=xdis,y=ydis,z=z,color=..level..),grddata,size=0.5,show.legend = F,binwidth=dhgt)+scale_color_gradientn("Terrain",colours = terrain.colors(20))
  p=p+geom_vector(aes(x=xdis,y=ydis,dx=u,dy=v),data=uv,skip=vec_skip,size=0.3)
  p=p+geom_label_contour(aes(x=xdis,y=ydis,z=z,label=..level..),grddata,label.size=0.5,size=2.5)
  p=p+scale_mag("Wind Speed",max_size=0.8,max = 4)+coord_fixed()
  p=p+xlab("X Dis. to origin(Km)")+ylab("Y Dis. to origin(Km)")
  p=p+labs(caption=paste0(myear,"-",mmonth,"-",mday," ",mhour,":00"))
  
  if(!is.na(witharem)){
    mapar=read.csv(file=paste0(emspath,witharem,"_AR_CLEM.MAP"),header=T)
    p=p+geom_polygon(aes(x=DISX,y=DISY,group=SNAM),mapar,size=0.8,color="red",fill=NA)
  }
  
  if(!is.na(withrec[1])){
    rec=data.frame(
      rec_x=distGeo(p1=GRID_ORIGIN,p2=c(withrec[1],GRID_ORIGIN[2]))/1000,
      rec_y=distGeo(p1=GRID_ORIGIN,p2=c(GRID_ORIGIN[1],withrec[2]))/1000
    )
    p=p+geom_point(aes(x=rec_x,y=rec_y),rec,size=2,shape=17,color="purple")
  }

  return(p)
}


stat_grid=function(pol="PM10",ifMet=F,lev="0",fun=mean){
  hours=c(0:23)
  days=seq(as.Date(StartDate),as.Date(EndDate),1)
  concdf=data.frame()
  pb=txtProgressBar(max=length(days),style=3)
  rawtmp=mapinfo
  for(dayid in 1:length(days)){
    xday=days[dayid]
    daystr=paste0("",xday)
    iyea=substr(daystr,1,4)
    imon=substr(daystr,6,7)
    iday=substr(daystr,9,10)
    for(xhour in hours){
      if(xhour<10) ihour=paste0("0",xhour) else ihour=paste0("",xhour)
      if(xday==days[1] & xhour==0){
        #print("Skip the first hour of simulation...")
      }else{
        if(ifMet){
          conname=paste0(iyea,"_m",imon,"_d",iday,"_",ihour,"00(UTC-0000)_L0",lev,"_1HR.",pol)
          raw_data=format_calpuff(mapinfo,read.table(file=paste0(metpath,conname),header=F,skip=5))
          raw_info=mapinfo
        }else{
          conname=paste0(iyea,"_m",imon,"_d",iday,"_",ihour,"00(utc+0000)_l0",lev,"_",pol,"_1hr_conc.dat")
          raw_data=read.table(file=paste0(conpath,conname),header=F,skip=6)
          names(raw_data)=c("xdis","ydis","z")
          raw_info=raw_data[,c(1:2)]
        }
        rawtmp=cbind(rawtmp,raw_data$z)
      }
    }
    setTxtProgressBar(pb,dayid)
  }
  close(pb)
  rawtmp=rawtmp[,c(-1,-2)]
  rawsta=apply(rawtmp,1,fun)
  return(data.frame(raw_info,z=rawsta))
}

plot_conc_fun=function(fun=mean,dhgt=50,maxconc=18,metonly=F,vec_skip=1){
  print("Processing wind field...")
  dconc=maxconc/10
  u_temp=stat_grid(pol="usp",ifMet = T,lev=1,fun=mean)
  v_temp=stat_grid(pol="vsp",ifMet = T,lev=1,fun=mean)
  uv=cbind(u_temp,v_temp$z)
  names(uv)=c("xdis","ydis","u","v")

  p=ggplot()
  if(!metonly){
    print("Processing concentration...")
    conc=stat_grid(pol=pol,ifMet = F,lev=0,fun=fun)
    names(conc)=c("xdis","ydis","conc")
    p=p+geom_contour_fill(aes(x=xdis,y=ydis,z=conc),data=conc,binwidth = dconc)
    p=p+scico::scale_fill_scico("Conc.",palette = "bilbao",breaks=seq(0,maxconc,dconc),limits=c(0,maxconc))+guides(fill = guide_colorsteps())
  }
  p=p+plottheme+theme_bw()
  p=p+geom_contour2(aes(x=xdis,y=ydis,z=z,color=..level..),grddata,size=0.5,show.legend = F,binwidth=dhgt)+scale_color_gradientn("Terrain",colours = terrain.colors(20))
  p=p+geom_vector(aes(x=xdis,y=ydis,dx=u,dy=v),data=uv,skip=vec_skip,size=0.3)
  p=p+geom_label_contour(aes(x=xdis,y=ydis,z=z,label=..level..),grddata,label.size=0.5,size=2.5)
  p=p+scale_mag("Wind Speed",max_size=0.8,max = 4)+coord_fixed()
  p=p+xlab("X Dis. to origin(Km)")+ylab("Y Dis. to origin(Km)")
  p=p+labs(caption=paste0("Model period"))
  
  if(!is.na(witharem)){
    mapar=read.csv(file=paste0(emspath,witharem,"_AR_CLEM.MAP"),header=T)
    p=p+geom_polygon(aes(x=DISX,y=DISY,group=SNAM),mapar,size=0.8,color="red",fill=NA)
  }
  
  if(!is.na(withrec[1])){
    rec=data.frame(
      rec_x=distGeo(p1=GRID_ORIGIN,p2=c(withrec[1],GRID_ORIGIN[2]))/1000,
      rec_y=distGeo(p1=GRID_ORIGIN,p2=c(GRID_ORIGIN[1],withrec[2]))/1000
    )
    p=p+geom_point(aes(x=rec_x,y=rec_y),rec,size=2,shape=17,color="purple")
  }
  
  return(p)
}

ext_conc=function(myear,mmonth,mday,mhour){
  conname=paste0(myear,"_m",mmonth,"_d",mday,"_",mhour,"00(utc+0000)_l00_",pol,"_1hr_conc.dat")
  metname=paste0(myear,"_M",mmonth,"_D",mday,"_",mhour,"00(UTC-0000)_L01_1HR")
  
  u_temp=format_calpuff(mapinfo,read.table(file=paste0(metpath,metname,".usp"),header=F,skip=5))
  v_temp=format_calpuff(mapinfo,read.table(file=paste0(metpath,metname,".vsp"),header=F,skip=5))
  uv=cbind(u_temp,v_temp$z)
  names(uv)=c("xdis","ydis","u","v")
  metonly=F
  if(!file.exists(paste0(conpath,conname))){
    metonly=T
    print("No conc file found for this timestep, extrat met only.")
  }

  if(!metonly){
    conc=read.table(file=paste0(conpath,conname),header=F,skip=6)
    names(conc)=c("xdis","ydis","conc")
  }
 
  if(!is.na(withrec[1])){
    rec=data.frame(
      rec_x=distGeo(p1=GRID_ORIGIN,p2=c(withrec[1],GRID_ORIGIN[2]))/1000,
      rec_y=distGeo(p1=GRID_ORIGIN,p2=c(GRID_ORIGIN[1],withrec[2]))/1000
    )
    
    delmet=abs(rec$rec_x-grddata$xdis)+abs(rec$rec_y-grddata$ydis)
    extpmt=which(delmet==min(delmet))

    if(min(delmet)>2*gdis) print("Warning: receptor seems to be outside of the domain.")
    temp=format_calpuff(mapinfo,read.table(file=paste0(metpath,paste0(myear,"_M",mmonth,"_D",mday,"_",mhour,"00(UTC-0000)_L01_1HR"),".deg"),header=F,skip=5))
    prcp=format_calpuff(mapinfo,read.table(file=paste0(metpath,paste0(myear,"_M",mmonth,"_D",mday,"_",mhour,"00(UTC-0000)_L00_1HR"),".prc"),header=F,skip=5))
    mixh=format_calpuff(mapinfo,read.table(file=paste0(metpath,paste0(myear,"_M",mmonth,"_D",mday,"_",mhour,"00(UTC-0000)_L00_1HR"),".mix"),header=F,skip=5))
    
    if(metonly){
      icn=-1  
    }else{
      delcon=abs(rec$rec_x-conc$xdis)+abs(rec$rec_y-conc$ydis)
      extpco=which(delcon==min(delcon))
      icn=conc[extpco,]$conc
    }
    iwu=uv[extpmt,]$u
    iwv=uv[extpmt,]$v
    iws=sqrt(iwu^2+iwv^2)
    idr=atan2(iwu/iws, iwv/iws)*180/pi
    itm=temp[extpmt,]$z
    ipc=prcp[extpmt,]$z
    imh=mixh[extpmt,]$z
    extval=data.frame(TTAG=paste0(myear,"-",mmonth,"-",mday,"_",mhour),
                      WSPD=iws,
                      WDIR=idr,
                      TEMP=itm,
                      PRCP=ipc,
                      MIXH=imh,
                      CONC=icn)
  }else{
    extval=data.frame()
  }
  return(extval)
}

par_plot=function(dayid){
  xday=days[dayid]
  daystr=paste0("",xday)
  iyea=substr(daystr,1,4)
  imon=substr(daystr,6,7)
  iday=substr(daystr,9,10)
  for(xhour in hours){
    if(xday==days[1] & xhour==0){
      #print("Skip the first hour of simulation...")
    }else{
      if(xhour<10) ihour=paste0("0",xhour) else ihour=paste0("",xhour)
      plotfile=paste0(pltpath,"CALPUFF_CONC_",iyea,"_",imon,"_",iday,"_",ihour,".",filetype)
      print(paste0("Processing ",plotfile,"..."))
      plt=plot_conc(iyea,imon,iday,ihour)
      if(trimpic) plt=plt+xlim(xrg)+ylim(yrg)
      ggsave(filename=plotfile,height=picsize*0.9,width=picsize,unit="mm",scale=picscale,dpi=picdpi)
    }
  }
}

seq_conc=function(){
  hours=c(0:23)
  days=seq(as.Date(StartDate),as.Date(EndDate),1)
  concdf=data.frame()
  pb=txtProgressBar(max=length(days),style=3)
  for(dayid in 1:length(days)){
    xday=days[dayid]
    daystr=paste0("",xday)
    iyea=substr(daystr,1,4)
    imon=substr(daystr,6,7)
    iday=substr(daystr,9,10)
    for(xhour in hours){
      if(xday==days[1] & xhour==0){
        #print("Skip the first hour of simulation...")
      }else{
        if(xhour<10) ihour=paste0("0",xhour) else ihour=paste0("",xhour)
        conpak=ext_conc(iyea,imon,iday,ihour)
        concdf=rbind(concdf,conpak)
      }
    }
    setTxtProgressBar(pb,dayid)
  }
  close(pb)
  return(concdf)
}

find_max_conc=function(pol="PM10"){
  hours=c(0:23)
  days=seq(as.Date(StartDate),as.Date(EndDate),1)
  maxcon=data.frame()
  pb=txtProgressBar(max=length(days),style=3)
  for(dayid in 1:length(days)){
    xday=days[dayid]
    daystr=paste0("",xday)
    iyea=substr(daystr,1,4)
    imon=substr(daystr,6,7)
    iday=substr(daystr,9,10)
    for(xhour in hours){
      if(xday==days[1] & xhour==0){
        #print("Skip the first hour of simulation...")
      }else{
        if(xhour<10) ihour=paste0("0",xhour) else ihour=paste0("",xhour)
        conname=paste0(conpath,iyea,"_m",imon,"_d",iday,"_",ihour,"00(utc+0000)_l00_",pol,"_1hr_conc.dat")
        if(file.exists(conname)){
          conc=read.table(file=conname,header=F,skip=6)
          names(conc)=c("xdis","ydis","conc")
          sorted_conc=sort(conc$conc,decreasing = T)
          maxval=sorted_conc[1]
          max95=sorted_conc[length(sorted_conc)*0.05]
          max99=sorted_conc[length(sorted_conc)*0.01]
          maxcon=bind_rows(maxcon,data.frame(TIME=paste0(iyea,"_m",imon,"_d",iday,"_",ihour),max=maxval,max95=max95,max99=max99))
        }
      }
    }
    setTxtProgressBar(pb,dayid)
  }
  close(pb)
  return(maxcon)
}

maxconc=find_max_conc()
station_loc=c(102.1155,26.90011)
concext=seq_conc()
write.csv(file=paste0(pltpath,"CALPUFF_REC.csv"),concext)

hours=c(0:23)
days=seq(as.Date(StartDate),as.Date(EndDate),1)

cl=makeCluster(mcores)
registerDoParallel(cl)
foreach(dayid=1:length(days),.export=c("par_plot","plot_conc","xrg","yrg","days","withrec","trimpic","MakeBreaks","pol"),.packages=c("ggplot2","metR","geosphere"),.verbose = F) %dopar% {
  par_plot(dayid)
}
stopCluster(cl)

pltlu=plot_lu(lucdata)
pltht=plot_ht(grddata)
conc_avg=plot_conc_fun()
pltmp=plot_map()
conc_max=plot_conc_fun(fun=max,maxconc=max(maxconc$max99))
if(trimpic){
  pltlu=pltlu+xlim(xrg)+ylim(yrg)
  pltht=pltht+xlim(xrg)+ylim(yrg)
  conc_avg=conc_avg+xlim(xrg)+ylim(yrg)
  conc_max=conc_max+xlim(xrg)+ylim(yrg)
}
ggsave(pltmp,filename=paste0(pltpath,"CALPUFF_DOMAIN.",filetype),height=picsize*0.9,width=picsize,unit="mm",scale=picscale,dpi=picdpi)
ggsave(pltlu,filename=paste0(pltpath,"CALPUFF_LU.",filetype),height=picsize*0.9,width=picsize,unit="mm",scale=picscale,dpi=picdpi)
ggsave(pltht,filename=paste0(pltpath,"CALPUFF_HT.",filetype),height=picsize*0.9,width=picsize,unit="mm",scale=picscale,dpi=picdpi)
ggsave(conc_avg,filename=paste0(pltpath,"CALPUFF_AVG.",filetype),height=picsize*0.9,width=picsize,unit="mm",scale=picscale,dpi=picdpi)
ggsave(conc_max,filename=paste0(pltpath,"CALPUFF_MAX.",filetype),height=picsize*0.9,width=picsize,unit="mm",scale=picscale,dpi=picdpi)