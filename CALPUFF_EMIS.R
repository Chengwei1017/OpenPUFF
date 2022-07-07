#CALPUFF EMIS TOOL by Chengwei Lu lcw@cdaes.cn
#20220219
#https://www.webqc.org/molecular-weight-of-HF.html
#datestr: YYYY-MM-DD-HH

library(ggplot)
library(geosphere)

GRID_ORIGIN=c(101.68,26.38)
GRID_DIR=c("E","N")

emis_path="E:\\обть\\calpuff\\CALPUFF6\\CALPUFFEM\\"

write_tab=function(shstr,shfile,append=T){
  write.table(shstr, shfile, col.names = FALSE,row.names = FALSE, quote = FALSE,append=append)
}

JulDate=function(indate){
  tmpdate=as.numeric(as.Date(indate)-as.Date(paste(substr(indate,1,4),"-01-01",sep=""))+1)
  if(nchar(tmpdate)==2) tmpdate=paste0("0",tmpdate)
  if(nchar(tmpdate)==1) tmpdate=paste0("00",tmpdate)
  return(paste0(substr(indate,1,4),tmpdate))
}

wstr=function(strin,slen=12){
  nlen=slen-nchar(strin)
  if(nlen<0)
    retstr=substr(strin,1,slen)
  else
    retstr=paste0(strin,paste0(rep(" ",nlen),collapse = ""))
  return(retstr)
}

time2str=function(datestr,endsec="0000"){
  strdate=substr(datestr,1,10)
  strhour=substr(datestr,12,13)
  jsday=JulDate(strdate)
  jsstr=paste0(
  substr(jsday,1,4)," ",
  substr(jsday,5,7)," ",
  strhour," ",endsec
  )
  return(jsstr)
}

gen_template=function(TEMPNAME,TYPE="AR"){
  headfile=paste0(emis_path,TEMPNAME,"_",TYPE,"_HEAD.csv")
  polifile=paste0(emis_path,TEMPNAME,"_",TYPE,"_POLI.csv")
  emisfile=paste0(emis_path,TEMPNAME,"_",TYPE,"_EMIS.csv")
  if(TYPE=="AR"){
    headdf=data.frame(SOURCENAME=NA,
                      UNIT=NA,
                      LAT1=NA,
                      LAT2=NA,
                      LAT3=NA,
                      LAT4=NA,
                      LON1=NA,
                      LON2=NA,
                      LON3=NA,
                      LON4=NA,
                      HT=NA,
                      ELEV=NA,
                      TEMPK=NA,
                      WEFF=NA,
                      REFF=NA,
                      SIGZ=NA)
    polidf=data.frame(POLNAME=NA,POLMWEIGHT=NA)
    emisdf=data.frame(STARTTIME  =NA,
                      ENDTIME    =NA,
                      SOURCENAME =NA,
                      POLNAME=NA)
  }
  write.csv(headdf,file=headfile)
  write.csv(polidf,file=polifile)
  write.csv(emisdf,file=emisfile)
  print(paste0("TEMPLATE FILES [",TEMPNAME,"] GENERATED."))
}

process_ar=function(emispack="TESTEM",emissday,emiseday){
  arhead=read.csv(paste0(emis_path,emispack,"_AR_HEAD.csv"),header=T)
  aremis=read.csv(paste0(emis_path,emispack,"_AR_EMIS.csv"),header=T)
  arpoli=read.csv(paste0(emis_path,emispack,"_AR_POLI.csv"),header=T)
  polstr=paste0(as.vector(sapply(paste0("'",as.vector(arpoli$POLNAME),"'"),wstr)),collapse = "")
  wgtstr=paste0(as.vector(arpoli$POLMWEIGHT),collapse = " ")
  nsrc=nrow(arhead)
  nse=nrow(arpoli)
  smdates=as.vector(unique(aremis$STARTTIME))
  emdates=as.vector(unique(aremis$ENDTIME))
  outfile=paste0(emispack,"_AR_CLEM.DAT")
  emfile=paste0(emis_path,outfile)

  
  headinfo=paste0(
    wstr("BAEMARB.DAT",16),
    wstr("2.1",16),
    wstr("Comments, times with seconds, time zone, coord info",64)
  )
  timeinfo=paste(
    time2str(emissday),time2str(emiseday,endsec="3600")
  )
  write_tab(headinfo,shfile=emfile,append=F)
  write_tab("1",shfile=emfile)
  write_tab(wstr("Emis file created by CET, CALPUFF EMIS TOOL.",80),shfile=emfile)
  write_tab(wstr("LCC",8),shfile=emfile)
  write_tab(paste0(wstr(paste0(GRID_ORIGIN[2],GRID_DIR[2]),16),wstr(paste0(GRID_ORIGIN[1],GRID_DIR[1]),16),wstr("30N",16),wstr("60N",16)),shfile=emfile)
  write_tab("0.0E+0 0.0E+0",shfile=emfile)
  write_tab(paste0(wstr("WGS-84",8),wstr("01-10-1996",12)),shfile=emfile)
  write_tab(wstr("KM",4),shfile=emfile)
  write_tab(wstr("UTC+0000",8),shfile=emfile)
  write_tab(timeinfo,shfile=emfile)
  write_tab(paste(nsrc,nse),shfile=emfile)
  write_tab(polstr,shfile=emfile)
  write_tab(wgtstr,shfile=emfile)
  for(xi in 1:nsrc){
    srcdf=arhead[xi,]
    write_tab(paste0(wstr(paste0("'",srcdf$SOURCENAME,"'"),16),
                     wstr(paste0("'",srcdf$UNIT,"'"),16),
                     " 0.0 0.0"),shfile=emfile)
  }
  emisloc=data.frame()
  for(xi in 1:length(smdates)){
    emdf=subset(aremis,STARTTIME==smdates[xi])
    if(smdates[xi]==-1){
      emistime=timeinfo
    }else{
      emistime=paste(
        time2str(smdates[xi]),time2str(emdates[xi],endsec="3600")
      )
    }
    
    
    
    write_tab(emistime,shfile=emfile)
    for(yi in 1:nrow(emdf)){
      iemis=aremis[yi,]
      ihead=subset(arhead,SOURCENAME==iemis$SOURCENAME)
      
      dist_x=distGeo(p1=data.frame(x=GRID_ORIGIN[1],y=GRID_ORIGIN[2]),
                        p2=data.frame(x=c(ihead$LON1,ihead$LON2,ihead$LON3,ihead$LON4),
                                y=c(GRID_ORIGIN[2],GRID_ORIGIN[2],GRID_ORIGIN[2],GRID_ORIGIN[2])))/1000
      
      dist_y=distGeo(p1=data.frame(x=GRID_ORIGIN[1],y=GRID_ORIGIN[2]),
                     p2=data.frame(x=c(GRID_ORIGIN[1],GRID_ORIGIN[1],GRID_ORIGIN[1],GRID_ORIGIN[1]),
                                   y=c(ihead$LAT1,ihead$LAT2,ihead$LAT3,ihead$LAT4)))/1000
      
      if(length(which(emisloc$SNAM==iemis$SOURCENAME))==0)
      emisloc=rbind(emisloc,data.frame(DISX=dist_x,DISY=dist_y,LATS=c(ihead$LAT1,ihead$LAT2,ihead$LAT3,ihead$LAT4),LONS=c(ihead$LON1,ihead$LON2,ihead$LON3,ihead$LON4),SNAM=iemis$SOURCENAME))
      
      write_tab(paste(wstr(paste0("'",iemis$SOURCENAME,"'"),16),
                      dist_x[1],
                      dist_x[2],
                      dist_x[3],
                      dist_x[4]),shfile=emfile)
      
      write_tab(paste(dist_y[1],
                      dist_y[2],
                      dist_y[3],
                      dist_y[4],
                      ihead$HT),shfile=emfile)
      
      write_tab(paste(ihead$ELEV,
                      ihead$TEMPK,
                      ihead$WEFF,
                      ihead$REFF,
                      ihead$SIGZ),shfile=emfile)
      
      write_tab(paste0(as.vector(t(iemis[,4:ncol(iemis)])),collapse = " "),shfile=emfile)
    }
  }
  
  mapfile=paste0(emis_path,emispack,"_AR_CLEM.MAP")
  write.csv(file=mapfile,emisloc)
  
  print(paste0("AR EMIS FILE [",outfile,"] GENERATED."))
}

process_ar("TESTEM","2021-01-01-00","2021-12-31-23")
