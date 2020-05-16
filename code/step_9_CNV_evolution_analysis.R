####get the cytoband information of the cells####
library(data.table)

x <- fread("http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/cytoBand.txt.gz", 
           col.names = c("chrom","chromStart","chromEnd","name","gieStain"))

x[, .(length = sum(chromEnd - chromStart)), 
   by = .(chrom, arm = substring(name, 1, 1)) ]


x[, .(start=min(chromStart),end=max(chromEnd)), 
  by = .(chrom, arm = substring(name, 1, 1)) ]


cytoband<-x[, .(start=min(chromStart),end=max(chromEnd)), 
            by = .(chrom, arm = substring(name, 1, 1)) ]
cytoband

write.csv(cytoband,file="hg19_cytoband.csv")
getwd()
setwd("/Users/peizhanchen/Documents/hu_scRNA/NC_Revised/osteosarcoma_project/4_CNV/20200505")
list.files()



####分析BC2####
data<-read.table(file="BC2.dat.txt",header=T)
arm<-vector()
for (i in 1:nrow(data)){
  #site<-data[i,c("chr","start","end")]
  #cyto<-cytoband[site$chr,]
  site<-data[i,c("chr","start","end")]
  cyto<-cytoband[chrom==site$chr,]
  
  arm[i]<-ifelse(site$start > cyto[2,]$start,"q",
              ifelse(site$end < cyto[1,]$end,"p","other"))
  
}

gainorloss<-ifelse(data$state=="1","LOH",
                   ifelse(data$state=="2","LOH",
                          ifelse(data$state=="3","neutral",
                                 ifelse(data$state=="4","Gain",
                                        ifelse(data$state=="5","Gain",
                                               ifelse(data$state=="6","Gain","")
                                               )))))
data$cytob<-arm
data$gainorloss<-gainorloss

library(stringr)
chr<-str_replace(data$chr,pattern = "chr",replacement = "")
infor<-paste(chr,data$cytob,sep = "")
infor<-paste(infor,data$gainorloss,sep="_")
infor
data$infor<-infor

tmp<-data$cytob=="other"|data$gainorloss=="neutral"
table(!tmp)

data<-data[!tmp,]


cell_group<-unique(data$cell_group_name)
length(cell_group)
data3<-data[data$cell_group_name==cell_group[1],]
data3<-data3[!duplicated(data3$infor),]
data2<-data3


for (i in 2:length(cell_group)){
data3<-data[data$cell_group_name==cell_group[i],]
data3<-data3[!duplicated(data3$infor),]
data2<-rbind(data2,data3)
}

data3<-data.frame(data2[,c("cell_group_name","infor")])
data3$infor2<-1
data4<-reshape2::dcast(data3,cell_group_name~infor,value.var = "infor2")
BC2group<-read.table(file="BC2.cell_groupings.txt",header=T)
BC2group_2<-grepl(pattern="references",BC2group$cell_group_name)
BC2group<-BC2group[!BC2group_2,]
####combine the cells and the cytoband information###
merge<-dplyr::left_join(x = BC2group,y=data4,by=c("cell_group_name"))
merge[is.na(merge)]<-0
nrow(merge)
percent<-colSums(merge[,-c(1:2)])/nrow(merge)
dim(merge)
percent
write.csv(percent,file="BC2_percent.csv")
write.csv(data4,file="BC2_cell_canonical.csv")





####分析BC3####
data<-read.table(file="BC3.dat.txt",header=T)
arm<-vector()
for (i in 1:nrow(data)){
  #site<-data[i,c("chr","start","end")]
  #cyto<-cytoband[site$chr,]
  site<-data[i,c("chr","start","end")]
  cyto<-cytoband[chrom==site$chr,]
  
  arm[i]<-ifelse(site$start > cyto[2,]$start,"q",
                 ifelse(site$end < cyto[1,]$end,"p","other"))
  
}

gainorloss<-ifelse(data$state=="1","LOH",
                   ifelse(data$state=="2","LOH",
                          ifelse(data$state=="3","neutral",
                                 ifelse(data$state=="4","Gain",
                                        ifelse(data$state=="5","Gain",
                                               ifelse(data$state=="6","Gain","")
                                        )))))
data$cytob<-arm
data$gainorloss<-gainorloss

library(stringr)
chr<-str_replace(data$chr,pattern = "chr",replacement = "")
infor<-paste(chr,data$cytob,sep = "")
infor<-paste(infor,data$gainorloss,sep="_")
infor
data$infor<-infor

tmp<-data$cytob=="other"|data$gainorloss=="neutral"
table(!tmp)

data<-data[!tmp,]


cell_group<-unique(data$cell_group_name)
length(cell_group)
data3<-data[data$cell_group_name==cell_group[1],]
data3<-data3[!duplicated(data3$infor),]
data2<-data3


for (i in 2:length(cell_group)){
  data3<-data[data$cell_group_name==cell_group[i],]
  data3<-data3[!duplicated(data3$infor),]
  data2<-rbind(data2,data3)
}

data3<-data.frame(data2[,c("cell_group_name","infor")])
data3$infor2<-1
data4<-reshape2::dcast(data3,cell_group_name~infor,value.var = "infor2")


####read in the cell_grouping information###
BC3group<-read.table(file="BC3.cell_groupings.txt",header=T)
BC3group_2<-grepl(pattern="references",BC3group$cell_group_name)
BC3group<-BC3group[!BC3group_2,]
####combine the cells and the cytoband information###
merge<-dplyr::left_join(x = BC3group,y=data4,by=c("cell_group_name"))
merge[is.na(merge)]<-0
nrow(merge)
percent<-colSums(merge[,-c(1:2)])/nrow(merge)
dim(merge)
percent
write.csv(percent,file="BC3_percent.csv")
write.csv(data4,file="BC3_cell_canonical.csv")

####分析BC5####
data<-read.table(file="BC5.dat.txt",header=T)
arm<-vector()
for (i in 1:nrow(data)){
  #site<-data[i,c("chr","start","end")]
  #cyto<-cytoband[site$chr,]
  site<-data[i,c("chr","start","end")]
  cyto<-cytoband[chrom==site$chr,]
  
  arm[i]<-ifelse(site$start > cyto[2,]$start,"q",
                 ifelse(site$end < cyto[1,]$end,"p","other"))
  
}

gainorloss<-ifelse(data$state=="1","LOH",
                   ifelse(data$state=="2","LOH",
                          ifelse(data$state=="3","neutral",
                                 ifelse(data$state=="4","Gain",
                                        ifelse(data$state=="5","Gain",
                                               ifelse(data$state=="6","Gain","")
                                        )))))
data$cytob<-arm
data$gainorloss<-gainorloss

library(stringr)
chr<-str_replace(data$chr,pattern = "chr",replacement = "")
infor<-paste(chr,data$cytob,sep = "")
infor<-paste(infor,data$gainorloss,sep="_")
infor
data$infor<-infor

tmp<-data$cytob=="other"|data$gainorloss=="neutral"
table(!tmp)

data<-data[!tmp,]


cell_group<-unique(data$cell_group_name)
length(cell_group)
data3<-data[data$cell_group_name==cell_group[1],]
data3<-data3[!duplicated(data3$infor),]
data2<-data3


for (i in 2:length(cell_group)){
  data3<-data[data$cell_group_name==cell_group[i],]
  data3<-data3[!duplicated(data3$infor),]
  data2<-rbind(data2,data3)
}

data3<-data.frame(data2[,c("cell_group_name","infor")])
data3$infor2<-1
data4<-reshape2::dcast(data3,cell_group_name~infor,value.var = "infor2")


####read in the cell_grouping information###
BC5group<-read.table(file="BC5.cell_groupings.txt",header=T)
BC5group_2<-grepl(pattern="references",BC5group$cell_group_name)
BC5group<-BC5group[!BC5group_2,]
####combine the cells and the cytoband information###
merge<-dplyr::left_join(x = BC5group,y=data4,by=c("cell_group_name"))
merge[is.na(merge)]<-0
nrow(merge)
percent<-colSums(merge[,-c(1:2)])/nrow(merge)
dim(merge)
percent
write.csv(percent,file="BC5_percent.csv")
write.csv(data4,file="BC5_cell_canonical.csv")


####分析BC6####
data<-read.table(file="BC6.dat.txt",header=T)
arm<-vector()
for (i in 1:nrow(data)){
  #site<-data[i,c("chr","start","end")]
  #cyto<-cytoband[site$chr,]
  site<-data[i,c("chr","start","end")]
  cyto<-cytoband[chrom==site$chr,]
  
  arm[i]<-ifelse(site$start > cyto[2,]$start,"q",
                 ifelse(site$end < cyto[1,]$end,"p","other"))
  
}

gainorloss<-ifelse(data$state=="1","LOH",
                   ifelse(data$state=="2","LOH",
                          ifelse(data$state=="3","neutral",
                                 ifelse(data$state=="4","Gain",
                                        ifelse(data$state=="5","Gain",
                                               ifelse(data$state=="6","Gain","")
                                        )))))
data$cytob<-arm
data$gainorloss<-gainorloss

library(stringr)
chr<-str_replace(data$chr,pattern = "chr",replacement = "")
infor<-paste(chr,data$cytob,sep = "")
infor<-paste(infor,data$gainorloss,sep="_")
infor
data$infor<-infor

tmp<-data$cytob=="other"|data$gainorloss=="neutral"
table(!tmp)

data<-data[!tmp,]


cell_group<-unique(data$cell_group_name)
length(cell_group)
data3<-data[data$cell_group_name==cell_group[1],]
data3<-data3[!duplicated(data3$infor),]
data2<-data3


for (i in 2:length(cell_group)){
  data3<-data[data$cell_group_name==cell_group[i],]
  data3<-data3[!duplicated(data3$infor),]
  data2<-rbind(data2,data3)
}

data3<-data.frame(data2[,c("cell_group_name","infor")])
data3$infor2<-1
data4<-reshape2::dcast(data3,cell_group_name~infor,value.var = "infor2")


####read in the cell_grouping information###
BC6group<-read.table(file="BC6.cell_groupings.txt",header=T)
BC6group_2<-grepl(pattern="references",BC6group$cell_group_name)
BC6group<-BC6group[!BC6group_2,]
####combine the cells and the cytoband information###
merge<-dplyr::left_join(x = BC6group,y=data4,by=c("cell_group_name"))
merge[is.na(merge)]<-0
nrow(merge)
percent<-colSums(merge[,-c(1:2)])/nrow(merge)
dim(merge)
percent
write.csv(percent,file="BC6_percent.csv")
write.csv(data4,file="BC6_cell_canonical.csv")




####分析BC10####
data<-read.table(file="BC10.dat.txt",header=T)
arm<-vector()
for (i in 1:nrow(data)){
  #site<-data[i,c("chr","start","end")]
  #cyto<-cytoband[site$chr,]
  site<-data[i,c("chr","start","end")]
  cyto<-cytoband[chrom==site$chr,]
  
  arm[i]<-ifelse(site$start > cyto[2,]$start,"q",
                 ifelse(site$end < cyto[1,]$end,"p","other"))
  
}

gainorloss<-ifelse(data$state=="1","LOH",
                   ifelse(data$state=="2","LOH",
                          ifelse(data$state=="3","neutral",
                                 ifelse(data$state=="4","Gain",
                                        ifelse(data$state=="5","Gain",
                                               ifelse(data$state=="6","Gain","")
                                        )))))
data$cytob<-arm
data$gainorloss<-gainorloss

library(stringr)
chr<-str_replace(data$chr,pattern = "chr",replacement = "")
infor<-paste(chr,data$cytob,sep = "")
infor<-paste(infor,data$gainorloss,sep="_")
infor
data$infor<-infor

tmp<-data$cytob=="other"|data$gainorloss=="neutral"
table(!tmp)

data<-data[!tmp,]


cell_group<-unique(data$cell_group_name)
length(cell_group)
data3<-data[data$cell_group_name==cell_group[1],]
data3<-data3[!duplicated(data3$infor),]
data2<-data3


for (i in 2:length(cell_group)){
  data3<-data[data$cell_group_name==cell_group[i],]
  data3<-data3[!duplicated(data3$infor),]
  data2<-rbind(data2,data3)
}

data3<-data.frame(data2[,c("cell_group_name","infor")])
data3$infor2<-1
data4<-reshape2::dcast(data3,cell_group_name~infor,value.var = "infor2")


####read in the cell_grouping information###
BC10group<-read.table(file="BC10.cell_groupings.txt",header=T)
BC10group_2<-grepl(pattern="references",BC10group$cell_group_name)
BC10group<-BC10group[!BC10group_2,]
####combine the cells and the cytoband information###
merge<-dplyr::left_join(x = BC10group,y=data4,by=c("cell_group_name"))
merge[is.na(merge)]<-0
nrow(merge)
percent<-colSums(merge[,-c(1:2)])/nrow(merge)
dim(merge)
percent
write.csv(percent,file="BC10_percent.csv")
write.csv(data4,file="BC10_cell_canonical.csv")




####分析BC11####
data<-read.table(file="BC11.dat.txt",header=T)
arm<-vector()
for (i in 1:nrow(data)){
  #site<-data[i,c("chr","start","end")]
  #cyto<-cytoband[site$chr,]
  site<-data[i,c("chr","start","end")]
  cyto<-cytoband[chrom==site$chr,]
  
  arm[i]<-ifelse(site$start > cyto[2,]$start,"q",
                 ifelse(site$end < cyto[1,]$end,"p","other"))
  
}

gainorloss<-ifelse(data$state=="1","LOH",
                   ifelse(data$state=="2","LOH",
                          ifelse(data$state=="3","neutral",
                                 ifelse(data$state=="4","Gain",
                                        ifelse(data$state=="5","Gain",
                                               ifelse(data$state=="6","Gain","")
                                        )))))
data$cytob<-arm
data$gainorloss<-gainorloss

library(stringr)
chr<-str_replace(data$chr,pattern = "chr",replacement = "")
infor<-paste(chr,data$cytob,sep = "")
infor<-paste(infor,data$gainorloss,sep="_")
infor
data$infor<-infor

tmp<-data$cytob=="other"|data$gainorloss=="neutral"
table(!tmp)

data<-data[!tmp,]


cell_group<-unique(data$cell_group_name)
length(cell_group)
data3<-data[data$cell_group_name==cell_group[1],]
data3<-data3[!duplicated(data3$infor),]
data2<-data3


for (i in 2:length(cell_group)){
  data3<-data[data$cell_group_name==cell_group[i],]
  data3<-data3[!duplicated(data3$infor),]
  data2<-rbind(data2,data3)
}

data3<-data.frame(data2[,c("cell_group_name","infor")])
data3$infor2<-1
data4<-reshape2::dcast(data3,cell_group_name~infor,value.var = "infor2")


####read in the cell_grouping information###
BC11group<-read.table(file="BC11.cell_groupings.txt",header=T)
BC11group_2<-grepl(pattern="references",BC11group$cell_group_name)
BC11group<-BC11group[!BC11group_2,]
####combine the cells and the cytoband information###
merge<-dplyr::left_join(x = BC11group,y=data4,by=c("cell_group_name"))
merge[is.na(merge)]<-0
nrow(merge)
percent<-colSums(merge[,-c(1:2)])/nrow(merge)
dim(merge)
percent
write.csv(percent,file="BC11_percent.csv")
write.csv(data4,file="BC11_cell_canonical.csv")








####分析BC17####
data<-read.table(file="BC17.dat.txt",header=T)
arm<-vector()
for (i in 1:nrow(data)){
  #site<-data[i,c("chr","start","end")]
  #cyto<-cytoband[site$chr,]
  site<-data[i,c("chr","start","end")]
  cyto<-cytoband[chrom==site$chr,]
  
  arm[i]<-ifelse(site$start > cyto[2,]$start,"q",
                 ifelse(site$end < cyto[1,]$end,"p","other"))
  
}

gainorloss<-ifelse(data$state=="1","LOH",
                   ifelse(data$state=="2","LOH",
                          ifelse(data$state=="3","neutral",
                                 ifelse(data$state=="4","Gain",
                                        ifelse(data$state=="5","Gain",
                                               ifelse(data$state=="6","Gain","")
                                        )))))
data$cytob<-arm
data$gainorloss<-gainorloss

library(stringr)
chr<-str_replace(data$chr,pattern = "chr",replacement = "")
infor<-paste(chr,data$cytob,sep = "")
infor<-paste(infor,data$gainorloss,sep="_")
infor
data$infor<-infor

tmp<-data$cytob=="other"|data$gainorloss=="neutral"
table(!tmp)

data<-data[!tmp,]


cell_group<-unique(data$cell_group_name)
length(cell_group)
data3<-data[data$cell_group_name==cell_group[1],]
data3<-data3[!duplicated(data3$infor),]
data2<-data3


for (i in 2:length(cell_group)){
  data3<-data[data$cell_group_name==cell_group[i],]
  data3<-data3[!duplicated(data3$infor),]
  data2<-rbind(data2,data3)
}

data3<-data.frame(data2[,c("cell_group_name","infor")])
data3$infor2<-1
data4<-reshape2::dcast(data3,cell_group_name~infor,value.var = "infor2")

data4<-data4[-c(9,10),]

####read in the cell_grouping information###
BC17group<-read.table(file="BC17.cell_groupings.txt",header=T)
BC17group_2<-grepl(pattern="references",BC17group$cell_group_name)
BC17group<-BC17group[!BC17group_2,]
####combine the cells and the cytoband information###
merge<-dplyr::left_join(x = BC17group,y=data4,by=c("cell_group_name"))
merge[is.na(merge)]<-0
nrow(merge)
percent<-colSums(merge[,-c(1:2)])/nrow(merge)
dim(merge)
percent
write.csv(percent,file="BC17_percent.csv")
write.csv(data4,file="BC17_cell_canonical.csv")






####分析BC20####
data<-read.table(file="BC20.dat.txt",header=T)
arm<-vector()
for (i in 1:nrow(data)){
  #site<-data[i,c("chr","start","end")]
  #cyto<-cytoband[site$chr,]
  site<-data[i,c("chr","start","end")]
  cyto<-cytoband[chrom==site$chr,]
  
  arm[i]<-ifelse(site$start > cyto[2,]$start,"q",
                 ifelse(site$end < cyto[1,]$end,"p","other"))
  
}

gainorloss<-ifelse(data$state=="1","LOH",
                   ifelse(data$state=="2","LOH",
                          ifelse(data$state=="3","neutral",
                                 ifelse(data$state=="4","Gain",
                                        ifelse(data$state=="5","Gain",
                                               ifelse(data$state=="6","Gain","")
                                        )))))
data$cytob<-arm
data$gainorloss<-gainorloss

library(stringr)
chr<-str_replace(data$chr,pattern = "chr",replacement = "")
infor<-paste(chr,data$cytob,sep = "")
infor<-paste(infor,data$gainorloss,sep="_")
infor
data$infor<-infor

tmp<-data$cytob=="other"|data$gainorloss=="neutral"
table(!tmp)

data<-data[!tmp,]


cell_group<-unique(data$cell_group_name)
length(cell_group)
data3<-data[data$cell_group_name==cell_group[1],]
data3<-data3[!duplicated(data3$infor),]
data2<-data3


for (i in 2:length(cell_group)){
  data3<-data[data$cell_group_name==cell_group[i],]
  data3<-data3[!duplicated(data3$infor),]
  data2<-rbind(data2,data3)
}

data3<-data.frame(data2[,c("cell_group_name","infor")])
data3$infor2<-1
data4<-reshape2::dcast(data3,cell_group_name~infor,value.var = "infor2")

#data4<-data4[-c(9,10),]

####read in the cell_grouping information###
BC20group<-read.table(file="BC20.cell_groupings.txt",header=T)
BC20group_2<-grepl(pattern="references",BC20group$cell_group_name)
BC20group<-BC20group[!BC20group_2,]
####combine the cells and the cytoband information###
merge<-dplyr::left_join(x = BC20group,y=data4,by=c("cell_group_name"))
merge[is.na(merge)]<-0
nrow(merge)
percent<-colSums(merge[,-c(1:2)])/nrow(merge)
dim(merge)
percent
write.csv(percent,file="BC20_percent.csv")
write.csv(data4,file="BC20_cell_canonical.csv")



####分析BC21####
data<-read.table(file="BC21.dat.txt",header=T)
arm<-vector()
for (i in 1:nrow(data)){
  #site<-data[i,c("chr","start","end")]
  #cyto<-cytoband[site$chr,]
  site<-data[i,c("chr","start","end")]
  cyto<-cytoband[chrom==site$chr,]
  
  arm[i]<-ifelse(site$start > cyto[2,]$start,"q",
                 ifelse(site$end < cyto[1,]$end,"p","other"))
  
}

gainorloss<-ifelse(data$state=="1","LOH",
                   ifelse(data$state=="2","LOH",
                          ifelse(data$state=="3","neutral",
                                 ifelse(data$state=="4","Gain",
                                        ifelse(data$state=="5","Gain",
                                               ifelse(data$state=="6","Gain","")
                                        )))))
data$cytob<-arm
data$gainorloss<-gainorloss

library(stringr)
chr<-str_replace(data$chr,pattern = "chr",replacement = "")
infor<-paste(chr,data$cytob,sep = "")
infor<-paste(infor,data$gainorloss,sep="_")
infor
data$infor<-infor

tmp<-data$cytob=="other"|data$gainorloss=="neutral"
table(!tmp)

data<-data[!tmp,]


cell_group<-unique(data$cell_group_name)
length(cell_group)
data3<-data[data$cell_group_name==cell_group[1],]
data3<-data3[!duplicated(data3$infor),]
data2<-data3


for (i in 2:length(cell_group)){
  data3<-data[data$cell_group_name==cell_group[i],]
  data3<-data3[!duplicated(data3$infor),]
  data2<-rbind(data2,data3)
}

data3<-data.frame(data2[,c("cell_group_name","infor")])
data3$infor2<-1
data4<-reshape2::dcast(data3,cell_group_name~infor,value.var = "infor2")

#data4<-data4[-c(9,10),]

####read in the cell_grouping information###
BC21group<-read.table(file="BC21.cell_groupings.txt",header=T)
BC21group_2<-grepl(pattern="references",BC21group$cell_group_name)
BC21group<-BC21group[!BC21group_2,]
####combine the cells and the cytoband information###
merge<-dplyr::left_join(x = BC21group,y=data4,by=c("cell_group_name"))
merge[is.na(merge)]<-0
nrow(merge)
percent<-colSums(merge[,-c(1:2)])/nrow(merge)
dim(merge)
percent
write.csv(percent,file="BC21_percent.csv")
write.csv(data4,file="BC21_cell_canonical.csv")






####分析BC22####
data<-read.table(file="BC22.dat.txt",header=T)
arm<-vector()
for (i in 1:nrow(data)){
  #site<-data[i,c("chr","start","end")]
  #cyto<-cytoband[site$chr,]
  site<-data[i,c("chr","start","end")]
  cyto<-cytoband[chrom==site$chr,]
  
  arm[i]<-ifelse(site$start > cyto[2,]$start,"q",
                 ifelse(site$end < cyto[1,]$end,"p","other"))
  
}

gainorloss<-ifelse(data$state=="1","LOH",
                   ifelse(data$state=="2","LOH",
                          ifelse(data$state=="3","neutral",
                                 ifelse(data$state=="4","Gain",
                                        ifelse(data$state=="5","Gain",
                                               ifelse(data$state=="6","Gain","")
                                        )))))
data$cytob<-arm
data$gainorloss<-gainorloss

library(stringr)
chr<-str_replace(data$chr,pattern = "chr",replacement = "")
infor<-paste(chr,data$cytob,sep = "")
infor<-paste(infor,data$gainorloss,sep="_")
infor
data$infor<-infor

tmp<-data$cytob=="other"|data$gainorloss=="neutral"
table(!tmp)

data<-data[!tmp,]


cell_group<-unique(data$cell_group_name)
length(cell_group)
data3<-data[data$cell_group_name==cell_group[1],]
data3<-data3[!duplicated(data3$infor),]
data2<-data3


for (i in 2:length(cell_group)){
  data3<-data[data$cell_group_name==cell_group[i],]
  data3<-data3[!duplicated(data3$infor),]
  data2<-rbind(data2,data3)
}

data3<-data.frame(data2[,c("cell_group_name","infor")])
data3$infor2<-1
data4<-reshape2::dcast(data3,cell_group_name~infor,value.var = "infor2")

#data4<-data4[-c(9,10),]

####read in the cell_grouping information###
BC22group<-read.table(file="BC22.cell_groupings.txt",header=T)
BC22group_2<-grepl(pattern="references",BC22group$cell_group_name)
BC22group<-BC22group[!BC22group_2,]
####combine the cells and the cytoband information###
merge<-dplyr::left_join(x = BC22group,y=data4,by=c("cell_group_name"))
merge[is.na(merge)]<-0
nrow(merge)
percent<-colSums(merge[,-c(1:2)])/nrow(merge)
dim(merge)
percent
write.csv(percent,file="BC22_percent.csv")
write.csv(data4,file="BC22_cell_canonical.csv")


####分析BC16####
data<-read.table(file="BC16.dat.txt",header=T)
arm<-vector()
for (i in 1:nrow(data)){
  #site<-data[i,c("chr","start","end")]
  #cyto<-cytoband[site$chr,]
  site<-data[i,c("chr","start","end")]
  cyto<-cytoband[chrom==site$chr,]
  
  arm[i]<-ifelse(site$start > cyto[2,]$start,"q",
                 ifelse(site$end < cyto[1,]$end,"p","other"))
  
}

gainorloss<-ifelse(data$state=="1","LOH",
                   ifelse(data$state=="2","LOH",
                          ifelse(data$state=="3","neutral",
                                 ifelse(data$state=="4","Gain",
                                        ifelse(data$state=="5","Gain",
                                               ifelse(data$state=="6","Gain","")
                                        )))))
data$cytob<-arm
data$gainorloss<-gainorloss

library(stringr)
chr<-str_replace(data$chr,pattern = "chr",replacement = "")
infor<-paste(chr,data$cytob,sep = "")
infor<-paste(infor,data$gainorloss,sep="_")
infor
data$infor<-infor

tmp<-data$cytob=="other"|data$gainorloss=="neutral"
table(!tmp)

data<-data[!tmp,]


cell_group<-unique(data$cell_group_name)
length(cell_group)
data3<-data[data$cell_group_name==cell_group[1],]
data3<-data3[!duplicated(data3$infor),]
data2<-data3


for (i in 2:length(cell_group)){
  data3<-data[data$cell_group_name==cell_group[i],]
  data3<-data3[!duplicated(data3$infor),]
  data2<-rbind(data2,data3)
}

data3<-data.frame(data2[,c("cell_group_name","infor")])
data3$infor2<-1
data4<-reshape2::dcast(data3,cell_group_name~infor,value.var = "infor2")

#data4<-data4[-c(9,10),]

####read in the cell_grouping information###
BC16group<-read.table(file="BC16.cell_groupings.txt",header=T)
BC16group_2<-grepl(pattern="references",BC16group$cell_group_name)
BC16group<-BC16group[!BC16group_2,]
####combine the cells and the cytoband information###
merge<-dplyr::left_join(x = BC16group,y=data4,by=c("cell_group_name"))
merge[is.na(merge)]<-0
nrow(merge)
percent<-colSums(merge[,-c(1:2)])/nrow(merge)
dim(merge)
percent
write.csv(percent,file="BC16_percent.csv")
write.csv(data4,file="BC16_cell_canonical.csv")

####draw the heatmap of the frequency of the chromosome deletion or amplification####
list.files()

bc2<-read.csv("./BC2/BC2_percent.csv",header = T)
colnames(bc2)<-c("site","BC2_frequency")

bc3<-read.csv("./BC3/BC3_percent.csv",header = T)
colnames(bc3)<-c("site","BC3_frequency")

bc5<-read.csv("./BC5/BC5_percent.csv",header = T)
colnames(bc5)<-c("site","BC5_frequency")

bc6<-read.csv("./BC6/BC6_percent.csv",header = T)
colnames(bc6)<-c("site","BC6_frequency")

bc10<-read.csv("./BC10/BC10_percent.csv",header = T)
colnames(bc10)<-c("site","BC10_frequency")

bc11<-read.csv("./BC11/BC11_percent.csv",header = T)
colnames(bc11)<-c("site","BC11_frequency")

bc17<-read.csv("./BC17/BC17_percent.csv",header = T)
colnames(bc17)<-c("site","BC17_frequency")

bc16<-read.csv("./BC16/BC16_percent.csv",header = T)
colnames(bc16)<-c("site","BC16_frequency")

bc20<-read.csv("./BC20/BC20_percent.csv",header = T)
colnames(bc20)<-c("site","BC20_frequency")

bc21<-read.csv("./BC21/BC21_percent.csv",header = T)
colnames(bc21)<-c("site","BC21_frequency")

bc22<-read.csv("./BC22/BC22_percent.csv",header = T)
colnames(bc22)<-c("site","BC22_frequency")

library(dplyr)

join1<-full_join(bc2,bc3,by="site")
join2<-full_join(bc5,bc6,by="site")
join3<-full_join(bc10,bc11,by="site")
join4<-full_join(bc16,bc17,by="site")
join5<-full_join(bc20,bc21,by="site")

join6<-full_join(join5,bc22,by="site")
join7<-full_join(join1,join2,by="site")
join8<-full_join(join3,join4,by="site")
join9<-full_join(join7,join8,by="site")
join10<-full_join(join9,join6,by="site")

join_total<-c("1p_Gain","1p_LOH","1q_Gain","1q_LOH","2p_Gain","2p_LOH","2q_Gain","2q_LOH",
              "3p_Gain","3p_LOH","3q_Gain","3q_LOH","4p_Gain","4p_LOH","4q_Gain","4q_LOH",
              "5p_Gain","5p_LOH","5q_Gain","5q_LOH","6p_Gain","6p_LOH","6q_Gain","6q_LOH",
              "7p_Gain","7p_LOH","7q_Gain","7q_LOH","8p_Gain","8p_LOH","8q_Gain","8q_LOH",
              "9p_Gain","9p_LOH","9q_Gain","9q_LOH","10p_Gain","10p_LOH","10q_Gain","10q_LOH",
              "11p_Gain","11p_LOH","11q_Gain","11q_LOH","12p_Gain","12p_LOH","12q_Gain","12q_LOH",
              "13p_Gain","13p_LOH","13q_Gain","13q_LOH","14p_Gain","14p_LOH","14q_Gain","14q_LOH",
              "15p_Gain","15p_LOH","15q_Gain","15q_LOH","16p_Gain","16p_LOH","16q_Gain","16q_LOH",
              "17p_Gain","17p_LOH","17q_Gain","17q_LOH","18p_Gain","18p_LOH","18q_Gain","18q_LOH",
              "19p_Gain","19p_LOH","19q_Gain","19q_LOH","20p_Gain","20p_LOH","20q_Gain","20q_LOH",
              "21p_Gain","21p_LOH","21q_Gain","21q_LOH","22p_Gain","22p_LOH","22q_Gain","22q_LOH")
joint_total<-as.data.frame(join_total)
#rownames(join_total)<-joint_total$join_total
?left_join
jo<-left_join(joint_total,join10,by=c("join_total"="site"))

jo[is.na(jo)]<-0
rownames(jo)<-jo$join_total

jo<-jo[,-1]

library(pheatmap)
colnames(jo)
ac<-data.frame(sampleID=colnames(jo))
rownames(ac)<-colnames(jo)
ac$sampleID<-factor(ac$sampleID,levels=colnames(jo))

?pheatmap
library(RColorBrewer)
pheatmap(jo,annotation_col = ac,cluster_rows = F,cluster_cols = F,annotation_names_col = T,show_colnames = F,
        fontsize_row = 6,color = colorRampPalette(brewer.pal(n = 7, name ="PuRd"))(100),fontsize = 8)


write.csv(jo,file="samples_cytoband_deficiency.csv")
