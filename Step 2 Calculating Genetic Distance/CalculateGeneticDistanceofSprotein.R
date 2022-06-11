# genetic distance for S protein
gd.s <- function(fn,vaccine,output_path,segment) {
  df <- read.csv(fn,colClasses = c("factor"))
  cb_df <- rbind(vaccine,df)
  season <- c(1,rep(2,length(df[,1])))
  cb_df1 <- cbind(cb_df[,1:3],season,cb_df[,4:length(df)])
  pv <- c()
  for (j in 5:length(cb_df1)) {
    proptable <- function(j){
      prop.table(xtabs(~season+cb_df1[,j],data=cb_df1,exclude = c("NA")),1)
    }
    prop_table <- proptable(j)
    ns <- colnames(prop_table)
    cnames=paste((j-4):(j-4),ns,sep=".")
    colnames(prop_table)=cnames
    if (length(prop_table[,1])!=1){
      pv = cbind(pv,prop_table)
    }
  }
  frq <- c()
  for (j in 5:length(cb_df1)) {
    frqtable <- function(j){
      xtabs(~season+cb_df1[,j],data=cb_df1,exclude = c("NA"))
    }
    frq_table <- frqtable(j)
    ns <- colnames(frq_table)
    cnames=paste((j-4):(j-4),ns,sep=".")
    colnames(frq_table)=cnames
    if (length(frq_table[,1])!=1){
      frq = cbind(frq,frq_table)
    }
  }
  table <- rbind(frq,pv[2,])
  table1 <- as.data.frame(t(table))
  colnames(table1) <- c("original","frq","prop")
  str_name <- strsplit(fn," ")
  
  sel_pv <- matrix(0,nrow=length(table1[,1]),ncol=length(table1[1,]))
  sel_pv <- table1
  for (i in 1:length(table1[,1])) {
    if(table1[i,1]==0) {
      sel_pv[i,] <- table1[i,]
    } else {
      sel_pv[i,] <- c(0,0,0)
    }
  }
  sel_pv = sel_pv[-(which(rowSums(sel_pv)==0)),]
  rn <- strsplit(rownames(sel_pv),"\\.")  
  site <- c()
  aa <- c()
  for (n in 1:length(sel_pv[,1])) {
    s <- rn[[n]][1]
    a <- rn[[n]][2]
    site <- c(site,s)
    aa <- c(aa,a)
  }
  sel_pv1 <- cbind(site,aa,sel_pv[,2:3])
  row.names(sel_pv1) <- 1:length(sel_pv[,1])
  out_fileName <- paste(str_name[[1]][1],str_name[[1]][2],str_name[[1]][3],
                        str_name[[1]][4],"mismatch sites.csv", sep =' ')
  write.csv(sel_pv1,file = paste(output_path,out_fileName,sep = "/"),row.names = FALSE)
  
  table_site <- sel_pv1
  table_site$site <- as.numeric(table_site$site)
  table_site$prop <- as.numeric(table_site$prop)
  cb_mismatch <- c()
  for (t in 1:length(type)) {
    data <- subset(table_site,site %in% segment[[t]])
    mismatch <- colSums(data[,4,drop=FALSE])
    cb_mismatch <- c(cb_mismatch,mismatch)
  }
  ss <- length(df[,1])
  cb_mismatch <- c(cb_mismatch,ss)
  names(cb_mismatch) <- c("NTD","RBD","All","sample size")
  return(cb_mismatch) 
}


# vaccine sequence
# set reference sequence path
ref_strain <- paste(ref_path,"/Wuhan_reference_",pr[18],".csv",sep="")
pre_seq <- read.csv(ref_strain,colClasses = c("factor"))

# testing sequence
# set input path
setwd(input_path)
files <- list.files(pattern = ".csv")
# set output path

# NTD
ntd <- 14:305
# RBD
rbd <- 319:541
# all
all <- 1:1273
segment <- list(ntd,rbd,all)

mm_table <- c()
for (f in 1:length(files)) {
  mismatch <- gd.s(files[f],pre_seq,outpath,segment)
  mm_table <- rbind(mm_table,mismatch)
}
row.names(mm_table) <- files

# collect VE data from literature
tp <- read.csv(file = "VE data.csv",colClasses = c("factor"))

ac_mm <- c()
for (i in 1:length(tp[,1])) {
  row <- match(tp$filenames[i],row.names(mm_table))
  extract_mismatch <- mm_table[row,]
  ac_mm <- rbind(ac_mm,extract_mismatch)
}
colnames(ac_mm) <- colnames(mm_table)
tp1 <- cbind(tp,ac_mm)
# set output path
write.csv(mm_table,paste(output_path,"S protein mismatch.csv",sep="/"))
write.csv(tp1,paste(output_path,"VE and S mismatch.csv",sep="/"),row.names = FALSE)

