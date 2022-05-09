# genetic distance for non-S proteins
gd.nonS <- <- function(fn,vaccine,outpath,prot) {
  df <- read.csv(fn,colClasses = c("factor"))
  cb_df <- rbind(vaccine,df)
  season <- c(1,rep(2,length(df[,1])))
  cb_df1 <- cbind(cb_df[,1:3],season,cb_df[,4:length(df)])
  pv <- c()
  for (j in 5:(length(cb_df1[1,]))) {
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
  for (j in 5:(length(cb_df1))) {
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
  
  if (length(sel_pv[,1]) != 0) {
    rn <- strsplit(rownames(sel_pv),"\\.")  
    site <- c()
    aa <- c()
    for (n in 1:length(sel_pv[,1])) {
      s <- rn[[n]][1]
      a<- rn[[n]][2]
      site <- c(site,s)
      aa <- c(aa,a)
    }
    sel_pv1 <- cbind(site,aa,sel_pv[,2:3])
    row.names(sel_pv1) <- 1:length(sel_pv[,1])
    out_fileName <- paste(str_name[[1]][1],str_name[[1]][2],str_name[[1]][3],
                          str_name[[1]][4],prot,"mismatch sites.csv", sep =' ')
    write.csv(sel_pv1,file = paste(outpath,out_fileName,sep = "/"),row.names = FALSE)
    table_site <- sel_pv1
    table_site$site <- as.numeric(table_site$site)
    table_site$prop <- as.numeric(table_site$prop)
    mismatch <- colSums(table_site[,4,drop=FALSE])
    ss <- length(df[,1])
    cb_mismatch <- c(mismatch,ss)
    names(cb_mismatch) <- c("All","sample size")
    return(cb_mismatch) 
  } else {
    ss <- length(df[,1])
    cb.mismatch <- c(0,ss)
    return(cb_mismatch)
  }
} 

# vaccine sequence
# set reference sequence path
pr1 <- pr[-18]
prot1 <- paste("SARS-CoV-2",pr1,"protein",sep=" ")
for (p in 1:length(pr1)) {
  ref.strain <- paste(ref_path,"/Wuhan_reference_",pr1[p],".csv",sep="")
  pre_seq <- read.csv(ref.strain,colClasses = c("factor"))
  # testing sequence
  # set input path
  setwd(paste(inputpath,prot1[p],sep = "/"))
  files <- list.files(pattern = ".csv")
  # set output path to export mismatch sites
  outpath <- paste(path,prot1[p],sep="/")
  dir.create(outpath)
  mm_table <- c()
  for (f in 1:length(files)) {
    mismatch <- gd.nonS(files[f],pre_seq,outpath,pr1[p])
    mm_table <- rbind(mm_table,mismatch)
  }
  row.names(mm_table) <- files
  colnames(mm_table) <- c(pr1[p],"sample size")
  # set path to save non-S protein mismatch summary tables
  write.csv(mm_table,file=paste(tarPath,paste(pr1[p],"mismatch.csv",sep=" "),sep="/"),row.names = TRUE)
}

setwd(tarPath)
mm_table <- read.csv(file = paste("pr1[1],"mismatch.csv",sep=" "),colClasses = c("factor"))
for (p in 2:length(pr1)) {
  mismatch <- read.csv(file = paste("pr1[p],"mismatch.csv",sep=" "),colClasses = c("factor"))
  mm_table <- cbind(mm_table,mismatch[,2])
}
colnames(mm_table)[1] <- c("filenames")
colnames(mm_table)[4:31] <- pr1[2:29]
names <- gsub("nsp1", "S", mm_table$filenames)
mm_table$filenames <- names
mm_table1 <- mm_table[,-3]
write.csv(mm_table1,file="non-S protein mismatch.csv",sep=" ",row.names = TRUE)

df <- read.csv("VE and S mismatch.csv", header = TRUE, fill = TRUE)
ac_mm <- c()
for (i in 1:length(df[,1])) {
  row <- match(df$filenames[i],mm_table1$filenames)
  extract_mismatch <- mm_table1[row,]
  ac_mm <- rbind(ac_mm,extract_mismatch)
}
colnames(ac_mm) <- colnames(mm_table1)
df1 <- cbind(df,ac_mm[,2:30])
write.csv(df1,"VE and whole genome mismatch.csv",row.names = FALSE)
