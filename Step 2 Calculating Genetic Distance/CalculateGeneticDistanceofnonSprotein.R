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
