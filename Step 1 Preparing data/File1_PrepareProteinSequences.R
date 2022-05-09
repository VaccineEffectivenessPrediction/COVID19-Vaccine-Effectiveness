# read fasta files

library("Biostrings")
read.fas <- function(fastaName,st,ed,protein,path) {
  seq1 <- readDNAStringSet(fastaName)
  seq2 <- DNAStringSet(seq1,start=st, end=ed)
  seq3 <- DNAStringSet(gsub("-","N",seq2))
  se <- translate(seq3, if.fuzzy.codon="X") 
  y <- strsplit(se@ranges@NAMES,"\\|")
  accession_number <- NULL
  date <- NULL
  strain_names <- NULL
  seq <- NULL
  for (i in 1:length(y)){
    accession_number <- c(accession_number,y[[i]][(length(y[[2]])-1)])
    date <- c(date,y[[i]][(length(y[[2]]))])
    strain_names <- c(strain_names,y[[i]][(length(y[[2]])-2)])
    seq <- c(seq,as.character(se[[i]]))
  }
  x <- cbind(accession_number,date,strain_names,seq)
  t <- strsplit(seq,"")
  sep_seq <- matrix(data="",ncol=length(t[[1]]),nrow=length(t))
  for (i in 1:length(t)){
    for (j in 1:length(t[[1]])) {
      sep_seq[i,j] <- t[[i]][j]
    }}
  sep_seq[sep_seq=="X"] <- NA
  sep_seq[sep_seq=="-"] <- NA
  colnames <- paste("V",1:length(t[[1]]),sep="")
  dimnames(sep_seq)=list(NULL,colnames)
  sel_seq <- cbind(accession_number,strain_names,date,sep_seq)
  sel_seq <- as.data.frame(sel_seq)
  sel_seq <- sel_seq[-1,] # Wuhan reference sequence used in multiple sequence alignment
  index_seq <- duplicated(sel_seq$strain_names)
  nodup_seq <- sel.seq[!index_seq,]
  seq_order <- nodup_seq[order(nodup_seq$date),]
  str_name <- strsplit(fastaName,".fas")
  out_fileName <- paste(str_name[[1]][1],"_",protein, ".csv", sep='')
  write.csv(seq_order,file=paste(path,out_fileName,sep = "/"),row.names = FALSE)
}

# set input path

# input fasta files
fasfiles <- dir(input_path)
setwd(input_path)
n <- 1:16
orf <- paste("nsp",n,sep = "")
pr <- c(orf,"nsp12-1","S","ORF3a","E","M","ORF6","ORF7a","ORF7b","ORF8","N","ORF10")
prot <- paste("SARS-CoV-2",pr,"protein",sep=" ")

s.all <- c(266,806,2720,8555,10055,10973,11843,12092,12686,13025,13442,13468,16237,18040,19621,20659,13442,21563,25393,26245,26523,27202,27394,27756,27894,28274,29558)
e.all <- c(805,2719,8554,10054,10972,11842,12091,12685,13024,13441,13480,16236,18039,19620,20658,21552,13468,25384,26220,26472,27191,27387,27759,27887,28259,29533,29674)

# set output path

# save data as csv files
for (p in 1:length(prot)) {
  outpath <- paste(path,prot[p],sep = "/")
  dir.create(outpath)
  s <- s.all[p]
  e <- e.all[p]
  protein <- c(pr[p])
  for (n in 1:length(fasfiles) {  
    read.fas(fasfiles[18],s,e,protein,outpath)
  }
}
