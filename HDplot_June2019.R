rm(list=ls())

library(vcfR)
library(ggplot2)
library(dplyr)
library(stringr)

setwd(dir="../../Desktop/PIRE REU/Genomics and Bioinformatics Workshop/Downloaded Files/")
#exampleLoci<-read.vcfR("example.vcf")
exampleLoci<-read.vcfR("BT-G80low40.recode.vcf")
#exampleLoci<-read.vcfR("GLC-g70low50mac3.recode.vcf")

exampleLoci@fix[,1]

#HDplot function
HDPlot<-function(vcfData){
  #set up results table
  HDplotTable<-as.data.frame(matrix(NA,nrow=dim(vcfData@gt)[1],ncol=10))
  colnames(HDplotTable)<-c("Chrom","Pos","depth_a","depth_b","ratio","num_hets","num_samples","het_perc","std","z")
  
  #get genotypes from vcf file
  genos<-extract.gt(vcfData, element = "GT", mask = FALSE, as.numeric = FALSE, return.alleles = FALSE, 
                    IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE)
  
  #get allele reads from vcf file
  reads<-extract.gt(vcfData, element = "AD", mask = FALSE, as.numeric = FALSE, return.alleles = FALSE, 
                    IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE)
    
  #replace . with 0
  reads<-gsub("\\.","0,0",reads)
  
  alleleReads<-apply(reads,2,function(x) str_split_fixed(x,",",2))
  alleleReads_1<-alleleReads[1:dim(reads)[1],]
  alleleReads_2<-alleleReads[dim(reads)[1]+1:dim(reads)[1],]
  #convert to numeric format
  alleleReads_1<-apply(alleleReads_1,2, function(x) as.numeric(x))
  alleleReads_2<-apply(alleleReads_2,2, function(x) as.numeric(x))
  #subset to heterozygous genotypes
  #make genotype matrix where heterozygotes are 1 and other genotypes are 0
  hetMatrix<-apply(genos,2,function(x) dplyr::recode(x,'0/0'=0,'1/1'=0,'./.'=0,'0/1'=1,'1/0'=1))
  #multiply read count matrices by heterozygote matrix to get read counts for heterozygous genotypes
  alleleReads_1_het<-alleleReads_1*hetMatrix
  alleleReads_2_het<-alleleReads_2*hetMatrix
  #rows are loci and columns are samples
  #sum reads per allele per locus for heterozygous samples
  A_reads<-apply(alleleReads_1_het,1,sum)
  B_reads<-apply(alleleReads_2_het,1,sum)
  totalReads<-A_reads+B_reads
  ratio<-A_reads/totalReads
  std<-sqrt(totalReads*0.5*0.5)
  z<- -(totalReads/2-A_reads)/std
  #get percent heterozygosity for each locus
  numHets<-apply(hetMatrix,1,sum)
  hetPerc<-numHets/dim(hetMatrix)[2]
  
  #assign results to HDplotTable
  HDplotTable$Chrom<-vcfData@fix[,1]
  HDplotTable$Pos<-vcfData@fix[,2]
  HDplotTable$Locus_ID
  HDplotTable$depth_a<-A_reads
  HDplotTable$depth_b<-B_reads
  HDplotTable$ratio<-ratio
  HDplotTable$num_hets<-numHets
  HDplotTable$num_samples<-dim(hetMatrix)[2]
  HDplotTable$het_perc<-hetPerc
  HDplotTable$std<-std
  HDplotTable$z<-z
  return(HDplotTable)
}


example_output<-HDPlot(exampleLoci)
write.csv(example_output, file="HDplot_results.csv")

high_het=example_output[which(example_output$het_perc>0.55),]
no_high_het=example_output[which(example_output$het_perc<0.55),]
no_hi_z=no_high_het[which(no_high_het$z>5),]
no_low_z=no_high_het[which(no_high_het$z< -5),]
blacklist=rbind(high_het,no_hi_z,no_low_z)
write.csv(blacklist, file="HDplot_results_filtered.csv")
write.table(blacklist[,1:2],file="hdplot_blacklist.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
#plot results
#deviation plot
ggplot()+geom_point(data=example_output,aes(x=het_perc,y=z),alpha=0.5)+xlim(0,1)
#ratio plot
ggplot()+geom_point(data=example_output,aes(x=het_perc,y=ratio),alpha=0.5)+xlim(0,1)
#ratio plot