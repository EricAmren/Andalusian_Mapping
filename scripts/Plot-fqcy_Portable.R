#Plot-Fqcy.R generates plots of CEW1snp-frequency and outputs 2 tables:
# -1 summarizing data by scaffold
# -1 list of scaffolds displaying mean(JU170snp-frequency)<0.1

args <- commandArgs(trailingOnly = TRUE)

#Required Library
library("ggplot2")

#Define a multiplot function to display frequency-plot and coverage plot on the same window
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
   require(grid)
    # Make a list from the ... arguments and plotlist
   plots <- c(list(...), plotlist)
    numPlots = length(plots)
    # If layout is NULL, then use 'cols' to determine layout
   if (is.null(layout)) {
     # Make the panel
     # ncol: Number of columns of plots
     # nrow: Number of rows needed, calculated from # of cols
     layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
   }
  
  if (numPlots==1) {
     print(plots[[1]])
  
   } else {
     # Set up the page
     grid.newpage()
     pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
  
     # Make each plot, in the correct location
     for (i in 1:numPlots) {
       # Get the i,j matrix positions of the regions that contain this subplot
       matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
  
       print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                       layout.pos.col = matchidx$col))
     }
   }
 }

#Read arguments
srcFile <- args[1]
x<-read.table(srcFile, dec="." , header=F)

#List the scaffolds present in the file
scaff=levels(x$V1)

#Create an empty table to store the summary of data
tabtest <- as.data.frame(setNames( c(replicate(3,numeric(0), simplify = F)), c("Position", "Coverage", "JU170_Frequency")))
colonne1=as.data.frame(setNames(replicate(1,character(0), simplify = F), "Contigs"))
tabfinal=data.frame(colonne1, tabtest)
levels(tabfinal$Contigs)=scaff

#Loop for all scaffolds of the assembly (at least those having enough snps so that they had a snp-frequency computed)
for (i in 1:length(scaff)) {
  outFile <- paste0(args[2], "_", scaff[i], ".pdf")
  pdf(outFile)
  
  #Generate Freq plots
  p1=qplot(x[x$V1==scaff[i],2],x[x$V1==scaff[i],8], xlim=c(0,4600000), ylim=c(0,1), xlab="SNP position", ylab="frequence of mapping-strain snp") + 
    stat_smooth(method="loess", se=FALSE) +
    ggtitle("Allele Frequency")
  #meanfreq=mean(x[x$V1==scaff[i],8])
  #p1<-a + geom_hline(yintercept=meanfreq, colour="red", size=1.5, linetype="dashed") 
  
  #Generate Coverage plots
  p2<-qplot(x[x$V1==scaff[i],2],x[x$V1==scaff[i],6], xlim=c(0,4600000), xlab="SNP position", ylab="Coverage") + stat_smooth(method="loess", se=FALSE) + ggtitle("coverage")
  
  multiplot(p1, p2, cols=1)
  dev.off()
  
  #Generate summary by scaffold of mean freq and mean coverage
  tabfinal[i,]=c(scaff[i],length(x[x$V1==scaff[i],1]), mean(x[x$V1==scaff[i],6]), mean(x[x$V1==scaff[i],8]))
}

# keep goodscaffs (CEW1snp-frequency<0,3)
tab_selection <- tabfinal[tabfinal$JU170_Frequency<=0.1,c(1,4)]

#Output tables
write.table(tabfinal, file=paste0(args[2], "_summary.tsv"), sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE, append=TRUE)
write.table(tab_selection, file=paste0(args[2], "_CandidateScaff.tsv"), sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE, append=TRUE)
