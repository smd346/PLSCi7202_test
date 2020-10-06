##########Uplift to V4 of genome and Manhattan Plot Creation for TWAS ionomics results############
#read in TWAS results
library(qqman)
setwd("C:/Users/Vanah/Documents/Cornell/Gorelab/TWAS_ionomics/TWAS_results/12.16.19/")
Tissues<-c("GRoot","GShoot","Kern","L3Base","L3Tip","LMAD","LMAN")
Traits<- c("B", "Ca", "Cu", "Fe", "K", "Mg", "Mn", "Mo", "Ni", "P", "Zn")

for(tissue in Tissues){
  #tissue <- Tissues[5]
  for(trait in Traits){
    results <- read.table(paste("C:/Users/Vanah/Documents/Cornell/Gorelab/TWAS_ionomics/TWAS_results/12.16.19/", "no_pos_",tissue,"_",trait,"_P5K25_TWAS_PEERgene_results.txt", sep =""), header = T)
    colnames(results)[1] <- "v3_gene_model"
    
    
    ##calc FDR, rank and percentile for each line
    results$FDR_adjusted <- p.adjust(results$pvalue, method = "fdr")
    
    results$rank <- rank(results$FDR_adjusted)
    
    results$percentile <- results$rank/nrow(results)*100
    
    #merge with file containing position
    results_with_pos <- read.table(paste("C:/Users/Vanah/Documents/Cornell/Gorelab/TWAS_ionomics/TWAS_results/12.16.19/",tissue,"_",trait,"_P5K25_TWAS_PEERgene_results.txt", sep =""), header = T)
    colnames(results_with_pos)[1]<- "v3_gene_model"
    results_with_pos <- merge(results, results_with_pos, by = "v3_gene_model")
    results_with_pos <- subset(results_with_pos[1:12])
    colnames(results_with_pos)[2] <- "pvalue"
    colnames(results_with_pos)[4] <- "r"
    #merging with position file results in losing 201 genes
    
    
    #uplift to v4 reference genome positions
    #need to merge the results_with_pos file so that we have all of the information in one file
    key_file <- read.delim("C:/Users/Vanah/Documents/Cornell/Gorelab/TWAS_ionomics/gene_model_xref_v4.txt", header = F)
    key_file <- key_file[c(1:8,11)] #get rid of unecessary columns
    colnames(key_file)[1] <- "v4_gene_model"
    colnames(key_file)[2] <- "v4_start_pos"
    colnames(key_file)[3] <- "v4_end_pos"
    colnames(key_file)[4] <- "v4_chr"
    colnames(key_file)[5] <- "v4_gm_genbank"
    colnames(key_file)[6] <- "v4_transcript"
    colnames(key_file)[7] <- "v4_transcript_acc"
    colnames(key_file)[8] <- "v4_locus"
    colnames(key_file)[9] <- "v3_gene_model"
    results_uplifted <- merge(results_with_pos, key_file, by = "v3_gene_model", all = F) #this will produce a table with both v3 and v4 gene models, v3 start position and all v4 positions
    colnames(results_uplifted)[9] <- "v3_chr"
    #colnames(results_uplifted)[13]<- "v4_gene_model"
    #colnames(results_uplifted)[12] <- "Start_pos_v3"
    #lose 4063 genes when you uplift
    
    
    ##write table with the adjusted p-values, FDR, rank, percentile and gene position
    write.table(results_uplifted, paste("C:/Users/Vanah/Documents/Cornell/Gorelab/TWAS_ionomics/TWAS_results/12.16.19/",tissue,"_",trait,"_P5K25_TWAS_PEERgene_results_FDR_with_pos_uplifted.txt", sep=""), col.names= T,row.names=FALSE,quote=F,sep="\t")
    
    #calc p-value threshold (FDR line to go on manhattan plot), set defaults for threshold if nothing passes
    #these are calculated using the file without position in order to ensure all genes are included in calculation
    FDR_data <- results
    FDR_lessThan_0.10 <- subset(FDR_data, FDR_data$FDR_adjusted < 0.10)
    FDR_greaterThan_0.10 <- subset(FDR_data, FDR_data$FDR_adjusted > 0.10)
    p_max <- max(FDR_lessThan_0.10$pvalue)
    p_min <- min(FDR_greaterThan_0.10$pvalue)
    FDR_cutoff <- (p_max+p_min)/2
    if(is.integer(FDR_cutoff)){
      FDR_line <- log(FDR_cutoff, 10)
      print(paste(tissue,"_",trait,"_",FDR_cutoff," is ", FDR_cutoff, sep =""))
    }else{FDR_line <- 5
    print(paste(tissue,"_",trait,"_FDR_cutoff"," is ", FDR_line, sep =""))}
    
    
    #make manhattan plot
    #make sure to plot using the position from v4 reference genome
    plot_file <- read.delim(paste("C:/Users/Vanah/Documents/Cornell/Gorelab/TWAS_ionomics/TWAS_results/12.16.19/",tissue,"_",trait,"_P5K25_TWAS_PEERgene_results_FDR_with_pos_uplifted.txt", sep=""), header = T) #reads in the file with FDR, position, and uplifted positions
    plot_file <- plot_file[grep("B73V4_ctg*", plot_file$v4_chr, invert = TRUE),] #remove the rows containing genes that have not been mapped to a chromosome; lose ____ genes #check that these genes are not significant
    plot_file <- plot_file[!(is.na(plot_file$v4_chr) | plot_file$v4_chr==""), ] #remove rows that don't contain a value for chromosome
    plot_file$v4_chr <- gsub(patter = "Chr*", replacement = "", x= plot_file$v4_chr)
    plot_file$v4_chr <- as.numeric(plot_file$v4_chr)
    pdf(paste(tissue,"_",trait,"_TWAS_manhattan_plot.pdf", sep ="")) #opens a file to save plot in
    manhattan(plot_file, chr = "v4_chr", bp = "v4_start_pos", p = "pvalue", col = c("blue4", "orange3"), suggestiveline = FDR_line)
    dev.off() #saves plot as pdf
    
    
    
    #make qqplot
    pdf(paste(tissue,"_",trait,"_TWAS_qqplot.pdf",sep ="")) #opens a file to save plot in
    qq(plot_file$pvalue)
    dev.off() #saves plot as pdf
    #save qqplot as file
  }
}
