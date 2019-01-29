#### Original author: Marc A. Beal
#### Modified: Matthew Meier
####In order to proceed the packages need to be installed
#source("http://bioconductor.org/biocLite.R")
#biocLite("SomaticSignatures")
#biocLite("VariantAnnotation")
#biocLite("Rsamtools")
#biocLite("deconstructSigs")

#######################
###Libraries to load###
#######################
library(SomaticSignatures)
library(VariantAnnotation)
library(Rsamtools)
library(reshape)
library(ggplot2)

###Give directory to required files in same directory (lacZ-normalized signatures, reference file, commoncalls1.txt)

	#args <- commandArgs(TRUE)

	#requiredFilesDirectory <- args[1]

	requiredFilesDirectory <- "./"

	fileName <- "Signature_input.txt"

	sigInput <- read.table(paste(requiredFilesDirectory, fileName, sep=""), header=T)

#Load reference sequence
fa_A = FaFile(paste(requiredFilesDirectory, "lacZ.fa", sep=""))



#Prepare VRange and mutation motifs
vr <- list()
df_vr <- list()
####Make sure group is in file 
for (i in 1:length(levels(sigInput$Group))) {

tempVRange <- VRanges(seqnames=rep("lacZ", nrow(sigInput[sigInput$Group==levels(sigInput$Group)[i],])),
	ranges = IRanges(sigInput[sigInput$Group==levels(sigInput$Group)[i],]$Position, end = sigInput[sigInput$Group==levels(sigInput$Group)[i],]$Position),
	ref = sigInput[sigInput$Group==levels(sigInput$Group)[i],]$Ref,
	alt = sigInput[sigInput$Group==levels(sigInput$Group)[i],]$Alt,
	sampleNames = rep(levels(sigInput$Group)[i], nrow(sigInput[sigInput$Group==levels(sigInput$Group)[i],])))

	vr[[i]] <- mutationContext(tempVRange, fa_A)
	df_vr[[i]] <- data.frame(t(data.frame(motifMatrix(vr[[i]], group = "sampleNames", normalize = TRUE))))

}

library(deconstructSigs)

for (i in 1:length(levels(sigInput$Group))) {
	colnames(df_vr[[i]]) <- colnames(signatures.cosmic)
}

###original
#lacZ_cosmic_signatures <- as.data.frame(t(read.table(paste(requiredFilesDirectory, "lacZ corrected COSMIC frequencies.txt", sep=""), sep="\t", header=FALSE)))
###original


lacZ_cosmic_signatures <- as.data.frame(t(read.table(paste(requiredFilesDirectory, "lacZ corrected COSMIC frequencies_with_controls(col31).txt", sep=""), sep="\t", header=FALSE)))

###original
#rownames(lacZ_cosmic_signatures) <- rownames(signatures.cosmic)
###original

rownames(lacZ_cosmic_signatures) <- c(rownames(signatures.cosmic), "Control")
colnames(lacZ_cosmic_signatures) <- colnames(signatures.cosmic)


plotSigs <- list()

for (i in 1:length(levels(sigInput$Group))) {
	plotSigs[[i]] <- whichSignatures(tumor.ref=df_vr[[i]], signatures.ref=lacZ_cosmic_signatures, sample.id=rownames(df_vr[[i]]))
}

###plotSignatures(plotSigs[[i]])


#Extract the signature data for each vcf file

allSigs <- matrix(ncol = 31, nrow = length(levels(sigInput$Group)))

###original
#allSigs <- matrix(ncol = 30, nrow = length(levels(sigInput$Group)))
###original

vcfNames <- vector()

for (i in 1:length(levels(sigInput$Group))) {
	allSigs[i,] <- as.numeric(plotSigs[[i]]$weights)
	vcfNames <- c(vcfNames, rownames(plotSigs[[i]]$weights))
}

rownames(allSigs) <- vcfNames
colnames(allSigs) <- names(plotSigs[[1]]$weights)


resultsToReturn <- list()

resultsToReturn$signatureSummary <- allSigs
resultsToReturn$plotSigs <- plotSigs
resultsToReturn$lacZ_cosmic_signatures <- lacZ_cosmic_signatures





################################################################################
################################################################################
#############################Plotter Function Starts############################
################################################################################
################################################################################


#inData comes from the data produced from the previous function
#vcfNumbers are the corresponding vcf files that you want to visualize that can be plotted when lacZSigs is False
#If you would like to plot existing signatures, set lacZSigs to TRUE and give the signature numbers as a vector (example: c(3, 4, 18))


plotSignatureResults <- function(inData, vcfNumbers, lacZSigs=F, whichSigs) {


if (lacZSigs == T) {
	lacZ_cosmic_signatures <- inData$lacZ_cosmic_signatures
	lacZ_cosmic_signatures$Signature <- row.names(lacZ_cosmic_signatures)
	lacZ_cosmic_signatures_melted <- melt(lacZ_cosmic_signatures, id.vars = "Signature" )
	tempmuts <-gsub("[A-Z]\\[", "", lacZ_cosmic_signatures_melted$variable)
	lacZ_cosmic_signatures_melted$muttype <- gsub("\\][A-Z]", "", tempmuts)
	lacZ_cosmic_signatures_melted$Signature_factor <- factor(lacZ_cosmic_signatures_melted$Signature, c(unique(lacZ_cosmic_signatures_melted$Signature)))

	listSigs <- vector()
	for (i in 1:length(whichSigs)) {
		listSigs <- c(listSigs, paste("Signature.", whichSigs[i], sep=""))
	}

	ggplot(lacZ_cosmic_signatures_melted[lacZ_cosmic_signatures_melted$Signature_factor %in% listSigs,], aes(x=variable, y=value)) + geom_bar(stat="identity", aes(fill=muttype), width=0.8) + scale_fill_manual(values = c("#1ebff0", "#000000", "#e62725", "#cbcacb", "#a1cf64", "#edc8c5")) + facet_grid(Signature_factor ~ muttype, scales="free_x") + theme(panel.spacing.x=unit(0, "lines"), panel.grid.major = element_line(colour="grey", size=0.5), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(breaks=seq(0,1,0.1), minor_breaks = NULL) + scale_x_discrete(breaks=NULL) + coord_cartesian(ylim=c(0,0.4))

}

if (lacZSigs == F) {
	vcfsToPlot <- matrix(ncol=96)
	colnames(vcfsToPlot) <- colnames(inData$plotSigs[[1]]$tumor)
	
	for (i in 1:length(vcfNumbers)) {	
		vcfsToPlot <- rbind(vcfsToPlot, inData$plotSigs[[vcfNumbers[i]]]$tumor)
		}
	vcfsToPlot <- as.data.frame(vcfsToPlot[-1,])

	vcfsToPlot$tissue <- rownames(vcfsToPlot)
	vcfsToPlot_melted <- melt(vcfsToPlot, id.vars = "tissue" )
	tempmuts <-gsub("[A-Z]\\[", "", vcfsToPlot_melted$variable)
	vcfsToPlot_melted$muttype <- gsub("\\][A-Z]", "", tempmuts)
	vcfsToPlot_melted$Tissue_factor <- factor(vcfsToPlot_melted$tissue, c(unique(vcfsToPlot_melted$tissue)))

	

	listvcfs <- vector()
	for (i in 1:length(vcfNumbers)) {
		listvcfs <- c(listvcfs, rownames(vcfsToPlot)[i])		
	}

	

	if (length(vcfNumbers) < 2) { print("vcfNumbers needs to include more than 1 data set to work properly") } else {
	ggplot(vcfsToPlot_melted[vcfsToPlot_melted$Tissue_factor %in% listvcfs,], aes(x=variable, y=value)) + geom_bar(stat="identity", aes(fill=muttype), width=0.8) + scale_fill_manual(values = c("#1ebff0", "#000000", "#e62725", "#cbcacb", "#a1cf64", "#edc8c5")) + facet_grid(Tissue_factor ~ muttype, scales="free_x") + theme(panel.spacing.x=unit(0, "lines"), panel.grid.major = element_line(colour="grey", size=0.5), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(breaks=seq(0,1,0.1), minor_breaks = NULL) + scale_x_discrete(breaks=NULL) + coord_cartesian(ylim=c(0,0.4))	
	}
}



}

################################################################################
################################################################################
##############################Plotter Function Ends#############################
################################################################################
################################################################################

#Use function to output plots to file

pdf(file=paste(requiredFilesDirectory, "Mutational_Profiles.pdf", sep=""))

plotSignatureResults(resultsToReturn, 1:length(rownames(resultsToReturn$signatureSummary)), lacZSigs=F)

dev.off()


################################################################################
################################################################################
####################Reconstruct & Pearson Coefficient Starts####################
################################################################################
################################################################################


reconstructSig <- function(inData) {

reconstructed <- list()
p.value <- vector()
cor <- vector()

	for (i in 1:nrow(inData$signatureSummary)) {
	reconstructed[[i]] <- inData$lacZ_cosmic_signatures * inData$signatureSummary[i,]
	result <- cor.test(as.vector(colSums(reconstructed[[i]])), as.vector(inData$plotSigs[[i]]$tumor), "two.sided", "pearson")
	
	p.value <- c(p.value, result$p.value)
	cor <- c(cor, result$estimate)

	}

resultsSummary <- rbind(p.value, cor)
colnames(resultsSummary) <- rownames(inData$signatureSummary)

return(resultsSummary)

}


################################################################################
################################################################################
#####################Reconstruct & Pearson Coefficient Ends#####################
################################################################################
################################################################################






################################################################################
################################################################################
############Correlation Between Ind. Sigs. and Mutation Data. Starts############
################################################################################
################################################################################

corTestSig <- function(inData) {

group <- list()

	for (i in 1:nrow(inData$signatureSummary)) {

	group[[i]] <- list()
	group[[i]]$groupName <- row.names(inData$signatureSummary)[i]
	
	tested <- vector()
	p.value <- vector()
	cor <- vector()
		
	corSigs <- as.integer(which(inData$signatureSummary[i,]>0))

		for (j in 1:length(corSigs)) {
		

		result <- cor.test(as.numeric(inData$lacZ_cosmic_signatures[corSigs[j],]), as.numeric(inData$plotSigs[[i]]$tumor), "two.sided", "pearson")		
		
		tested <- c(tested, paste("Signature.", corSigs[j], sep=""))
		p.value <- c(p.value, result$p.value)
		cor <- c(cor, result$estimate)

		}

	group[[i]]$result <- rbind(tested, p.value, cor)

	}

return(group)

}


################################################################################
################################################################################
#############Correlation Between Ind. Sigs. and Mutation Data Ends##############
################################################################################
################################################################################



###Remove zero sigs, merge summary with pearson (from reconstructSigs), print out raw summary

rawSummary <- resultsToReturn$signatureSummary[, colSums(resultsToReturn$signatureSummary != 0) > 0]

Residual <- 1 - t(t(rowSums(rawSummary)))

rawSummary <- cbind(rawSummary, Residual)

colnames(rawSummary)[ncol(rawSummary)] = "Residual"

rawSummary <- cbind(rawSummary, t(t(reconstructSig(resultsToReturn)[2,])))

colnames(rawSummary)[ncol(rawSummary)] = "reconstructSigPearsonCoefficient"

rawSummary <- t(rawSummary)


#Write table of raw signature summary
write.table(rawSummary, file=paste(requiredFilesDirectory, "Signature_raw_summary.txt", sep=""), sep="\t", col.names=NA)



###Write all correlations

correlationResults <- corTestSig(resultsToReturn)

sink(paste(requiredFilesDirectory, "Signature_Pearson_correlations_between_sigs_and_mut_profiles.txt", sep=""))
print(correlationResults)
sink()


#Find largest residual and filter out signatures with contributions below
largestResidual <- max(Residual)

finalSummary <- resultsToReturn$signatureSummary[, colSums(resultsToReturn$signatureSummary != 0) > 0]

finalSummary <- ifelse(finalSummary<largestResidual,0,finalSummary)



#now remove anything with pearson correlation below 0.5

for(i in 1:nrow(finalSummary)) {

tempData <- t(correlationResults[[i]]$result)
nonImportantSigs <- as.vector(t(tempData[tempData[,3]<=0.5,]))
nonImportantSigs <- nonImportantSigs[grep("Signature", nonImportantSigs)]
nonImportantSigs <- gsub("Signature.31", "Control", nonImportantSigs)
finalSummary[i,][nonImportantSigs] <- 0

#nonImportantSigs <- as.vector(tempData[tempData[,3]<=0.5,][,1])
#nonImportantSigs <- tempData[grep("Signature", tempData[tempData[,3]<=0.5,])]
#ifelse(is.na(nonImportantSigs), finalSummary <- finalSummary, finalSummary[i,][nonImportantSigs] <- 0)



}


#write final results that pass the thresholds

finalSummary <- cbind(finalSummary, t(t(reconstructSig(resultsToReturn)[2,])))

colnames(finalSummary)[ncol(finalSummary)] = "reconstructSigPearsonCoefficient"

finalSummary <- t(finalSummary)

write.table(finalSummary, file=paste(requiredFilesDirectory, "Signature_final_summary_above_thresholds.txt", sep=""), sep="\t", col.names=NA)





###Summarize output
sink(file=paste(requiredFilesDirectory, "README - Summary of Signature Data.txt", sep=""))
cat(paste("Analysis was completed on ", Sys.time(),sep=""))
cat("\n")
cat("\n")
cat("\n")
cat("The raw results from signature deconstruction can be found in Signature_raw_summary.txt.\n")
cat("\n")
cat("\n")
cat("The final result can be found in Signature_final_summary_above_thresholds.txt.\n")
cat("\n")
cat("\n")
cat(paste("Signatures that had a contribution below the largest residual of ", largestResidual, " were removed and reported as 0.\n"))
cat("\n")
cat("Furthermore, if the Pearson coefficient between the signature and mutation data was below 0.5, that signature was also removed and reported as 0.\n")
cat("The raw P.C. information is reported in the file Signature_Pearson_correlations_between_sigs_and_mut_profiles.txt")
cat("\n")
cat("\n")
cat("Best of luck with your data interpretation.")
sink()



























################################################################################
################################################################################
##############Reconstruct, Min Removal & Pearson Coefficient Starts#############
################################################################################
################################################################################

removeSig <- function(inData) {

group <- list()

	for (i in 1:nrow(inData$signatureSummary)) {

	group[[i]] <- list()
	group[[i]]$groupName <- row.names(inData$signatureSummary)[i]
	
	removed <- vector()
	p.value <- vector()
	cor <- vector()
		
	corSigs <- as.integer(which(inData$signatureSummary[i,]>0))

		for (j in 1:length(corSigs)) {
		reconstructed <- inData$lacZ_cosmic_signatures[-corSigs[j],] * inData$signatureSummary[i,-corSigs[j]]
		result <- cor.test(as.vector(colSums(reconstructed)), as.vector(inData$plotSigs[[i]]$tumor), "two.sided", "pearson")		
		
		removed <- c(removed, paste("Signature.", corSigs[j], sep=""))
		p.value <- c(p.value, result$p.value)
		cor <- c(cor, result$estimate)

		}

	group[[i]]$result <- rbind(removed, p.value, cor)

	}

return(group)

}






################################################################################
################################################################################
###############Reconstruct, Min Removal & Pearson Coefficient Ends##############
################################################################################
################################################################################






################################################################################
################################################################################
##############################Function Starts Here##############################
################################################################################
################################################################################

reDeconstruct <- function(vcfDirectory, requiredFilesDirectory) {

#Load reference sequence
fa_A = FaFile(paste(requiredFilesDirectory, "lacZ.fa", sep=""))


vcfFiles <- list.files(vcfDirectory)
vr <- list()
df_vr <- list()

for (i in 1:length(vcfFiles)) {
	vr[[i]] <- mutationContext(readVcfAsVRanges(paste(vcfDirectory, vcfFiles[i], sep=""), paste("Genome", i, sep="")), fa_A)
	df_vr[[i]] <- data.frame(t(data.frame(motifMatrix(vr[[i]], group = "sampleNames", normalize = TRUE))))
	
}

library(deconstructSigs)

for (i in 1:length(vcfFiles)) {
	colnames(df_vr[[i]]) <- colnames(signatures.cosmic)
}

lacZ_cosmic_signatures <- as.data.frame(t(read.table(paste(requiredFilesDirectory, "lacZ corrected COSMIC frequencies.txt", sep=""), sep="\t", header=FALSE)))


rownames(lacZ_cosmic_signatures) <- rownames(signatures.cosmic)
colnames(lacZ_cosmic_signatures) <- colnames(signatures.cosmic)


plotSigs <- list()

for (i in 1:length(vcfFiles)) {
	plotSigs[[i]] <- whichSignatures(tumor.ref=df_vr[[i]], signatures.ref=lacZ_cosmic_signatures, sample.id=rownames(df_vr[[i]]))
}

#Extract the signature data for each vcf file

allSigs <- matrix(ncol = 30, nrow = length(vcfFiles))
vcfNames <- vector()

for (i in 1:length(vcfFiles)) {
	allSigs[i,] <- as.numeric(plotSigs[[i]]$weights)
	vcfNames <- c(vcfNames, rownames(plotSigs[[i]]$weights))
}

rownames(allSigs) <- vcfNames
colnames(allSigs) <- names(plotSigs[[1]]$weights)


resultsToReturn <- list()

resultsToReturn$signatureSummary <- allSigs
resultsToReturn$plotSigs <- plotSigs
resultsToReturn$lacZ_cosmic_signatures <- lacZ_cosmic_signatures


sigToRemove <- vector()

for (i in 1:nrow(resultsToReturn$signatureSummary)) {
	sigToRemove <- c(sigToRemove, as.integer(which(resultsToReturn$signatureSummary[i,] == min(resultsToReturn$signatureSummary[i,resultsToReturn$signatureSummary[i,]>0]))))

}

plotSigs2 <- list()

for (i in 1:length(vcfFiles)) {
	plotSigs2[[i]] <- whichSignatures(tumor.ref=df_vr[[i]], signatures.ref=lacZ_cosmic_signatures[-sigToRemove[i],], sample.id=rownames(df_vr[[i]]))
}


allSigs2 <- matrix(ncol = 30, nrow = length(vcfFiles))
vcfNames2 <- vector()

for (i in 1:length(vcfFiles)) {
	allSigs2[i,] <- as.numeric(plotSigs2[[i]]$weights)
	vcfNames2 <- c(vcfNames, rownames(plotSigs2[[i]]$weights))
}


rownames(allSigs2) <- vcfNames2
colnames(allSigs2) <- names(plotSigs2[[1]]$weights)


resultsToReturn2 <- list()

resultsToReturn2$signatureSummary <- allSigs
resultsToReturn2$plotSigs <- plotSigs
resultsToReturn2$lacZ_cosmic_signatures <- lacZ_cosmic_signatures


}




################################################################################
################################################################################
###############################Function Ends Here###############################
################################################################################
################################################################################


