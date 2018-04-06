args <- commandArgs(TRUE)

# TO RUN ON MANY FILES, USE THE FOLLOWING COMMAND IN BASH:
# for i in $(seq -w 0 297) ; do echo $i; cov=$(echo "CNV_${i}.1000_windows.median_TGR_coverage.bed"); acgh=$(echo "CNV_${i}.combined_aCGH_values.bed"); echo $cov; echo $acgh ; Rscript plotting_script.R $cov $acgh CNV_${i} CNV_${i}_original CNV_${i}_CNVnator CNV_${i}_manta; done

library('ggplot2')
#library('reshape2')
library('tidyr')
library('dplyr')

coverage <- read.csv(args[1], header=F,sep="\t")

# rename columns
#chr
#start
#end
#Depth_for_BigBlueMouse
#Depth_for_MutaMouse13
#Depth_for_MutaMouse14
#Depth_for_MutaMouse17
#Depth_for_MutaMouse18
#Depth_for_MutaMouse19
#Depth_for_MutaMouseEuroFemale
#Depth_for_MutaMouseEuroMale
#Depth_for_gptDeltaMouse
#Depth_for_lacZplasmidMouse
colnames(coverage) <- c("chr","start","end", "BigBlue Mouse", "MutaMouse 1", "MutaMouse 2", "MutaMouse 3", "MutaMouse 4", "MutaMouse 5", "MutaMouse Covance Female", "MutaMouse Covance Male", "gpt delta mouse", 'lacZ plasmid mouse')
tidy_coverage <- gather(coverage, mouse, value, -chr, -start, -end)
tidy_coverage$type <- "NGS Coverage"

acgh <- read.csv(args[2], header=F,sep="\t")
# rename columns
#chr
#start
#end
#MutaMouse 13
#MutaMouse 14
#MutaMouse 17
#MutaMouse 18
#MutaMouse 19
#MutaMouse Euro Female
#MutaMouse Euro Male
colnames(acgh) <- c("chr","start","end","MutaMouse 1", "MutaMouse 2", "MutaMouse 3", "MutaMouse 4", "MutaMouse 5", "MutaMouse Covance Female", "MutaMouse Covance Male")
tidy_array <- gather(acgh, mouse, value, -chr, -start, -end)
tidy_array$type <- "aCGH"

tidy_data <- bind_rows(tidy_array,tidy_coverage)

# names to plot
mouse_names <- c("MutaMouse 1",  "MutaMouse 2", "MutaMouse 3", "MutaMouse 4", "MutaMouse 5" )
mouse_names_covance <- c("MutaMouse 1",  "MutaMouse 2", "MutaMouse 3", "MutaMouse 4", "MutaMouse 5", "MutaMouse Covance Female", "MutaMouse Covance Male" )
cnvLine <- read.csv(args[4], header=F, sep="\t")
colnames(cnvLine) <- c("chr", "start", "end", "name")
cnvLine$type <- "aCGH"

mantaColor="red"
cnvNatorColor="black"

if(file.exists(args[5])) {
        cnvNator <- tryCatch(read.csv(args[5], header=F, sep="\t"), error=function(e) { cnvNator <- NULL })
        if (!is.null(cnvNator)) {
        print(cnvNator)
        colnames(cnvNator) <- c("chr", "start", "end")
        cnvNator$type <- "NGS Coverage"
        }
}

if(file.exists(args[5])) {
        manta <-  tryCatch(read.csv(args[6], header=F, sep="\t"), error=function(e) { manta <- NULL })
        if (!is.null(manta)) {
        print(manta)
        colnames(manta) <- c("chr", "start", "end")
        manta$type <- "NGS Coverage"
        }
}

if (is.null(cnvNator)) {
        cnvNatorColor="white"
        cnvNator <- cnvLine
        cnvNator$type <- "NGS Coverage"
}

if (is.null(manta)) {
        mantaColor="white"
        manta <- cnvLine
        manta$type <- "NGS Coverage"
}

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# ggplot(data_melted, aes(x=V2,y=value,color=variable)) + geom_point()

ggplot(tidy_data[tidy_data$mouse %in% mouse_names,], aes(x=start,y=value,color=mouse)) +
        geom_point(size=1, alpha=0.4) +
        facet_grid(type ~ . , scales = "free") +
        ### geom_segment(data=cnvLine, aes(x=cnvLine$start, xend=cnvLine$end, y=0, yend=0), color="black", linetype="solid", size=1, alpha=0.4) +
        geom_vline(data=cnvLine, aes(xintercept=cnvLine$start), color="black", linetype="dotted", size=1, alpha=1)  +
        geom_vline(data=cnvLine, aes(xintercept=cnvLine$end), color="black", linetype="dotted", size=1, alpha=1) +
        geom_segment(data=cnvNator, aes(x=start, xend=end, y=-10, yend=-10), color=cnvNatorColor, size=1, alpha=0.2)  +
        geom_segment(data=manta, aes(x=start, xend=end, y=-20, yend=-20), color=mantaColor, size=1, alpha=1)  +
        theme_bw() + theme(legend.title = element_blank()) +  xlab(cnvLine$chr) + ylab("") # +  scale_colour_manual(values=cbPalette)

ggsave(paste(args[3],".png", sep=""), plot = last_plot(), device = "png", height=5, width=7)
ggsave(paste(args[3],".pdf", sep=""), plot = last_plot(), device = "pdf", height=5, width=7)

# ggplot(tidy_data[tidy_data$mouse %in% mouse_names_covance,], aes(x=start,y=value,color=mouse)) + geom_point(size=1, alpha=0.4) + facet_grid(type ~ . , scales = "free") + geom_segment(data=cnvLine, aes(x=cnvLine$start, xend=cnvLine$end, y=0, yend=0), color="black", linetype="solid", size=3, alpha=0.4) + theme_bw() +  scale_colour_manual(values=cbPalette) + theme(legend.title = element_blank()) +  xlab("Position") + ylab("")

#ggsave(paste(args[3],".covnace.png", sep=""), plot = last_plot(), device = "png")
#ggsave(paste(args[3],".covance.pdf", sep=""), plot = last_plot(), device = "pdf")

#ggplot(tidy_data, aes(x=start,y=value,color=mouse)) + geom_point(size=1, alpha=0.4) + facet_grid(type ~ . , scales = "free") + geom_segment(data=cnvLine, aes(x=cnvLine$start, xend=cnvLine$end, y=0, yend=0), color="black", linetype="solid", size=3, alpha=0.4) + theme_bw() +  scale_colour_manual(values=cbPalette) + theme(legend.title = element_blank()) +  xlab("Position") + ylab("")

#ggsave(paste(args[3],".all.png", sep=""), plot = last_plot(), device = "png")
#ggsave(paste(args[3],".all.pdf", sep=""), plot = last_plot(), device = "pdf")



