args <- commandArgs(TRUE)

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

# ggplot(data_melted, aes(x=V2,y=value,color=variable)) + geom_point()
ggplot(tidy_data[tidy_data$mouse %in% mouse_names,], aes(x=start,y=value,color=mouse)) + geom_point(size=0.1, alpha=0.4) + facet_grid(type ~ . , scales = "free")


ggsave(paste(args[3],".png", sep=""), plot = last_plot(), device = "png")
ggsave(paste(args[3],".pdf", sep=""), plot = last_plot(), device = "pdf")


