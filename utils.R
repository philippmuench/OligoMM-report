# all collected timepoints
timePointsVector <- c(0, 4, 9, 14, 18, 23, 30, 37, 44, 49, 53, 58, 63, 67, 72, 79)

# converts AF to major AF
minorAfToMajorAf <- function(x) if (x > 0) max(1 - x, x) else NA

# day to phase
binDaysByPhase <- function(dat){
	require(dplyr)
	dat <- dat %>% replace(. == 0, "pre-treatment") %>% 
		replace(. == 4, "post-treatment") %>% 
		replace(. == 9, "short-recovery") %>% 
		replace(. == 14, "short-recovery") %>% 
		replace(. == 18, "post-treatment") %>% 
		replace(. == 23, "short-recovery") %>% 
		replace(. == 30, "short-recovery") %>%
		replace(. == 37, "long-recovery") %>% 
		replace(. == 44, "long-recovery") %>% 
		replace(. == 49, "long-recovery") %>% 
		replace(. == 53, "post-treatment") %>% 
		replace(. == 58, "short-recovery") %>% 
		replace(. == 63, "short-recovery") %>% 
		replace(. == 67, "post-treatment") %>% 
		replace(. == 72, "short-recovery") %>% 
		replace(. == 79, "short-recovery")
	return(dat)
}

# day to phase
binDaysByPhaseGroup <- function(dat){
	require(dplyr)
	dat <- dat %>% replace(. == 0, 1) %>% 
		replace(. == 4, 1) %>% 
		replace(. == 9, 1) %>% 
		replace(. == 14, 2) %>% 
		replace(. == 18, 2) %>% 
		replace(. == 23, 3) %>% 
		replace(. == 30, 4) %>%
		replace(. == 37, 1) %>% 
		replace(. == 44, 2) %>% 
		replace(. == 49, 3) %>% 
		replace(. == 53, 3) %>% 
		replace(. == 58, 5) %>% 
		replace(. == 63, 6) %>% 
		replace(. == 67, 3) %>% 
		replace(. == 72, 7) %>% 
		replace(. == 79, 8)
	return(dat)
}

# Mouse ID to type
translateMouseIdToReplicateGroup <- function(dat) {
	require(dplyr)
	dat <- dat %>% replace(. == "1683", "Full") %>% 
		replace(. == "1688", "Full") %>% 
		replace(. == "1692", "Full") %>% 
		replace(. == "1699", "Full") %>% 
		replace(. == "1681", "Replicate 1") %>% 
		replace(. == "1684", "Replicate 2") %>% 
		replace(. == "1686", "Replicate 1") %>% 
		replace(. == "1690", "Replicate 2") %>% 
		replace(. == "1693", "Replicate 1") %>% 
		replace(. == "1694", "Replicate 2") %>% 
		replace(. == "1697", "Replicate 1") %>% 
		replace(. == "1698", "Replicate 2")
	return(dat)
}

# mouse ID to treatment group
translateMouseIdToTreatmentGroup <- function(dat) {
	require(dplyr)
	dat <- dat %>% replace(. == "1683", "Water") %>% 
		replace(. == "1688", "Ciprofloxacin") %>% 
		replace(. == "1692", "Tetracyclin") %>% 
		replace(. == "1699", "Vancomycin") %>% 
		replace(. == "1681", "Water") %>% 
		replace(. == "1684", "Water") %>% 
		replace(. == "1686", "Ciprofloxacin") %>% 
		replace(. == "1690", "Ciprofloxacin") %>% 
		replace(. == "1693", "Tetracyclin") %>% 
		replace(. == "1694", "Tetracyclin") %>% 
		replace(. == "1697", "Vancomycin") %>% 
		replace(. == "1698", "Vancomycin")
	return(dat)
}

# genome ID to genome name
translateGenomeIdToFullName <- function(dat) {
	require(dplyr)
	dat <- dat %>% replace(. == "yl44", "A. muciniphila") %>% 
		replace(. == "i48", "B. caecimuris")  %>% 
		replace(. == "yl27", "M. intestinale")  %>% 
		replace(. == "yl45", "T. muris")  %>% 
		replace(. == "yl2", "B. longum")  %>% 
		replace(. == "kb1", "E. faecalis")  %>% 
		replace(. == "kb18", "A. muris")  %>% 
		replace(. == "yl32", "C. clostridioforme")  %>% 
		replace(. == "yl31", "F. plautii")  %>% 
		replace(. == "yl58", "B. coccoides")  %>% 
		replace(. == "i49", "L. reuteri")  %>% 
		replace(. == "i46", "C. innocuum")
	return(dat)
}


getGene <- function(x){
	pos <- as.numeric(x[9])
	y <- gff[which(with(gff, start <=  pos & end >= pos)),12][1]
	if (length(y) == 0)
		y <- NA
	y
}

getGeneByPosition <- function(x, gff.df,  pos.column = 2, chr.column = 1){
	pos  <- as.numeric(as.matrix(x[pos.column]))
	chr  <- as.character(as.matrix(x[chr.column]))
	gff.chr <- gff.df[which(gff.df$chr == chr),]
	gene <- gff[which(with(gff, start <=  pos & end >= pos)),]$product
	return(gene)
}

# Mouse ID to time group (?) 
# TODO: unclear how to group these
translateMouseIdToClaudiaGroup <- function(dat) {
	require(dplyr)
	dat <- dat %>% replace(. == "1423", "14") %>% 
		replace(. == "1425", "14") %>% 
		replace(. == "1607", "16") %>% 
		replace(. == "1612", "16") %>% 
		replace(. == "1660", "16") %>% 
		replace(. == "1664", "16") %>% 
		replace(. == "1750", "17") %>% 
		replace(. == "1753", "17") %>% 
		replace(. == "1779", "17") %>% 
		replace(. == "1783", "17") %>% 
		replace(. == "1789", "17") %>%
		replace(. == "1801", "18") %>%
		replace(. == "1807", "18") %>%
		replace(. == "1814", "18") %>%
		replace(. == "1815", "18") %>%
		replace(. == "1877", "18") %>%
		replace(. == "1880", "18") %>%
		replace(. == "1881", "18") %>%
		replace(. == "1882", "18") %>%
		replace(. == "1883", "18") %>%
		replace(. == "1884", "18") %>%
		replace(. == "1885", "18") %>%
		replace(. == "I_cc", "I") %>%
		replace(. == "I_mi", "I")
	return(dat)
}

omm_colors = c("yl2" = "#6a51a3","i48" = "#9e9ac8", "yl58" = "#cbc9e2",
yl44 = "#ce1256", "yl27" = "#525252", "yl45" = "#969696", "kb1" = "#969696", "kb18" = "#969696", "yl32" = "#74c476", "i46" = "#238b45","i49" ="#ffffb2", "yl31"="#88419d")

sample_colors = c("stably colonized" = "#74c476", "20 OMM mix" = "#d7b5d8", "40 OMM mix" = "#df65b0", "80 OMM mix" = "#ce1256", "2nd generation" = "#6a51a3", "40 cecal content" = "#41b6c4",  "20 cecal content" = "#a1dab4", "Tag 0, von cecal content" = "#cccccc", "Tag 0, vor 20 tage"=  "#969696" )

variant_colors = c("missense_variant" = "#b2e2e2", "synonymous_variant" = "#66c2a4", "intragenic_variant" = "#fdcc8a", "stop_gained" = "#b30000", "stop_lost&splice_region_variant" = "#252525",  "non_coding_transcript_variant" = "#810f7c")
