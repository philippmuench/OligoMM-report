# all collected timepoints
timePointsVector <- c(0, 4, 9, 14, 18, 23, 30, 37, 44, 49, 53, 58, 63, 67, 72, 79)

# converts AF to major AF
minorAfToMajorAf <- function(x) if (x > 0) max(1 - x, x) else NA

bugcolors <- list(bug = c("A. muciniphila" = "red", "B. caecimuris"  = "blue", "B. coccoides" = "yellow",
	 "C. clostridioforme"  = "purple", "C. innocuum" = "orange", "F. plautii"= "brown",  "L. reuteri" ="black",
	"M. intestinale" = "grey", "T. muris" = "green"), fixed = c("TRUE"= "red", "FALSE" = "white"))

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
binDaysByLongPhase <- function(dat){
  require(dplyr)
  dat <- dat %>% replace(. == 0, "pre-treatment") %>% 
    replace(. == 4, "after-first-treatment") %>% 
    replace(. == 9, "first-short-recovery") %>% 
    replace(. == 14, "first-short-recovery") %>% 
    replace(. == 18, "after-second-treatment") %>% 
    replace(. == 23, "second-short-recovery") %>% 
    replace(. == 30, "second-short-recovery") %>%
    replace(. == 37, "long-recovery") %>% 
    replace(. == 44, "long-recovery") %>% 
    replace(. == 49, "long-recovery") %>% 
    replace(. == 53, "after-third-treatment") %>% 
    replace(. == 58, "third-short-recovery") %>% 
    replace(. == 63, "third-short-recovery") %>% 
    replace(. == 67, "after-forth-treatment") %>% 
    replace(. == 72, "forth-short-recovery") %>% 
    replace(. == 79, "forth-short-recovery")
  return(dat)
}



# day to phase
binDaysByShortPhase <- function(dat){
  require(dplyr)
  dat <- dat %>% replace(. == 0, FALSE) %>% 
    replace(. == 4, TRUE) %>% 
    replace(. == 9, FALSE) %>% 
    replace(. == 14, FALSE) %>% 
    replace(. == 18, TRUE) %>% 
    replace(. == 23, FALSE) %>% 
    replace(. == 30, FALSE) %>%
    replace(. == 37, FALSE) %>% 
    replace(. == 44, FALSE) %>% 
    replace(. == 49, FALSE) %>% 
    replace(. == 53, TRUE) %>% 
    replace(. == 58, FALSE) %>% 
    replace(. == 63, FALSE) %>% 
    replace(. == 67, TRUE) %>% 
    replace(. == 72, FALSE) %>% 
    replace(. == 79, FALSE)
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
		replace(. == "1698", "Replicate 2") %>% 
	  replace(. == "1682", "qPCR 1") %>% 
	  replace(. == "1685", "qPCR 2") %>% 
	  replace(. == "1687", "qPCR 1") %>% 
	  replace(. == "1689", "qPCR 2") %>% 
	  replace(. == "1691", "qPCR 1") %>% 
	  replace(. == "1695", "qPCR 2") %>% 
	  replace(. == "1696", "qPCR 1") %>% 
	  replace(. == "1700", "qPCR 2")
	return(dat)
}

# mouse ID to treatment group
translateMouseIdToTreatmentGroup <- function(dat) {
	require(dplyr)
	dat <- dat %>% replace(. == "1683", "Water") %>% 
		replace(. == "1691", "Tetracyclin") %>% 
	  replace(. == "1692", "Tetracyclin") %>% 
	  replace(. == "1693", "Tetracyclin") %>% 
	  replace(. == "1694", "Tetracyclin") %>% 
	  replace(. == "1695", "Tetracyclin") %>% 
		replace(. == "1681", "Water") %>% 
	  replace(. == "1682", "Water") %>% 
	  replace(. == "1683", "Water") %>% 
		replace(. == "1684", "Water") %>% 
	  replace(. == "1685", "Water") %>% 
		replace(. == "1686", "Ciprofloxacin") %>% 
		replace(. == "1687", "Ciprofloxacin") %>% 
	  replace(. == "1688", "Ciprofloxacin") %>% 
	  replace(. == "1689", "Ciprofloxacin") %>% 
	  replace(. == "1690", "Ciprofloxacin") %>% 
		replace(. == "1696", "Vancomycin") %>% 
		replace(. == "1697", "Vancomycin") %>% 
	  replace(. == "1698", "Vancomycin") %>% 
	  replace(. == "1699", "Vancomycin") %>% 
		replace(. == "1700", "Vancomycin")
	return(dat)
}

# for reseq samples
translateSampletoMouse <- function(dat) {
	require(dplyr)
	dat <- dat %>% replace(. == "DR1", "1681") %>% 
		replace(. == "DR8", "1691")  %>% 
		replace(. == "DR11", "1696")  %>% 
		replace(. == "DR3", "1681")  %>% 
		replace(. == "DR6", "1687")  %>% 
		replace(. == "DR10", "1691")  %>% 
		replace(. == "DR12", "1696")  %>% 
		replace(. == "DR7", "1687")  %>% 
		replace(. == "DR4", "1687")  %>% 
		replace(. == "DR13", "1696")
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
		replace(. == "ecol", "Mt1B1")  %>% 
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

variant_colors = c("missense_variant" = "#b2e2e2", "synonymous_variant" = "#66c2a4", "intragenic_variant" = "#fdcc8a", "stop_gained" = "#b30000", "stop_lost&splice_region_variant" = "#252525",  "non_coding_transcript_variant" = "#810f7c", "unknown" = "#bdbdbd")

variant_shapes =  c("missense_variant" = 15, "synonymous_variant" = 16, "intragenic_variant" = 17, "stop_gained" = 3, "stop_lost&splice_region_variant" = 8,  "non_coding_transcript_variant" = 6, "unknown" = 12)

orf_shapes =  c("coding" = 15, "non-coding" = 3)



theme_pmuench <- function(base_size = 11, base_family = "") 
{
  half_line <- base_size/2
  theme(
    line = element_line(colour = "black", size = 0.5, 
                        linetype = 1, lineend = "butt"), 
    rect = element_rect(fill = "white", colour = "black",
                        size = 0.5, linetype = 1),
    text = element_text(family = base_family, face = "plain",
                        colour = "black", size = base_size,
                        lineheight = 0.9,  hjust = 0.5,
                        vjust = 0.5, angle = 0, 
                        margin = margin(), debug = FALSE), 
    
    axis.line = element_line(colour = "black", size = 0.5),
    axis.text = element_text(size = rel(0.8), colour = "black"),
    axis.text.x = element_text(margin = margin(t = 0.8*half_line/2), 
                               vjust = 1), 
    axis.text.y = element_text(margin = margin(r = 0.8*half_line/2),
                               hjust = 1),
    axis.ticks = element_line(colour = "black"), 
    axis.ticks.length = unit(half_line/2, "pt"), 
    axis.title.x = element_text(margin = margin(t = 0.8 * half_line,
                                                b = 0.8 * half_line/2)),
    axis.title.y = element_text(angle = 90, 
                                margin = margin(r = 0.8 * half_line,
                                                l = 0.8 * half_line/2)),
    
    legend.background = element_rect(colour = NA), 
    legend.margin = unit(0.2, "cm"), 
    legend.key = element_rect(fill = "white", colour = "white"),
    legend.key.size = unit(1, "lines"), 
    legend.key.height = NULL,
    legend.key.width = NULL, 
    legend.text = element_text(size = rel(0.8)),
    legend.text.align = NULL,
    legend.title = element_text(hjust = 0), 
    legend.title.align = NULL, 
    legend.position = "right", 
    legend.direction = NULL,
    legend.justification = "center", 
    legend.box = NULL, 
    
    panel.background = element_rect(fill = "white", colour = NA),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.margin = unit(half_line, "pt"), panel.margin.x = NULL, 
    panel.margin.y = NULL, panel.ontop = FALSE, 
    
    strip.background = element_rect(colour = "white", fill = "white"),
    strip.text = element_text(colour = "black", face="bold", size = rel(0.8)),
    strip.text.x = element_text(margin = margin(t = half_line,
                                                b = half_line)), 
    strip.text.y = element_text(angle = -90, 
                                margin = margin(l = half_line, 
                                                r = half_line)),
    strip.switch.pad.grid = unit(0.1, "cm"),
    strip.switch.pad.wrap = unit(0.1, "cm"), 
    
    plot.background = element_rect(colour = "white"), 
    plot.title = element_text(size = rel(1.2), 
                              margin = margin(b = half_line * 1.2)),
    plot.margin = margin(half_line, half_line, half_line, half_line),
    complete = TRUE)
}
