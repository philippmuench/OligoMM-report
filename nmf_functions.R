

# 
#' @export
plotSamples <- function(s, normalize = FALSE, remove.sample.names = F, title = "") {
  h <- data.matrix(s$samples)
  # if(normalize)
  #  h = h / rowSums(h)
  w_df <- reshape2::melt(h, varnames = c("sample", "haplotypes"))
  w_df$haplotypes = factor(w_df$haplotypes)
  
  set.seed(42)
  palette <- randomcoloR::distinctColorPalette(decomposed$nHapotypes)
  
  p <- ggplot2::ggplot(w_df)
  p <- p + ggplot2::geom_bar(aes_string(x = "sample", y = "value", fill = "haplotypes"),
                             color = "black", size = 0.3, stat = "identity", position = "stack")
   p <- p + ggplot2::scale_fill_manual(values = palette) 
  p <- p + ggplot2::coord_flip() + ggplot2::theme_bw()
  p <- p + ggplot2::xlab("") + ggplot2::ylab("Haplotype Contribution")
  p <- p + ggplot2::theme_bw() + theme(panel.border = element_blank(),
                                       panel.grid.major = element_blank(),
                                       panel.grid.minor = element_blank(),
                                       axis.line = element_line(color = "black"))
  p <- p + ggplot2::ggtitle(title)
  p <- p + ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
  if(remove.sample.names)
    p <- p + ggplot2::theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) 
  
  return(p)
}

plotSample <- function(s, normalize = T, sample = "reseq 1696") {
  h <- data.matrix(s$samples)
  if(normalize)
    h = h / rowSums(h)
  w_df <- reshape2::melt(h, varnames = c("sample", "signature"))
  w_df$signature <- factor(w_df$signature)
  w_df <- w_df[which(w_df$sample == sample),]
  
  set.seed(42)
  palette <- randomcoloR::distinctColorPalette(decomposed$nHapotypes)
  
  p <- ggplot(w_df, aes (x = sample, y = value, fill = signature))
  p <- p + geom_bar(size = 0, color = "black", stat = "identity",
                    position = "stack")
  p <- p + theme_bw() + theme(panel.border = element_blank(),
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              axis.line = element_line(color = "black"))
  p <- p + scale_fill_manual(values = palette) + theme_minimal() 
  p <- p + xlab("") + ylab("Haplotype Contribution")
  p <- p + scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
  return(p)
}


plotSamplesByGroup <- function(s, m, normalize = T, title = "") {
  h <- data.matrix(s$samples)
  if(normalize)
    h = h / rowSums(h)
  w_df <- reshape2::melt(h, varnames = c("sample", "signature"))
  w_df$signature <- factor(w_df$signature)
  
  set.seed(42)
  palette <- randomcoloR::distinctColorPalette(decomposed$nHapotypes)
  
  w_df$group <- m[match(w_df$sample, m$sample),]$group
  w_df$day <- m[match(w_df$sample, m$sample),]$day
  w_df$study <- m[match(w_df$sample, m$sample),]$study
  w_df$mouse <- m[match(w_df$sample, m$sample),]$mouse
  
  w_df <- w_df[which(!is.na(w_df$mouse)),] 
  
  p <- ggplot(w_df, aes (x = day, y = value, fill = signature))
  p <- p + geom_bar(size = 0, color = "black", stat = "identity",
                    position = "stack")
  p <- p + facet_wrap(group ~  mouse,  scales = "free", shrink = T, ncol = 3,
                      drop =T)
  p <- p + theme_bw() + theme(panel.border = element_blank(),
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              axis.line = element_line(color = "black"))
  p <- p + geom_vline(xintercept = c(4, 18, 53, 67))
  p <- p + ggtitle(title)
  p <- p + scale_fill_manual(values = palette) + theme_minimal() 
  p <- p + xlab("") + ylab("Haplotype Contribution")
  p <- p + scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
  return(p)
}

plotSamplesByGroup2 <- function(s, m, normalize = T, title = "") {
    library(dplyr)
    h <- data.matrix(s$samples)
    if(normalize)
      h = h / rowSums(h)
    w_df <- reshape2::melt(h, varnames = c("sample", "signature"))
    w_df$signature <- factor(w_df$signature)
    
    set.seed(42)
    palette <- randomcoloR::distinctColorPalette(decomposed$nHapotypes)
    
    w_df$group <- m[match(w_df$sample, m$sample),]$group
    w_df$day <- m[match(w_df$sample, m$sample),]$day
    w_df$study <- m[match(w_df$sample, m$sample),]$study
    w_df$mouse <- m[match(w_df$sample, m$sample),]$mouse
    w_df <- w_df[which(!is.na(w_df$mouse)),] 
    
    w_df2 <- w_df %>% group_by(day, mouse) %>% mutate(Nor = value/sum(value))
    
    p <- ggplot(w_df2, aes (x = day, y = Nor, fill = signature))
    p <- p + geom_area(color = "black", size = 0.1)
    
    p <- p + facet_wrap(group ~  mouse,  scales = "free", shrink = T, ncol = 3,
                        drop =T)
    p <- p + theme_bw() + theme(panel.border = element_blank(),
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                axis.line = element_line(color = "black"))
    p <- p + geom_vline(xintercept = c(4, 18, 53, 67))
    p <- p + ggtitle(title)
    p <- p + scale_fill_manual(values = palette) + theme_minimal() 
    p <- p + xlab("") + ylab("Haplotype Contribution")
    p <- p + scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
    
    return(p)
}

# 
#' @export
plotHaplotypeMap <- function(sigs) {
  df <- data.matrix(sigs$signatures)
  pal <- wesanderson::wes_palette("Zissou1", 21, type = "continuous")
  ha <- ComplexHeatmap::Heatmap(df,
                                name = "contribution",
                                border = TRUE,
                                column_title = "Haplotype contributions",
                                show_row_names = F,
                                cluster_rows = F,
                                show_column_dend = F,
                                show_row_dend = F,
                                cluster_columns = T,
                                col = pal)
  
  return(ha)
}

# 
#' @export
plotNumberHaplotyes <- function(gof) {
  m <- reshape2::melt(gof, id.vars = c("NumberHaplotyes", "Replicate"),
                      measure.vars = c("ExplainedVariance"), variable.name = "stat")
  p <- ggplot2::ggplot(m, ggplot2::aes_string(x = "NumberHaplotyes", y = "value", group = "NumberHaplotyes"))
  p <- p + ggplot2::stat_summary(fun.y = mean, colour = "red", size = 2.5, geom = "point")
  p <- p + ggplot2::geom_point(color = "black", shape = 3)
  p <- p + ggplot2::ylim(0, 1)
  p <- p + ggplot2::facet_wrap(~stat, nrow = 2, scales = "free")
  p <- p + ggplot2::theme_bw() + ggplot2::xlab("Number of Haplotyes") + ggplot2::ylab("explained variance")
  return(p)
}


# 
#' @export
plotSamples <- function(s, normalize = FALSE, remove.sample.names = F, title = "") {
  h <- data.matrix(s$samples)
  if(normalize)
    h = h / rowSums(h)
  w_df <- reshape2::melt(h, varnames = c("sample", "haplotypes"))
  w_df$haplotypes = factor(w_df$haplotypes)
  
  set.seed(42)
  palette <- randomcoloR::distinctColorPalette(decomposed$nHapotypes)
  
  p <- ggplot2::ggplot(w_df)
  p <- p + ggplot2::geom_bar(aes_string(x = "sample", y = "value", fill = "haplotypes"),
                             color = "black", size = 0.3, stat = "identity", position = "stack")
  p <- p + ggplot2::scale_fill_manual(values = palette) 
  p <- p + ggplot2::coord_flip() + ggplot2::theme_bw()
  p <- p + ggplot2::xlab("") + ggplot2::ylab("Haplotype Contribution")
  p <- p + ggplot2::theme_bw() + theme(panel.border = element_blank(),
                                       panel.grid.major = element_blank(),
                                       panel.grid.minor = element_blank(),
                                       axis.line = element_line(color = "black"))
  p <- p + ggplot2::ggtitle(title)
  p <- p + ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
  if(remove.sample.names)
    p <- p + ggplot2::theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) 
  
  return(p)
}

# 
#' @export
plotSample <- function(s, normalize = T, sample = "reseq 1696", title = "") {
  h <- data.matrix(s$samples)
  if(normalize)
    h = h / rowSums(h)
  w_df <- reshape2::melt(h, varnames = c("sample", "signature"))
  w_df$signature <- factor(w_df$signature)
  w_df <- w_df[which(w_df$sample == sample),]
  
  set.seed(42)
  palette <- randomcoloR::distinctColorPalette(decomposed$nHapotypes)
  
  p <- ggplot(w_df, aes (x = sample, y = value, fill = signature))
  p <- p + geom_bar(size = 0, color = "black", stat = "identity",
                    position = "stack")
  p <- p + theme_bw() + theme(panel.border = element_blank(),
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              axis.line = element_line(color = "black"))
  p <- p + scale_fill_manual(values = palette) + theme_minimal() 
  p <- p + xlab("") + ylab("Haplotype Contribution")
  p <- p + ggplot2::ggtitle(title)
  p <- p + scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
  return(p)
}


plotSamplesByTime <- function(s, m, normalize = T, title = "") {
  h <- data.matrix(s$samples)
  if(normalize)
    h = h / rowSums(h)
  w_df <- reshape2::melt(h, varnames = c("sample", "signature"))
  w_df$signature <- factor(w_df$signature)
  
  set.seed(42)
  palette <- randomcoloR::distinctColorPalette(decomposed$nHapotypes)
  
  w_df$group <- m[match(w_df$sample, m$sample),]$group
  w_df$day <- m[match(w_df$sample, m$sample),]$day
  w_df$study <- m[match(w_df$sample, m$sample),]$study
  w_df$mouse <- m[match(w_df$sample, m$sample),]$mouse
  
  w_df <- w_df[which(!is.na(w_df$mouse)),] 
  
  p <- ggplot(w_df, aes (x = day, y = value, fill = signature))
  p <- p + geom_bar(size = 0, color = "black", stat = "identity",
                    position = "stack")
  p <- p + facet_wrap(group ~  mouse,  scales = "free", shrink = T, ncol = 3,
                      drop =T)
  p <- p + theme_bw() + theme(panel.border = element_blank(),
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              axis.line = element_line(color = "black"))
  p <- p + geom_vline(xintercept = c(4, 18, 53, 67))
  p <- p + ggtitle(title)
  p <- p + scale_fill_manual(values = palette) + theme_minimal() 
  p <- p + xlab("") + ylab("Haplotype Contribution")
  p <- p + scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
  return(p)
}


# 
#' @export
plotHaplotypeAnnotation <- function(decomposed, omm_snp_annotation,
                                    hide_hyp = F, sig_threshold = 0.005){
  
  sig <- data.matrix(decomposed$signatures)
  res <- list()
  for (i in seq_along(1:ncol(sig))){
   # message(i)
    one_sig <- as.data.frame(sig[,i])
    one_sig$annotation <- omm_snp_annotation[match(rownames(one_sig),
                                                   omm_snp_annotation$id),]$description
    if (hide_hyp){
      one_sig <- one_sig[which(one_sig$annotation != "outside ORFs" &
                                 one_sig$annotation != "hypothetical protein"),]
    }
    
    colnames(one_sig) <-c("value", "annotation")
    one_sig <- one_sig[order(-one_sig$value),] 
    one_sig$id <- rownames(one_sig)
    one_sig <- one_sig[which(one_sig$value > sig_threshold),]
    one_sig$sig_num <- i
    res[[i]] <- one_sig
  }
  
  all_sig <- do.call(rbind, res)
  set.seed(42)
  palette <- randomcoloR::distinctColorPalette(decomposed$nHapotypes)
  
  p <- ggplot(all_sig, aes(x= reorder(annotation, value), y = value, fill = factor(sig_num)))
  p <- p + geom_bar(stat = "identity")
  p <- p + coord_flip() + theme_minimal() + xlab("") 
  p <- p + scale_fill_manual(values = palette) 
  p <- p + facet_wrap(. ~ sig_num)
  p <- p + theme(text = element_text(size = 6))
  
  return(p)
}