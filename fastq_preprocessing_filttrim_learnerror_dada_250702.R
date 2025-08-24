######################
# Initialize
######################
rm(list = ls())


######################
# Library
######################
# BiocManager::install('dada2')
library(dada2)
library(tidyverse)
library(dplyr)
library(Biostrings)
library(ShortRead)
library(phyloseq)
library(ggplot2)
library(ggpubr)
library(vegan)


######################
# Set directory
######################
setwd('/Users/User/Desktop/CRC_restart/')


######################
# function
######################
# is_interleaved_fastq: Check fastq files
is_interleaved_fastq <- function(fastq_file) {
  # Read all lines (compressed .gz files are supported by readLines)
  lines <- readLines(fastq_file)
  
  # FASTQ format: 4 lines per read → extract only header lines
  header_lines <- lines[seq(1, length(lines), 4)]
  
  # Count number of forward and reverse read headers
  r1_count <- sum(grepl(" 1:N:", header_lines))
  r2_count <- sum(grepl(" 2:N:", header_lines))
  
  message("Number of R1 headers: ", r1_count)
  message("Number of R2 headers: ", r2_count)
  
  # Return TRUE only if both R1 and R2 headers are present
  return(r1_count > 0 && r2_count > 0)
}

# quality_profile_png: Save quality profile images
quality_profile_png <- function(fq, dir){
  sample_names <- basename(fq)
  png_file <- file.path(dir, paste0(sample_names, '_QC.png'))
  
  walk2(.x = fq, .y = png_file, ~{
    p <- plotQualityProfile(.x)
    ggsave(filename = .y, plot = p)
  })
}

# mean_quality_by_cycle: Calculate mean quality per cycle for one FASTQ file
mean_quality_by_cycle <- function(fq) {
  pdata <- suppressMessages(plotQualityProfile(fq)$data)
  
  pdata %>%
    group_by(Cycle) %>%
    summarise(avg_quality = weighted.mean(Score, Count)) %>%
    ungroup()
}
# mean_quality_by_cycle <- function(fq) {
#   pdata <- suppressMessages(ShortRead::qa(fq))
#   pdata_percycle_quality <- pdata[['perCycle']]$quality
#   
#   pdata_percycle_quality %>%
#     group_by(Cycle) %>%
#     summarise(avg_quality = weighted.mean(Score, Count)) %>%
#     ungroup()
# }

q_lowq_cumulative_steps <- function(merged_quality, qcut = 30, cycle_cut = 200, plot = TRUE, label_cycles = NULL) {
  # Subset to only cycles >= cycle_cut
  plot_df <- subset(merged_quality, Cycle >= cycle_cut)
  
  # Create a logical matrix: TRUE where quality < qcut, FALSE otherwise
  is_low_q <- plot_df[, -1, drop = FALSE] < qcut
  is_low_q[is.na(is_low_q)] <- FALSE
  
  # For each sample, mark as TRUE once it goes below qcut (cumulative minimum)
  # cummax ensures that once a sample is below qcut, it stays TRUE for all subsequent cycles
  cum_lowq_matrix <- apply(is_low_q, 2, cummax)
  
  # For each cycle, count how many samples have EVER gone below qcut up to that point
  cum_lowq_n <- rowSums(cum_lowq_matrix)
  
  # Assemble a result data.frame with Cycle and cumulative low Q sample count
  res_df <- data.frame(
    Cycle = plot_df$Cycle,
    Cumulative_LowQ_Samples = cum_lowq_n
  )
  
  # For each cycle, get the list of samples that have ever gone below qcut up to that point
  sample_names <- colnames(cum_lowq_matrix)
  cumulative_lowq_lists <- apply(cum_lowq_matrix, 1, function(x) sample_names[as.logical(x)])
  # cumulative_lowq_lists: a list, each element is a character vector of sample names for that cycle
  
  # Find the cycles where the cumulative number of low Q samples increases (i.e., step points)
  step_idx <- c(1, which(diff(res_df$Cumulative_LowQ_Samples) > 0) + 1)
  step_df <- res_df[step_idx, ]
  
  # Plot stepwise change in cumulative low Q samples, if requested
  if (plot) {
    library(ggplot2)
    p <- ggplot(res_df, aes(x = Cycle, y = Cumulative_LowQ_Samples)) +
      geom_step(size = 1.1, direction = "hv") +
      labs(
        title = sprintf("Cumulative samples with Q < %d (Cycle >= %d)", qcut, cycle_cut),
        x = "Cycle",
        y = sprintf("Cumulative samples < Q%d", qcut)
      ) +
      theme_bw()
    
    # If label_cycles is provided, show sample counts as labels on the plot at those cycles
    if (!is.null(label_cycles)) {
      label_points <- res_df[res_df$Cycle %in% label_cycles, ]
      p <- p + geom_text(
        data = label_points,
        aes(label = Cumulative_LowQ_Samples),
        vjust = -0.5, color = "red", size = 5
      )
    }
    print(p)
  }
  
  # Return step points, full summary, and cumulative low Q sample lists per cycle
  return(list(
    step_df = step_df,
    summary_df = res_df,
    cumulative_lowq_lists = cumulative_lowq_lists
  ))
}


######################
# Code
######################
# Read tumor fastq files
fnFs_t <- list.files(path = 'Tumor_rm_primer_fastq/', pattern = '_R1_001_cut.fastq.gz', full.names = T)
fnFs_n <- list.files(path = 'Normal_rm_primer_fastq/', pattern = '_R1_001_cut.fastq.gz', full.names = T)


# Check interleaved status for all files
# Tumor
# interleaved_status_t <- sapply(fnFs_t, is_interleaved_fastq)
# result_t <- data.frame(File = basename(fnFs_t),
#                        Interleaved = interleaved_status_t,
#                        row.names = NULL)
# 
# table(result_t$Interleaved) # FALSE
# 
# Normal
# interleaved_status_n <- sapply(fnFs_n, is_interleaved_fastq)
# result_n <- data.frame(File = basename(fnFs_n),
#                        Interleaved = interleaved_status_n,
#                        row.names = NULL)
# 
# table(result_n$Interleaved) # FALSE


# Visualizing the quality profiles of the forward reads and save
quality_profile_png(fq = fnFs_n, dir = 'Normal_rm_primer_fastq_QC/')
quality_profile_png(fq = fnFs_t, dir = 'Tumor_rm_primer_fastq_QC/')


# Get quality summaries for each file
# quality_profiles_t <- setNames(lapply(fnFs_t, mean_quality_by_cycle), nm = basename(fnFs_t))
# quality_profiles_n <- setNames(lapply(fnFs_n, mean_quality_by_cycle), nm = basename(fnFs_n))
# saveRDS(object = quality_profiles_t, file = 'quality_profiles_tumor_20250702.rds')
# saveRDS(object = quality_profiles_n, file = 'quality_profiles_normal_20250702.rds')
quality_profiles_n <- readRDS('quality_profiles_normal_20250702.rds')
quality_profiles_t <- readRDS('quality_profiles_tumor_20250702.rds')


# Merge all summaries by Cycle position
for (i in seq_along(quality_profiles_t)) {
  colnames(quality_profiles_t[[i]])[2] <- names(quality_profiles_t)[i]
}
for (i in seq_along(quality_profiles_n)) {
  colnames(quality_profiles_n[[i]])[2] <- names(quality_profiles_n)[i]
}
merged_quality_t <- Reduce(function(x, y) full_join(x, y, by = 'Cycle'), quality_profiles_t)
merged_quality_n <- Reduce(function(x, y) full_join(x, y, by = 'Cycle'), quality_profiles_n)

# Visualization of where the largest drop in sample quality occurs
check_t <- q_lowq_cumulative_steps(merged_quality = merged_quality_t, qcut = 30, cycle_cut = 200, label_cycles = c(235, 236))
check_n <- q_lowq_cumulative_steps(merged_quality = merged_quality_n, qcut = 30, cycle_cut = 200, label_cycles = c(239, 240))

remain_samples_t <- setdiff(basename(fnFs_t), check_t$cumulative_lowq_lists[[36]])
remain_samples_t <- file.path(paste0('Tumor_rm_primer_fastq/', remain_samples_t))
remain_samples_n <- setdiff(basename(fnFs_n), check_n$cumulative_lowq_lists[[40]])
remain_samples_n <- file.path(paste0('Normal_rm_primer_fastq/', remain_samples_n))

# filterAndTrim
t_out <- filterAndTrim(fwd = remain_samples_t, filt = file.path(paste0('Tumor_rm_primer_filtrim_fastq/', basename(remain_samples_t))),
                       truncLen = 235, maxN = 0, maxEE = 2, truncQ = 2, rm.phix = T, compress = T, multithread = F) # On Windows set multithread = FALSE

n_out <- filterAndTrim(fwd = remain_samples_n, filt = file.path(paste0('Normal_rm_primer_filtrim_fastq/', basename(remain_samples_n))), 
                       truncLen = 235, maxN = 0, maxEE = 2, truncQ = 2, rm.phix = T, compress = T, multithread = F)


# Learn to error rate
fnFs_t_sel <- list.files(path = 'Tumor_rm_primer_filtrim_fastq/', pattern = 'R1_001_cut.fastq.gz', full.names = T)
fnFs_n_sel <- list.files(path = 'Normal_rm_primer_filtrim_fastq/', pattern = 'R1_001_cut.fastq.gz', full.names = T)

err_t <- learnErrors(fls = fnFs_t_sel, multithread = T, nbases = 1e8)
plotErrors(err_t, nominalQ = T)
err_n <- learnErrors(fls = fnFs_n_sel, multithread = T, nbases = 1e8)
plotErrors(err_n, nominalQ = T)


# dada2
dada_t <- dada(derep = fnFs_t_sel, err = err_t, multithread = T)
names(dada_t) <- paste0('t_', basename(names(dada_t)))
dada_n <- dada(derep = fnFs_n_sel, err = err_n, multithread = T)
names(dada_n) <- paste0('n_', basename(names(dada_n)))


# Save
saveRDS(dada_t, file = 'dada_t_2025_07_02.rds')
saveRDS(dada_n, file = 'dada_n_2025_07_02.rds')
