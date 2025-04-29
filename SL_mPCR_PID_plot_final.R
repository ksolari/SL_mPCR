setwd("~/Documents/R")
library("PopGenUtils")
library(vcfR)
library(adegenet)
library(dplyr)
library(tidyr)
library(magrittr)
library(tibble)
library(ggplot2)

vcfPK <- read.vcfR( file = "143SNP_WGSn33_mPCR_n5_56_merge_PKonly.recode.vcf", verbose = FALSE )
my_genindPK <- vcfR2genind(vcfPK)
pid_calc(obj = my_genindPK)
pid_permute(my_genindPK, nrep = 1000)
pid_merged

####edited from popgenutils - https://github.com/nikostourvas/PopGenUtils/blob/master/R/pid.R

#' Calculate Probability of Identity (PID) & PIDsibs
#'
#' @param obj a diploid genind object
#'
#' @return a list with four elements reporting
#' PID by locus,
#' PID combinded for all loci,
#' PIDsibs by locus,
#' PIDsibs combined for all loci
#' @export
#'
#' @author Nikolaos Tourvas
#' @description This function calculates PID & PIDsibs by locus, as well as the
#' combined statistic for all loci. The current implementation allows only for
#' diploid data.
#'
#' @examples
#' \dontrun{
#' data(nancycats) # load the "nancycats" genind obj from package adegenet
#' pid <- pid_calc(obj = nancycats)
#' }
#' @import adegenet
#' @import dplyr
#' @import tidyr
#' @importFrom magrittr %>%
#'
pid_calc <- function(obj){
  ## Split obj by locus
  obj_by_loc <- seploc(obj)
  
  ## Calc freqs
  freq_list <- lapply(obj_by_loc, makefreq, quiet=TRUE)
  freq_means_list <- lapply(freq_list, colSums, na.rm=T)
  freq_means_list_perc <- lapply(freq_means_list, function(x){
    x / sum(x) # Values correct - validated via GenAlEx
  })
  
  # Calculate Probability of Identity by locus
  pid_by_loc <- unlist(lapply(freq_means_list_perc, function(x){
    2 * (sum(x^2)^2) - sum(x^4)
  })
  )
  
  # Calculate PIDsibs by locus
  pidsibs_by_loc <- unlist(lapply(freq_means_list_perc, function(x){
    0.25 + (0.5 * sum(x^2)) +
      (0.5 * (sum(x^2)^2) ) -
      (0.25 * sum(x^4))
  })
  )
  
  # Calculate Combined PID & PIDsibs
  pid_comb <- prod(pid_by_loc)
  pidsibs_comb <- prod(pidsibs_by_loc)
  
  return(list("pid_by_loc"=pid_by_loc,
              "pid_comb"=pid_comb,
              "pidsibs_by_loc"=pidsibs_by_loc,
              "pidsibs_comb"=pidsibs_comb))
}



# sample without replacement
# the goal is to create a function similar to genotype_curve
# from poppr package but for pid

#' Calculate Probability of Identity (PID) & PIDsibs for different
#' combinations of loci
#'
#' @param obj a diploid genind object
#' @param nrep the number of permutations to perform
#'
#' @return a list with three elements, reporting
#' the raw data.frame from the permutations performed,
#' the median values of PID & PIDsibs from 1 to the total number of loci,
#' a boxplot plot
#' @export
#'
#' @author Nikolaos Tourvas
#' @examples
#' \dontrun{
#' data(nancycats) # load the "nancycats" genind obj from package adegenet
#' pid_perm <- pid_permute(obj = nancycats)
#' }
#'
#' @import adegenet
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import ggplot2
#' @importFrom magrittr %>%
#pid_permute <- function(obj, nrep=1000)
  
pid_calc <- pid_calc(my_genindPK) # function within the function
  
  maxloc <- nLoc(my_genindPK) - 1
  pid <- list()
  pidsibs <- list()
  for(i in 1:maxloc){
    pid[[i]] <- replicate(1000, prod(sample(pid_calc[[1]], size = i)))
    pidsibs[[i]] <- replicate(1000, prod(sample(pid_calc[[3]], size = i)))
  }
  
  # melt to df
  pid_df <- reshape2::melt(pid)
  pid_df$L1 <- factor(pid_df$L1)
  
  pidsibs_df <- reshape2::melt(pidsibs)
  pidsibs_df$L1 <- factor(pidsibs_df$L1)
  
  pid_merged <- cbind("PID"=pid_df[, 1], "PIDsibs"=pidsibs_df)
  pid_merged <- pid_merged %>%
    rename(PIDsibs=PIDsibs.value,
           loci=PIDsibs.L1)
  
  # make the df tidy!
  pid_merged <- pid_merged %>%
    gather(key = statistic, value = value, -loci)
  
  # calculate some basic statistics
  stats <- pid_merged %>%
    group_by(loci, statistic) %>%
    summarise(median = median(value))
  # spread the data.frame
  stats <- stats %>%
    pivot_wider(names_from = statistic,
                values_from = median)
  # add a line for the the values using all loci
  levels(stats$loci) <- c(levels(stats$loci),
                          nLoc(my_genindPK))
  stats <- add_row(ungroup(stats),
                   loci = as.character(nLoc(my_genindPK)),
                   PID = pid_calc[["pid_comb"]],
                   PIDsibs = pid_calc[["pidsibs_comb"]]
  )
  cbPaletteB = c("PID" = "black","PIDsibs" = "grey")
  
  # boxplot ggplot2
  plot <- ggplot(pid_merged, aes(y=value, x=loci, color=statistic)) +
    scale_color_manual(name="statistic", values=cbPaletteB) +
    geom_boxplot() +
    theme_bw() +
    xlab("Number of loci") +
    ylab("Probability of Identity")
 plot 
 

