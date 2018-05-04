#############################################################################
#
# File:    get_relative_error_stats.r
# Author:  Alex Stivala
# Created: July 2011
#
# get_relative_error_stats.r - get stats on relative errors from
#                              compare_pairwise_with_subsets.py output files
#                              and output as LaTeX table
#
# $Id: get_relative_error_stats.r 778 2011-10-03 04:48:19Z astivala $
#############################################################################


rtabs <- c('compare_individual_with_subsets.rtab',
           'compare_pairwise_oneonly.rtab',
           'compare_pairwise_with_subsets.rtab')
# correpsonds to above files NB also used in if statements, careful!
descriptions <- c('individual only', 'significant pairwise',  'all pairwise')



cat('\\begin{tabular}{lrrrr}\n')
cat('\\hline\n')

cat('Data used &  $\\Delta\\mbox{VHT}$ computations  & Mean Error \\% & Num Error $> 10\\%$ \\\\\n')
#cat('Data used &  $\\Delta\\mbox{VHT}$ computations  & Mean Error \\% & Num Error $> 10\\%$ & Max Error \\%\\\\\n')
cat('\\hline\n')

for (i in 1:length(rtabs)) {
  tab <- read.table(rtabs[i], header=T, stringsAsFactors=F)
  numUpgrades <- length(unique(unlist(strsplit(tab$Subset, ','))))
  numDeltaVHT <- NULL
  if (descriptions[i] == 'individual only') {
    numDeltaVHT <- numUpgrades
  }
  else if (descriptions[i] == 'all pairwise') {
    numDeltaVHT <- numUpgrades*(numUpgrades-1)/2 + numUpgrades
  }
  else {
    # FIXME this is really dodgy, just hardcoding these here, should
    # find some way to get it out of data properly
    if (length(which(unique(unlist(strsplit(tab$Subset, ','))) == "03_03_0101")) == 1) {
      # Chicago Regional data
      numDeltaVHT <- numUpgrades + 1
    }
    else if (length(which(unique(unlist(strsplit(tab$Subset, ','))) == "ber01")) == 1) {
      # Berlin Center data
      numDeltaVHT <- numUpgrades + 11
    }
  }
  cat(descriptions[i], ' & ',  numDeltaVHT  , ' & ',
      format(mean(abs(tab$RelDiff))*100, digits=2,nsmall=1), ' & ')
   cat(length(tab$RelDiff[abs(tab$RelDiff)*100 > 10]),' \\\\\n ')
#  cat(length(tab$RelDiff[abs(tab$RelDiff)*100 > 10]),' & ')
#  cat(format(max(abs(tab$RelDiff))*100, digits=2,nsmall=1), '\\\\\n')

}


rtab <- read.table('compare_pairwise_with_subsets.rtab', header=T, stringsAsFactors=F)
numUpgrades <- length(unique(unlist(strsplit(tab$Subset, ','))))
MAXSUBSETSIZE <- numUpgrades
for (i in 3:MAXSUBSETSIZE) {
  numDeltaVHT <- sum(choose(numUpgrades, 1:i))
  relDiffColName <- paste("RelDiff", i, sep="")
  cat(paste("all subsets size $\\leq$", i), '  & ',
      numDeltaVHT, ' & ',
      format(mean(abs(rtab[[relDiffColName]][!is.na(rtab[[relDiffColName]])]))*100, digits=2,nsmall=1), ' & ')

  cat(length(which(abs(rtab[[relDiffColName]][!is.na(rtab[[relDiffColName]])])*100 > 10)) ,' \\\\\n ')
#  cat(format(max(abs(rtab[[relDiffColName]][!is.na(rtab[[relDiffColName]])]))*100, digits=2,nsmall=1), '\\\\\n')
#  cat(format(max(abs(rtab[[relDiffColName]][!is.na(rtab[[relDiffColName]])]))*100, digits=2,nsmall=1), '\\\\\n')
}

cat('\\hline\n')
cat('\\end{tabular}\n')

