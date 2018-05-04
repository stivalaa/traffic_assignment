###############################################################################
#
# plottestspeedup.r - plot speedup graph from mkrtab.sh output
#
# File:    plottestspeedup.r
# Author:  Alex Stivala
# Created: May 2009
#
# Plot speedup graph from mkrtab.sh output
#
# Usage:
#       R --vanilla --slave -f plottestspeedup.r --args rtabfile
#
#    rtabfile is name of 3-colmn (runNum, threads, seconds) .rtab
#      file created by mkrtab.sh
#    output is PostScript file named foo.eps where foo.rtab was the input file
#
# 
# Requires the gplots package from CRAN to draw error bars.
#
# $Id: plotspeedup.r 395 2011-06-22 04:53:26Z astivala $
# 
###############################################################################

library(gplots)

rtabfile <- commandArgs(trailingOnly=TRUE)
timetab <- read.table(rtabfile,header=TRUE)
timetab <- timetab[sort.list(timetab$threads),] # sort by threads ascending

maxthreads = max(timetab$threads)



# time for "0 threads", the baseline (NB not the multithread program on 1
# thread, but a version compiled with no thread code at all)
# (speedup is relative to this )
basetime = mean(subset(timetab, threads==0)$seconds)


# EPS suitable for inserting into LaTeX
postscript(sub('[.]rtab$','.eps',rtabfile),
           onefile=FALSE,paper="special",horizontal=FALSE, 
           width = 9, height = 6)

x <- subset(timetab, runNum==1)$threads
means <- c()
stdev <- c()
ciw <- c()
for (i in x) { 
  # FIXME this code is ugly and inefficient, should use sapply or something
  means <-  c(means, mean(basetime/subset(timetab, threads==i)$seconds))
  thisstdev <-  sqrt(var(basetime/subset(timetab, threads==i)$seconds))
  stdev <-  c(stdev, thisstdev)
  n <- length(subset(timetab, threads==i)$seconds)
  ciw <- c(ciw, qt(0.975, n) * thisstdev / sqrt(n))
}

means
stdev
ciw

plotCI(x, y=means, uiw=ciw, xlab="threads", ylab="speedup" )
lines(x,means)
dev.off()
warnings()
