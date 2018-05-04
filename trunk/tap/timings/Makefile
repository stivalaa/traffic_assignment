##############################################################################
#
# File      : Makefile
# Author    : Alex Stivala (astivala)
# Created   : February 2009
#
# $Id: Makefile 418 2011-06-27 23:34:49Z astivala $
#
# The default target "all" will make graphs of timings (speedups)
# The "clean" target will remove graphs and rtab files
#
# Requires R and the gplots package from CRAN to draw error bars.
# 
##############################################################################


R       = R --vanilla --slave

RTABS = carlton.rtab carlton_10runs.rtab tango_10run.rtab tango.rtab
GRAPHS = $(RTABS:.rtab=.eps)

all: graphs

graphs: $(GRAPHS)

clean:
	rm -f $(RTABS) $(GRAPHS)

%.eps: %.rtab plotspeedup.r 
	$(R) -f plotspeedup.r --args $<

%.rtab: %.out mktimertab.sh
	./mktimertab.sh < $< > $@

