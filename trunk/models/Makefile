##############################################################################
#
# File      : Makefile
# Author    : Alex Stivala (astivala)
# Created   : September 2011
#
# $Id: Makefile 761 2011-09-27 07:08:03Z astivala $
#
##############################################################################


###############################################################################
#
# generate tables used in doc
# 
###############################################################################

all: chicago_heuristic_results.tex chicago_independent_schedule_results.tex \
     chicago_independent_npv_schedule_results.tex \
	 chicago_heuristic_npv_results.tex \
     algorithm_npv_summary.tex

chicago_greedy_output.txt: greedy_heuristic.sh
	./greedy_heuristic.sh > $@

chicago_greedy_npv_output.txt: greedy_heuristic_npv.sh
	./greedy_heuristic_npv.sh > $@

chicago_heuristic_results.tex: chicago_greedy_output.txt ../scripts/schedule_output_to_tex.py
	schedule_output_to_tex.py ../tap/ChicagoRegional_mods.txt < chicago_greedy_output.txt > $@

chicago_heuristic_npv_results.tex: chicago_greedy_npv_output.txt ../scripts/schedule_output_to_tex.py
	schedule_output_to_tex.py ../tap/ChicagoRegional_mods.txt < chicago_greedy_npv_output.txt > $@

chicago_independent_schedule_output.txt: run_chicago_regional_indepdent_eclipse_model.sh
	run_chicago_regional_indepdent_eclipse_model.sh > $@

chicago_independent_schedule_results.tex: chicago_independent_schedule_output.txt ../scripts/schedule_output_to_tex.py
	schedule_output_to_tex.py ../tap/ChicagoRegional_mods.txt < chicago_independent_schedule_output.txt > $@

chicago_independent_npv_schedule_results.tex chicago_independent_npv_schedule_output.txt ../scripts_schedule_output_to_tex.py:
	schedule_output_to_tex.py ../tap/ChicagoRegional_mods.txt < chicago_independent_npv_schedule_output.txt > $@

algorithm_npv_summary.tex: chicago_heuristic_results.tex  chicago_independent_npv_schedule_results.tex	chicago_independent_schedule_results.tex verification/chicago_greedy_schedule_verify.out verification_indpendent_npv/tap_npv_schedule_verify.o1360711 verification_independent/tap_schedule_verify.o1358875 mksummarytab.sh
	./mksummarytab.sh > $@

clean:
	rm chicago_heuristic_results.tex chicago_greedy_output.txt
	rm chicago_independent_schedule_results.tex chicago_indepndent_schedule_output.txt
	rm chicago_indpendent_npv_schedule_results.tex chicago_indpendent_npv_schedule_output_txt
	rm algorithm_npv_summary.tex
