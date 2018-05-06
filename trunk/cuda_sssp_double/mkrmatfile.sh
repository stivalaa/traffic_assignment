#!/bin/sh

# make the graph matrix file from Chicago Regional flow output. The cost on each edge
# is the link cost on that link from TAP solution so is what the TAP is doing
# at each iteration (when it has solved to this precision)

gunzip -c ../tap/results/ChicagoRegional_flows.txt.gz | ../scripts/flow2rmat.py  > rmat.txt.chicagoregionalflows


