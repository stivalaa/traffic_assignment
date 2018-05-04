#!/bin/bash


uname -a >&2
date >&2



echo  "threads = 0" >&2
/usr/bin/time ./tap_frankwolfe -i 10000  ../testdata/ChicagoRegional/ChicagoRegional_net.txt  ../testdata/ChicagoRegional/ChicagoRegional_trips.txt  > /dev/null 

for threads in 1 2 3 4 5 6 7 8
do
  echo  "threads = $threads" >&2
  /usr/bin/time ./tap_frankwolfe_pthread -i 10000 -n $threads ../testdata/ChicagoRegional/ChicagoRegional_net.txt  ../testdata/ChicagoRegional/ChicagoRegional_trips.txt  > /dev/null 
done

times

