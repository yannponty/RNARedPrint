#!/bin/sh
for pid in `ps -ef |  grep tw-heuristic | awk {'print $2'}`
do
	TIME=`ps -ef | grep tw-heuristic | awk -F':' {'print $3'} | head -1`
	echo $TIME
	if [ $TIME -ge 5 ]   
	then
		kill -15 $pid
		echo "killing $pid running time $TIME"
	fi
done
echo `date`

