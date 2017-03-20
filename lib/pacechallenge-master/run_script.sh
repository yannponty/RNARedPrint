#!/bin/bash
#
###################
FOLDER=./test/easy
# without slash at the end! PLEASE!!
################### 

PROGRAM=./tw-heuristic
LOG=/tmp/`date +'%F-%H%M%S'`_log.csv
for grfile in $FOLDER/*.gr;
do
	file="${grfile%%.gr}"
	echo $grfile ">" $file.td
	printf "$file.gr, " >> $LOG
	DT=`date '+%H:%M:%S'`
	printf "$DT, " >> $LOG
	$PROGRAM < $grfile > $file.td	
	if [ $? -eq 0 ] 
	then
		echo "done"
		printf "C," >> $LOG
		
	else
		printf "T," >> $LOG
		echo "terminated"
		
	fi
	LOUTPUT=`head $file.td -n1`
	echo $LOUTPUT
	LOUTPUT=`head $file.td -n1 | awk {'print $4'}` 
	echo $LOUTPUT
	printf "$LOUTPUT, " >> $LOG
	DT=`date '+%H:%M:%S'`
	printf "$DT\n" >> $LOG
	LOUTPUT=0
done


