===============================
tw-heuristic sequential track !
===============================
To create an executable just do "make".  
tw-heuristic will be created
```
$ make
```
=======================================
Steps to create our experiment results
=======================================

run_script.sh
This is script we have used to test on our machine 

If one runs, this creates a "*_log.csv" file by executing all the testcases 
inside the "test/easy" folder

To run on other folder e.g medium, hard  and random. One needs to edit
the line 4 to suitable folder appropriately. 

==============================
kill_script.sh

Then the kill script should be scheduled in crontab!
by doing crontab -e. 

"/home/rajz/c/pace/bitbucket/tw-heuristic/kill_script.sh" 
should be replaced appropriately before starting run_script.

e.g 
*/1 * * * * /home/rajz/c/pace/bitbucket/tw-heuristic/kill_script.sh >> /home/rajz/c/pace/killlog

The default time is 2 mins for the easy testcases.
While running on other folder medium, hard, random,etc One needs to edit
the line 6 to suitable 5, 59, 10, respectively.

=========================================================================
Remarks
========
1. We ran using the updated test bed. For some test cases our program is
is printing only the s-line and nothing after that. 
2. However if we run our program on those testcases manually it is 
working perfectly as it should be. 
=========================================================================
