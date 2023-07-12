#!/bin/bash
srun -p queue -t 1 -N 2 -n 96 conway1 1>ConsoleOut.txt 2>>err.log &

echo "Submit successful, you can see file  'err.log' to see your jobId"
echo "console output is in file 'ConsoleOut.txt'"
echo "And use command 'squeue' to see the status of your job"

