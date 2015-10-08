fabio
export i=$1 #ignore the first $i lines 
export jobFile=$2
export nodeFile=$OAR_NODEFILE
export jobId=${OAR_JOB_ID}
export sshCmd=oarsh

echo "START $i master" $HOSTNAME "jobid" $jobId "jobfile" $jobFile "jobnode" $nodeFile
echo -

export j=0;
declare -A processes;
export ncpu=`eval cat $nodeFile | wc -l`
cpu=(`cat $nodeFile`)

while read cmd
do
   if [ $i -le 0 ]
   then
      if [ 0 -lt $ncpu  ]
      then
         $sshCmd ${cpu[$j]} "export pmc_job_id=${jobId}.${j}; ${cmd}" &
         pid=$!
         echo "cpu=[${cpu[$j]}] pid=[$pid] cmd=[${cmd}]"
         processes[${j}]="${pid}"
         let ncpu-=1
         let j+=1
      fi
   fi
   let i-=1
done < $jobFile

echo -

finished=1
while [ "${finished}" != "0" ]
do
   finished=0
   i=$j
   let i-=1
   while [ $i -ge 0 ] 
   do
      if [ "${processes[${i}]}" != "" ]
      then
         #ps aux | grep ${processes[${i}]}
         let finished+=`ps -p ${processes[${i}]} -o pid= -o comm= | wc -l`
      fi
      let i-=1
   done
   sleep 30
done

echo -
echo "END master" $HOSTNAME "jobid" $jobId "jobfile" $jobFile "jobnode" $nodeFile

