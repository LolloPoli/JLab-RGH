#!/bin/bash
cd $RGH_SIM_WORK_DIR
i=1

export NITERATIONS=5000
export NEVENTS=10000

while [ $i -le $NITERATIONS ]
do
echo "$i > $PWD/rgh_submits_list/submit$i.sh"
echo
cp job.sh rgh_jobs_list/job$i.sh
cp submit.sh rgh_submits_list/submit$i.sh
sed -i "s;MCINDEX=0;MCINDEX=$i;g" rgh_jobs_list/job$i.sh
sed -i "s;NEVENTS=100;NEVENTS=$NEVENTS;g" rgh_jobs_list/job$i.sh
sed -i "s;\$RGH_SIM_HOME/job.sh;bash \$RGH_SIM_HOME/rgh_jobs_list/job$i.sh;g" rgh_submits_list/submit$i.sh
sbatch rgh_submits_list/submit$i.sh
((i++))
done
