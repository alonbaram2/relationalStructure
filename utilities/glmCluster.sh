
#!/bin/bash

sub=$1
glm=$2
roi=$3
queue=$4
iBlock=$5

scriptsDir=/home/fs0/abaram/scripts/NN2020/utilities
cat ${scriptsDir}/clusterWrapper.m | sed "s/sub-00/${sub}/g" | sed "s/GLM1/${glm}/g" | sed "s/surf/${roi}/g" | sed "s/999/${iBlock}/g"> ${scriptsDir}/${sub}_${glm}_${roi}_${iBlock}.m;
jobID1=`fsl_sub -N ${sub}${iBlock}${glm}${roi}${queue} -q ${queue}.q matlab -nodisplay -nosplash \< ${scriptsDir}/${sub}_${glm}_${roi}_${iBlock}.m`;
echo $jobID1
