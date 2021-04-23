#!/bin/bash
# Submit the first job and save the JobID as JOBONE
JOBONE=$(qsub md_VAR_435ns.pbs)
# Submit the second job, use JOBONE as depend, save JobID
JOBTWO=$(qsub -W depend=afterok:$JOBONE md_VAR_460ns.pbs)
# Submit the third job, use JOBTWO as depend, save JobID
JOBTHREE=$(qsub -W depend=afterok:$JOBTWO md_VAR_485ns.pbs)
# Submit the fourth job, use JOBTHREE as depend, save JobID
JOBFOUR=$(qsub -W depend=afterok:$JOBTHREE md_VAR_510ns.pbs)
