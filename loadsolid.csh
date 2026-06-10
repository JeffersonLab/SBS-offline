#clean
module purge
#load standard Hall A software
module use /apps/modulefiles
module use /group/halla/modulefiles
module load analyzer/1.7.12
module list
# load SBS software
setenv SBS_REPLAY /work/halla/solid/efuchey/software/SBS-replay_SOLIDtestbeam/
setenv DB_DIR $SBS_REPLAY/DB
setenv OUT_DIR /volatile/halla/solid/efuchey/replay_test
setenv LOG_DIR /volatile/halla/solid/efuchey/replay_test
setenv ANALYZER_CONFIGPATH $SBS_REPLAY/replay
source /work/halla/solid/efuchey/software/SBS-offline-install/bin/sbsenv.csh
