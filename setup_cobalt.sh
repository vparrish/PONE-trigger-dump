# /bin/bash

# setting environment variables
HOME=$(pwd)

source $HOME/src/envs/tga/bin/activate
if [ ! -d "build" ]; then
    mkdir build
fi
if [ ! -d "outdir" ]; then
    mkdir outdir
fi
# install simweights if not installedcd
if ! python -c "import simweights" >/dev/null 2>&1; then
    pip install simweights
fi


eval $(/cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/setup.sh)
/cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/RHEL_7_x86_64/metaprojects/icetray/v1.8.2/env-shell.sh 

