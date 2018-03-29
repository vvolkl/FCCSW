
#!/bin/sh -u
# set FCCSWBASEDIR to the directory containing this script
export FCCSWBASEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
weekday=`date +%a`
source /cvmfs/fcc.cern.ch/sw/views/releases/0.9.1/x86_64-slc6-gcc62-opt/setup.sh
export PYTHIA8_XML=/cvmfs/sft.cern.ch/lcg/releases/MCGenerators/pythia8/226-0ca72/x86_64-slc6-gcc62-opt/share/Pythia8/xmldoc/
export PYTHIA8DATA=$PYTHIA8_XML
