#!/bin/csh

 source /afs/cern.ch/sw/lcg/external/gcc/4.6/x86_64-slc6/setup.csh

 setenv ROOTSYS /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.10/x86_64-slc6-gcc46-opt/root

 echo "  "
 echo "You have choosen the next ROOTSYS version:"
 echo "$ROOTSYS"

if( $?LD_LIBRARY_PATH ) then
  setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:$ROOTSYS/lib
else
 echo "No LD_LIBRARY_PATH at all - check $LD_LIBRARY_PATH"
endif

 setenv MyRoot $ROOTSYS/bin/root

 echo " "
 echo "Your ROOT executable is called MyRoot (alias)"
 echo "$MyRoot"
 echo " "

