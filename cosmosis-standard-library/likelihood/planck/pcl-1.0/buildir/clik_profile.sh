# this code cannot be run directly
# do 'source /home/ktanidis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/bin/clik_profile.sh' from your sh shell or put it in your profile

function addvar () {
local tmp="${!1}" ;
tmp="${tmp//:${2}:/:}" ; tmp="${tmp/#${2}:/}" ; tmp="${tmp/%:${2}/}" ;
export $1="${2}:${tmp}" ;
} 

if [ -z "${PATH}" ]; then 
PATH=/home/ktanidis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/bin
export PATH
else
addvar PATH /home/ktanidis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/bin
fi
if [ -z "${PYTHONPATH}" ]; then 
PYTHONPATH=/home/ktanidis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/lib/python2.7/site-packages
export PYTHONPATH
else
addvar PYTHONPATH /home/ktanidis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/lib/python2.7/site-packages
fi
if [ -z "${LD_LIBRARY_PATH}" ]; then 
LD_LIBRARY_PATH=
export LD_LIBRARY_PATH
else
addvar LD_LIBRARY_PATH 
fi
if [ -z "${LD_LIBRARY_PATH}" ]; then 
LD_LIBRARY_PATH=/home/ktanidis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/lib
export LD_LIBRARY_PATH
else
addvar LD_LIBRARY_PATH /home/ktanidis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/lib
fi
if [ -z "${LD_LIBRARY_PATH}" ]; then 
LD_LIBRARY_PATH=
export LD_LIBRARY_PATH
else
addvar LD_LIBRARY_PATH 
fi
if [ -z "${LD_LIBRARY_PATH}" ]; then 
LD_LIBRARY_PATH=/usr1/local/lib
export LD_LIBRARY_PATH
else
addvar LD_LIBRARY_PATH /usr1/local/lib
fi
CLIK_DATA=/home/ktanidis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/share/clik
export CLIK_DATA

CLIK_PLUGIN=basic,ffp6_foreground,pep_cib
export CLIK_PLUGIN

