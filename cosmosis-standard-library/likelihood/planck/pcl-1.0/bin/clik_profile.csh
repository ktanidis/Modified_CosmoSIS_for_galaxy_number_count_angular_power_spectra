# this code cannot be run directly
# do 'source /home/ktanidis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/bin/clik_profile.csh' from your csh shell or put it in your profile

 

if !($?PATH) then
setenv PATH /home/ktanidis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/bin
else
set newvar=$PATH
set newvar=`echo ${newvar} | sed s@:/home/ktanidis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/bin:@:@g`
set newvar=`echo ${newvar} | sed s@:/home/ktanidis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/bin\$@@` 
set newvar=`echo ${newvar} | sed s@^/home/ktanidis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/bin:@@`  
set newvar=/home/ktanidis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/bin:${newvar}                     
setenv PATH /home/ktanidis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/bin:${newvar} 
endif
if !($?PYTHONPATH) then
setenv PYTHONPATH /home/ktanidis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/lib/python2.7/site-packages
else
set newvar=$PYTHONPATH
set newvar=`echo ${newvar} | sed s@:/home/ktanidis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/lib/python2.7/site-packages:@:@g`
set newvar=`echo ${newvar} | sed s@:/home/ktanidis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/lib/python2.7/site-packages\$@@` 
set newvar=`echo ${newvar} | sed s@^/home/ktanidis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/lib/python2.7/site-packages:@@`  
set newvar=/home/ktanidis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/lib/python2.7/site-packages:${newvar}                     
setenv PYTHONPATH /home/ktanidis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/lib/python2.7/site-packages:${newvar} 
endif
if !($?LD_LIBRARY_PATH) then
setenv LD_LIBRARY_PATH 
else
set newvar=$LD_LIBRARY_PATH
set newvar=`echo ${newvar} | sed s@::@:@g`
set newvar=`echo ${newvar} | sed s@:\$@@` 
set newvar=`echo ${newvar} | sed s@^:@@`  
set newvar=:${newvar}                     
setenv LD_LIBRARY_PATH :${newvar} 
endif
if !($?LD_LIBRARY_PATH) then
setenv LD_LIBRARY_PATH /home/ktanidis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/lib
else
set newvar=$LD_LIBRARY_PATH
set newvar=`echo ${newvar} | sed s@:/home/ktanidis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/lib:@:@g`
set newvar=`echo ${newvar} | sed s@:/home/ktanidis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/lib\$@@` 
set newvar=`echo ${newvar} | sed s@^/home/ktanidis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/lib:@@`  
set newvar=/home/ktanidis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/lib:${newvar}                     
setenv LD_LIBRARY_PATH /home/ktanidis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/lib:${newvar} 
endif
if !($?LD_LIBRARY_PATH) then
setenv LD_LIBRARY_PATH 
else
set newvar=$LD_LIBRARY_PATH
set newvar=`echo ${newvar} | sed s@::@:@g`
set newvar=`echo ${newvar} | sed s@:\$@@` 
set newvar=`echo ${newvar} | sed s@^:@@`  
set newvar=:${newvar}                     
setenv LD_LIBRARY_PATH :${newvar} 
endif
if !($?LD_LIBRARY_PATH) then
setenv LD_LIBRARY_PATH /usr1/local/lib
else
set newvar=$LD_LIBRARY_PATH
set newvar=`echo ${newvar} | sed s@:/usr1/local/lib:@:@g`
set newvar=`echo ${newvar} | sed s@:/usr1/local/lib\$@@` 
set newvar=`echo ${newvar} | sed s@^/usr1/local/lib:@@`  
set newvar=:${newvar}                     
setenv LD_LIBRARY_PATH :${newvar} 
endif
setenv CLIK_DATA /home/ktanidis/cosmosis/cosmosis-standard-library/likelihood/planck/plc-1.0/share/clik

setenv CLIK_PLUGIN basic,ffp6_foreground,pep_cib

