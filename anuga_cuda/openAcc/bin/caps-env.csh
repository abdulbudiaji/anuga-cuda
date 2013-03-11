#!/bin/csh

setenv CAPS_HOME /home/devres/tools/CAPSCompilers-3.3.2
setenv HMPP_HOME ${CAPS_HOME}
setenv HMPP_INC_PATH ${HMPP_HOME}/include
setenv HMPP_LIB_PATH ${HMPP_HOME}/lib
setenv HMPP_BIN_PATH ${HMPP_HOME}/bin
setenv HMPP_MAN_PATH ${HMPP_HOME}/man

if (${?PATH} == 1 ) then
  setenv PATH ${HMPP_BIN_PATH}:${HMPP_HOME}/licenses:${PATH}
else
  setenv PATH ${HMPP_BIN_PATH}:${HMPP_HOME}/licenses
endif

if (${?LD_LIBRARY_PATH} == 1 ) then
  setenv LD_LIBRARY_PATH ${HMPP_LIB_PATH}:${LD_LIBRARY_PATH}
else
  setenv LD_LIBRARY_PATH ${HMPP_LIB_PATH}
endif

if (${?MANPATH} == 1 ) then
  setenv MANPATH ${HMPP_MAN_PATH}:${HMPP_HOME}/doc/hmpprt-doxygen/cxx/man:${MANPATH}
else
  setenv MANPATH ${HMPP_MAN_PATH}:${HMPP_HOME}/doc/hmpprt-doxygen/cxx/man:
endif

