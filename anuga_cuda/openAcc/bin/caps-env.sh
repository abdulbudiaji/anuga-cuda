#!/bin/sh -e

export CAPS_HOME=/home/devres/tools/CAPSCompilers-3.3.2
export HMPP_HOME=${CAPS_HOME}
export HMPP_INC_PATH=${HMPP_HOME}/include
export HMPP_LIB_PATH=${HMPP_HOME}/lib
export HMPP_BIN_PATH=${HMPP_HOME}/bin
export HMPP_MAN_PATH=${HMPP_HOME}/man

if [ "${PATH:-unset}" = "unset" ]; then
    PATH=${HMPP_BIN_PATH}:${HMPP_HOME}/licenses
else
    PATH=${HMPP_BIN_PATH}:${HMPP_HOME}/licenses:${PATH}
fi
export PATH

if [ "${LD_LIBRARY_PATH:-unset}" = "unset" ]; then
    LD_LIBRARY_PATH=${HMPP_LIB_PATH}
else
    LD_LIBRARY_PATH=${HMPP_LIB_PATH}:${LD_LIBRARY_PATH}
fi
export LD_LIBRARY_PATH

if [ "${MANPATH:-unset}" = "unset" ]; then
    MANPATH=${HMPP_MAN_PATH}:${HMPP_HOME}/doc/hmpprt-doxygen/cxx/man:
else
    MANPATH=${HMPP_MAN_PATH}:${HMPP_HOME}/doc/hmpprt-doxygen/cxx/man:${MANPATH}
fi
export MANPATH
