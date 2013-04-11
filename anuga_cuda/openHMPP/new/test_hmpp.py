import ctypes
import sys

flags = sys.getdlopenflags()
sys.setdlopenflags(flags | ctypes.RTLD_GLOBAL)

from hmpp_python_glue import hmpp_python_test
sys.setdlopenflags(flags)

hmpp_python_test()
