

# Import kernel path and necessary directories
from anuga_cuda.config import *


# Import GPU_domain Class
from anuga_cuda.gpu_domain_advanced import CUDA_advanced_domain
from anuga_cuda.gpu_domain_basic import CUDA_basic_domain
GPU_domain = CUDA_advanced_domain


# Import utilities
from utilities.sort_domain import *
from utilities.utility import *


