# Import kernel path and necessary directories
from anuga_cuda.config import *


# Import GPU_domain Class
from anuga_cuda.gpu_domain_advanced import GPU_domain


# Import utilities
from anuga_cuda.testing_domains.sort_domain import *
from anuga_cuda.testing_domains.utility import *


# Import Channel1 domain
from anuga_cuda.testing_domains.channel1 import *
# Import Channel3 domain
from anuga_cuda.testing_domains.channel3 import *
# Import Merimbula domain
from anuga_cuda.testing_domains.merimbula.generate_domain import *
# Import Cairns domain
from anuga_cuda.testing_domains.cairns.runcairns import *
