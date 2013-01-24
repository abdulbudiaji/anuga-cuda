# Import kernel path and necessary directories
from anuga_cuda.config import *


# Import GPU_domain Class
from anuga_cuda.gpu_domain import GPU_domain


# Import utilities
from anuga_cuda.testing_domains.sort_domain import *
from anuga_cuda.testing_domains.utility import *


# Import Channel1 domain
from anuga_cuda.testing_domains.channel1 import generate_domain as c1_domain
# Import Channel3 domain
from anuga_cuda.testing_domains.channel3 import generate_domain as c3_domain
# Import Merimbula domain
from anuga_cuda.testing_domains.merimbula.generate_domain import \
            domain_create as mer_domain
# Import Cairns domain
from anuga_cuda.testing_domains.cairns.runcairns import \
            generate_domain as cairns_domain
