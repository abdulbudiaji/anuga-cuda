

# Import kernel path and necessary directories
from anuga_cuda.config import archM, codeM, testing_domain_path, \
        merimbula_dir, cairns_dir, kernel_path, kernel_block_configuration,\
        balance_stream, cf_central_stream, extra_1_stream, \
        extra_2_sw_stream, extra_velocity_2_stream, \
        extra_2_edge_swb2_stream, extra_2_limit_vertex_stream, \
        extra_2_limit_edge_stream, evaluate_seg_ref_stream, \
        evaluate_seg_dir_1_stream, evaluate_seg_dir_2_stream, \
        gravity_wb_stream, interpolate_stream, manning_flat_stream, \
        manning_sloped_stream, protect_sw_stream, protect_swb2_stream, \
        saxpy_stream, set_boundary_from_edge_stream, \
        update_centroid_stream, update_stream 


# Import utilities
from utilities import *


# Import GPU_domain Class
from anuga_cuda.gpu_domain_advanced import CUDA_advanced_domain
from anuga_cuda.gpu_domain_basic import CUDA_basic_domain
GPU_domain = CUDA_advanced_domain



