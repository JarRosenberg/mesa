include ../make/defaults-module.mk

# Build

MODULE_NAME := turb
SRCS := public/turb_lib.f90 \
        public/turb_def.f90 \
        private/fingering_modes.f90 \
        private/magnetic_diffusion.f90 \
        private/mlt.f90 \
        private/parasite_model_matrices.f90 \
        private/parasite_model.f90 \
        private/tdc_support.f90 \
        private/tdc.f90 \
        private/thermohaline.f90 \
        private/semiconvection.f90 \
        private/velocity_limiter.f90
SRCS_CHECK := test/src/test_turb.f90
INTERNAL_DEPENDS_ON := const utils auto_diff num chem forum
# Test dependency
INTERNAL_DEPENDS_ON += math
EXTERNAL_DEPENDS_ON := lapack95 lapack blas
BINTYPE := lib

# Testing

CHECK_RESULTS_GOLDEN := test/test_output

# Install

MODULES := turb_def.mod turb_lib.mod
INSTALL_INCLUDES :=

include $(MAKE_DIR)/Makefile
