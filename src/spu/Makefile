
# Should we do anything at all?
CELLTOP = 
ifdef CELLTOP

    # Targets
    PROGRAM_spu := runner_spu
    LIBRARY_embed64 := librunner_spu.a

    # Objects
    OBJS = runner_spu.o

    # Local Defines
    IMPORTS = -lm -lmisc -lsimdmath

    # Compiler flags and defines
    CC_OPT_LEVEL := -O3
    CFLAGS := -g -std=gnu99
    CFLAGS := $(CFLAGS) -DVECTORIZE
    CFLAGS := $(CFLAGS) -DSORTED_INTERACTIONS
    # CFLAGS := $(CFLAGS) -DVEC2
    # CFLAGS := $(CFLAGS) -DEXPLICIT_POTENTIALS
    # CFLAGS := $(CFLAGS) -DUSE_DOUBLES
    LDFLAGS := -Wl,-q

    # Footer
    CELL_TOP = $(CELLTOP)/sdk
    include $(CELL_TOP)/buildutils/make.footer
    
else

all:
all-am:
clean:
install:

endif
