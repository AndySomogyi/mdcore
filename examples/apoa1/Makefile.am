
# This file is part of mdcore.
# Coypright (c) 2010 Pedro Gonnet (gonnet@maths.ox.ac.uk)
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Automake stuff
AUTOMAKE_OPTIONS=gnu

# Add the source directory and debug to CFLAGS
AM_CFLAGS = -g -I../../src -Wall $(OPENMP_CFLAGS)

# Set the linker flags if metis is present
AM_LDFLAGS = $(METIS_LDFLAGS) $(FFTW3_LDFLAGS)

# Set-up the library
bin_PROGRAMS = apoa1_noverlet apoa1_pseudoverlet

# Build a CUDA-enabled version too?
if HAVECUDA
    bin_PROGRAMS += apoa1_cuda_noverlet apoa1_cuda_verlet
endif

# Build a MPI-enabled version too?
if HAVEMPI
    bin_PROGRAMS += apoa1_mpi_noverlet apoa1_mpi_pseudoverlet
endif

# Sources for apoa1
apoa1_noverlet_SOURCES = apoa1.c
apoa1_noverlet_CFLAGS =  -DFPTYPE_SINGLE $(AM_CFLAGS)
apoa1_noverlet_LDADD = ../../src/.libs/libmdcore_single.a


# Using pairwise Verlet lists of the second kind
# Sources for apoa1_pseudoverlet
apoa1_pseudoverlet_SOURCES = apoa1.c
apoa1_pseudoverlet_CFLAGS = -DFPTYPE_SINGLE -DENGINE_FLAGS="engine_flag_verlet_pseudo" $(AM_CFLAGS)
apoa1_pseudoverlet_LDADD = ../../src/.libs/libmdcore_single.a


# Using CUDA and Verlet lists
# Sources for apoa1_cuda_verlet
apoa1_cuda_verlet_SOURCES = apoa1.c
apoa1_cuda_verlet_CFLAGS = -DFPTYPE_SINGLE -DWITH_CUDA -DENGINE_FLAGS="(engine_flag_cuda | engine_flag_verlet_pseudo)" $(AM_CFLAGS)
apoa1_cuda_verlet_LDADD = ../../src/.libs/libmdcore_single_cuda.a $(CUDA_LIBS)


# Using CUDA 
# Sources for apoa1_cuda
apoa1_cuda_noverlet_SOURCES = apoa1.c
apoa1_cuda_noverlet_CFLAGS = -DFPTYPE_SINGLE -DWITH_CUDA -DENGINE_FLAGS="engine_flag_cuda" $(AM_CFLAGS)
apoa1_cuda_noverlet_LDADD = ../../src/.libs/libmdcore_single_cuda.a $(CUDA_LIBS)


# Using MPI
# Sources for apoa1
apoa1_mpi_noverlet_SOURCES = apoa1.c
apoa1_mpi_noverlet_CFLAGS =  -DFPTYPE_SINGLE -DWITH_MPI $(AM_CFLAGS)
apoa1_mpi_noverlet_LDADD = ../../src/.libs/libmdcore_single_mpi.a

# Sources for apoa1_pseudoverlet
apoa1_mpi_pseudoverlet_SOURCES = apoa1.c
apoa1_mpi_pseudoverlet_CFLAGS = -DFPTYPE_SINGLE -DWITH_MPI -DENGINE_FLAGS="engine_flag_verlet_pseudo" $(AM_CFLAGS)
apoa1_mpi_pseudoverlet_LDADD = ../../src/.libs/libmdcore_single_mpi.a

