
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

# Set the linker flags with any extras
AM_LDFLAGS = $(METIS_LDFLAGS) $(FFTW3_LDFLAGS)

# Set-up the library
bin_PROGRAMS = jac_noverlet jac_pseudoverlet jac_double_pseudoverlet

# Build a CUDA/MPI-enabled version too?
if HAVECUDA
bin_PROGRAMS += jac_cuda_noverlet jac_cuda_verlet
if HAVEMPI
bin_PROGRAMS += jac_mpi_cuda_noverlet
endif
else
if HAVEMPI
bin_PROGRAMS += jac_mpi_noverlet jac_mpi_pseudoverlet
endif
endif

# Sources for jac
jac_noverlet_SOURCES = jac.c
jac_noverlet_CFLAGS =  -DFPTYPE_SINGLE $(AM_CFLAGS)
jac_noverlet_LDADD = ../../src/.libs/libmdcore_single.a


# Using pairwise Verlet lists of the second kind
# Sources for jac_pseudoverlet
jac_pseudoverlet_SOURCES = jac.c
jac_pseudoverlet_CFLAGS = -DFPTYPE_SINGLE -DENGINE_FLAGS="engine_flag_verlet_pseudo" $(AM_CFLAGS)
jac_pseudoverlet_LDADD = ../../src/.libs/libmdcore_single.a

# Sources for jac_double_pseudoverlet
jac_double_pseudoverlet_SOURCES = jac.c
jac_double_pseudoverlet_CFLAGS = -DFPTYPE_DOUBLE -DENGINE_FLAGS="engine_flag_verlet_pseudo" $(AM_CFLAGS)
jac_double_pseudoverlet_LDADD = ../../src/.libs/libmdcore.a


# Using CUDA and Verlet lists
# Sources for jac_cuda_verlet
jac_cuda_verlet_SOURCES = jac.c
jac_cuda_verlet_CFLAGS = -DFPTYPE_SINGLE -DWITH_CUDA -DENGINE_FLAGS="(engine_flag_cuda | engine_flag_verlet_pseudo)" $(AM_CFLAGS)
jac_cuda_verlet_LDADD = ../../src/.libs/libmdcore_single_cuda.a $(CUDA_LIBS)


# Using CUDA 
# Sources for jac_cuda
jac_cuda_noverlet_SOURCES = jac.c
jac_cuda_noverlet_CFLAGS = -DFPTYPE_SINGLE -DWITH_CUDA -DENGINE_FLAGS="engine_flag_cuda" $(AM_CFLAGS)
jac_cuda_noverlet_LDADD = ../../src/.libs/libmdcore_single_cuda.a $(CUDA_LIBS)


# Using MPI
# Sources for jac
jac_mpi_noverlet_SOURCES = jac.c
jac_mpi_noverlet_CFLAGS =  -DFPTYPE_SINGLE -DWITH_MPI $(AM_CFLAGS)
jac_mpi_noverlet_LDADD = ../../src/.libs/libmdcore_single_mpi.a

# Sources for jac_cuda
jac_mpi_cuda_noverlet_SOURCES = jac.c
jac_mpi_cuda_noverlet_CFLAGS =  -DFPTYPE_SINGLE -DWITH_MPI -DWITH_CUDA $(AM_CFLAGS)
jac_mpi_cuda_noverlet_LDADD = ../../src/.libs/libmdcore_single_cuda_mpi.a $(CUDA_LIBS)

# Sources for jac_pseudoverlet
jac_mpi_pseudoverlet_SOURCES = jac.c
jac_mpi_pseudoverlet_CFLAGS = -DFPTYPE_SINGLE -DWITH_MPI -DENGINE_FLAGS="engine_flag_verlet_pseudo" $(AM_CFLAGS)
jac_mpi_pseudoverlet_LDADD = ../../src/.libs/libmdcore_single_mpi.a

