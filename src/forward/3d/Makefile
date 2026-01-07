#
# © 2025. Triad National Security, LLC. All rights reserved.
#
# This program was produced under U.S. Government contract 89233218CNA000001
# for Los Alamos National Laboratory (LANL), which is operated by
# Triad National Security, LLC for the U.S. Department of Energy/National Nuclear
# Security Administration. All rights in the program are reserved by
# Triad National Security, LLC, and the U.S. Department of Energy/National
# Nuclear Security Administration. The Government is granted for itself and
# others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
# license in this material to reproduce, prepare derivative works,
# distribute copies to the public, perform publicly and display publicly,
# and to permit others to do so.
#
# Author:
#    Kai Gao, kaigao@lanl.gov
#


# paths
bindir = $(HOME)/bin
objdir = $(bindir)/obj/owl_modeling3
moddir = $(bindir)/obj/owl_modeling3

# dependencies
common = \
	module_su.o \
	module_grid.o \
	module_fdcoef.o \
	module_parameters.o \
	module_model.o \
	module_source_receiver.o \
	module_utility.o \
	module_data_processing.o \
	module_anisotropy.o

object1 = \
	interface.o \
	module_vars.o \
	module_cfspml.o \
	module_wavefield.o \
	module_boundary_saving.o \
	submodule_forward_modeling.o

object2 = \
	interface.o \
	module_vars.o \
	module_cfspml.o \
	module_wavefield.o \
	module_boundary_saving.o \
	submodule_forward_modeling.o

object3 = \
	interface.o \
	module_vars.o \
	module_cfspml.o \
	module_wavefield.o \
	module_boundary_saving.o \
	submodule_forward_modeling.o

objectf = \
	module_modeling.o \
	main_forward.o

objc = $(addprefix $(objdir)/c/, $(common))
obj1 = $(addprefix $(objdir)/1/, $(object1))
obj2 = $(addprefix $(objdir)/2/, $(object2))
obj3 = $(addprefix $(objdir)/3/, $(object3))
objf = $(addprefix $(objdir)/f/, $(objectf))

obj = $(objc) $(obj1) $(obj2) $(obj3) $(objf)

# targets
exec = $(bindir)/owl_modeling3

all: makedir $(exec)

# options
# _fdorder*_ can be 2, 4, ..., 14, 16
flitdir = $(HOME)/src/flit

include $(flitdir)/src/Makefile.in
inc = $(base_inc) -I$(flitdir)/lib
lflags = $(flitdir)/lib/libflit.a $(base_lflags)
fflags = $(base_fflags) \
	-D_dim3_ -D_forward_ -D_fdorder16_ \
	# -O0 -g -check bounds -fpe:0 -init=snan

# compile
$(objdir)/c/%.o : $(PWD)/../../lib/%.f90
	$(fc) -o $@ -c $(fflags) $(inc) $<

$(objdir)/1/%.o : $(PWD)/../../solver/acoustic_iso_3d/%.f90
	$(fc) -o $@ -c $(fflags) $(inc) $<

$(objdir)/2/%.o : $(PWD)/../../solver/elastic_vhtiort_3d/%.f90
	$(fc) -o $@ -c $(fflags) $(inc) $<

$(objdir)/3/%.o : $(PWD)/../../solver/elastic_tti_3d/%.f90
	$(fc) -o $@ -c $(fflags) $(inc) $<

$(objdir)/f/%.o : $(PWD)/%.f90
	$(fc) -o $@ -c $(fflags) $(inc) $<

$(objdir)/f/%.o : $(PWD)/../%.f90
	$(fc) -o $@ -c $(fflags) $(inc) $<

# link
$(exec) : $(obj)
	$(fc) -o $@ $^ $(lflags) $(inc)

# make directory
makedir:
	-@mkdir -p $(bindir)
	-@mkdir -p $(objdir)
	-@mkdir -p $(objdir)/c
	-@mkdir -p $(objdir)/f
	-@mkdir -p $(objdir)/1
	-@mkdir -p $(objdir)/2
	-@mkdir -p $(objdir)/3
	-@mkdir -p $(moddir)

# clean
clean:
	-@rm -rf $(objdir) $(moddir)
