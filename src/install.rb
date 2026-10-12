#
# © 2025-2026. Triad National Security, LLC. All rights reserved.
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


#
# Build the CPU executables of OWL:
#     ruby install.rb                  build
#     ruby install.rb clean            rebuild from scratch
# Options:
#     fdorder=2, 4, ..., 16            order of the finite-difference stencils (default 16)
#     hdf5io=yes or no                 HDF5 I/O of models, wavefield snapshots and inversion models (default no)
#     flitdir=<directory>              the FLIT directory (default ~/src/flit)
# e.g. ruby install.rb clean fdorder=8 hdf5io=yes. The options take effect on all files only with clean.
#

fdorder = "16"
hdf5io = "no"
flitdir = File.join(Dir.home, "src", "flit")
ARGV.each do |arg|
	key, value = arg.split("=", 2)
	case key
	when "clean"
	when "fdorder"
		fdorder = value
	when "hdf5io"
		hdf5io = value
	when "flitdir"
		flitdir = File.expand_path(value.to_s)
	else
		abort "install.rb: unknown option #{arg}"
	end
end
abort "install.rb: fdorder must be one of 2, 4, ..., 16" unless (2..16).step(2).map(&:to_s).include?(fdorder)
abort "install.rb: hdf5io must be yes or no" unless ["yes", "no"].include?(hdf5io)
abort "install.rb: #{flitdir}/src/Makefile.in not found (flitdir)" unless File.file?(File.join(flitdir, "src", "Makefile.in"))
opts = "fdorder=#{fdorder} hdf5io=#{hdf5io} flitdir='#{flitdir}'"

for dir in ["forward/2d",
            "forward/3d",
            "inverse/2d",
            "inverse/3d"]

	if ARGV.include?("clean")

		system(" cd #{dir}; make clean; make #{opts} &")

	else

		system(" cd #{dir}; make #{opts} & ")

	end

end
