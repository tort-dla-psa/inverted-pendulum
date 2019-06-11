extra = 
cxx := g++
cxxflags := -ffast-math -std=c++17 $(extra)
incdir = include
srcdir = src
libdir = libs
bindir = bin

IO = $(libdir)/IO.flag

sim_tg:= $(bindir)/sim
sim_src:= $(srcdir)/*.cpp
sim_inc:= -I $(incdir)
sim_inc_files:= $(incdir)/*.h

DEFAULT: $(sim_tg)

prepare:
	mkdir -p $(libdir) $(incdir) $(srcdir) $(bindir)

$(IO):
	$(MAKE) -C libs/UnixIO-cpp/ all
	touch $(IO)

$(sim_tg): $(IO) $(sim_src) $(sim_inc_files)
	$(cxx) $(cxxflags) $(sim_src) $(sim_inc) \
		$(libdir)/UnixIO-cpp/build/file.o \
		$(libdir)/UnixIO-cpp/build/socket_op.o \
		-I $(libdir)/UnixIO-cpp/include \
		-l pthread -lpugixml \
		-o $@

all: prepare $(IO) $(sim_tg)

clean:
	$(MAKE) -C $(libdir)/UnixIO-cpp/ clean
	rm $(IO) $(sim_tg)

