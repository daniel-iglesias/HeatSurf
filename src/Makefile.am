bin_PROGRAMS = beam_dump
beam_dump_SOURCES = beam_dump.cpp cone.cpp cylinder.cpp element.cpp \
	geometry.cpp node.cpp ogive.cpp particle.cpp ring.cpp simulation.cpp twoplates.cpp

# set the include path found by configure
INCLUDES = $(all_includes) -I/usr/include/gsl/ -I/usr/include/vtk-5.6

# the library search path.
beam_dump_LDFLAGS = $(all_libraries) 
noinst_HEADERS = cone.h cylinder.h element.h geometry.h node.h ogive.h \
	particle.h ring.h simulation.h twoplates.h
beam_dump_LDADD = -L/usr/lib/vtk-5.6 -lgsl -lgslcblas -lvtkCommon \
	-lvtksys -lvtkFiltering -lvtkGraphics -lvtkverdict -lvtkHybrid \
	-lvtkImaging  \
	-lvtksqlite -lvtkmetaio -lvtkDICOMParser -lvtkexoIIc -lvtkIO -lvtkRendering \
	-lvtkWidgets 
