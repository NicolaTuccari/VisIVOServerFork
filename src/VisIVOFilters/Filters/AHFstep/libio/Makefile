#**
# Makefile for libio.a, a small library doing input/output stuff for
# AMIGA.
#
# version 0.1 (17.08.2006), Steffen Knollmann (sknollmann@aip.de)
#**

# If ${OPTIMIZE}, ${CCFLAGS} and/or ${DEFINEFLAGS} are set, they will be
# used for the compilation
CFLAGS+=${OPTIMIZE} ${CCFLAGS} ${DEFINEFLAGS}

# If a different MAKE is specified use that one, otherwise try your luck
# with just make; maybe it will work.
MAKE?=make

# The same for ar
AR?=ar

# PHONY targets
.PHONY: all clean

OBJS = io_logging.o \
       io_file.o \
       io_amiga.o \
       io_amiga_header.o \
       io_ares.o \
       io_ares_header.o \
       io_gadget.o \
       io_gadget_header.o \
       io_mgadget.o \
		 io_deva.o\
		 io_deva_header.o\
		 io_tipsy.o\
		 io_tipsy_header.o\
       io_util.o \
       io_parameter.o

# Catch empty make call
all:
	${MAKE} libio.a

# Cleaning up
clean:
	rm -f libio.a $(OBJS)

# The library itself
libio.a: $(OBJS)
	${AR} -r libio.a $(OBJS)

# The logging stuff
io_logging.o: io_logging.h io_logging_defs.h io_logging.c
	$(CC) $(CFLAGS) -c -o io_logging.o io_logging.c

io_file.o: io_file.h io_file_aux.h \
           io_amiga.h io_amiga_def.h \
           io_amiga_header.h io_amiga_header_def.h \
           io_gadget.h io_gadget_def.h \
           io_gadget_header.h io_gadget_header_def.h \
			  io_deva.h io_deva_def.h \
			  io_deva_header.h io_deva_header_def.h \
			  io_tipsy.h io_tipsy_def.h \
			  io_tipsy_header.h io_tipsy_header_def.h \
           io_logging.h \
           io_file.c
	$(CC) $(CFLAGS) -c -o io_file.o io_file.c

io_amiga.o: io_amiga.h \
            io_amiga_header_def.h io_amiga_def.h io_file.h io_file_aux.h \
            io_logging.h io_util.h \
            io_amiga.c
	$(CC) $(CFLAGS) -c -o io_amiga.o io_amiga.c

io_amiga_header.o: io_amiga_header.h \
                   io_amiga_header_def.h io_amiga_def.h io_logging.h \
                   io_util.h \
                   io_amiga_header.c
	$(CC) $(CFLAGS) -c -o io_amiga_header.o io_amiga_header.c

io_ares.o: io_ares.h \
            io_ares_header_def.h io_ares_def.h io_file.h io_file_aux.h \
            io_logging.h io_util.h \
            io_ares.c
	$(CC) $(CFLAGS) -c -o io_ares.o io_ares.c

io_ares_header.o: io_ares_header.h \
                   io_ares_header_def.h io_ares_def.h io_logging.h \
                   io_util.h \
                   io_ares_header.c
	$(CC) $(CFLAGS) -c -o io_ares_header.o io_ares_header.c

io_gadget.o: io_gadget.h \
             io_gadget_header_def.h io_gadget_def.h io_file.h io_file_aux.h \
             io_logging.h io_util.h \
             io_gadget.c
	$(CC) $(CFLAGS) -c -o io_gadget.o io_gadget.c

io_gadget_header.o: io_gadget_header.h \
                    io_gadget_header_def.h io_gadget_def.h io_logging.h \
                    io_util.h \
                    io_gadget_header.c
	$(CC) $(CFLAGS) -c -o io_gadget_header.o io_gadget_header.c

io_mgadget.o: io_mgadget.h \
              io_mgadget_def.h \
              io_gadget_def.h io_gadget_header_def.h io_gadget.h \
              io_file.h io_file_aux.h io_logging.h io_util.h \
              io_mgadget.c
	$(CC) $(CFLAGS) -c -o io_mgadget.o io_mgadget.c

io_deva.o: io_deva.h \
             io_deva_header_def.h io_deva_def.h io_file.h io_file_aux.h \
             io_logging.h io_util.h \
             io_deva.c
	$(CC) $(CFLAGS) -c -o io_deva.o io_deva.c

io_deva_header.o: io_deva_header.h \
                    io_deva_header_def.h io_deva_def.h io_logging.h \
                    io_util.h \
                    io_deva_header.c
	$(CC) $(CFLAGS) -c -o io_deva_header.o io_deva_header.c

io_tipsy.o: io_tipsy.h \
             io_tipsy_header_def.h io_tipsy_def.h io_file.h io_file_aux.h \
             io_logging.h io_util.h \
             io_tipsy.c
	$(CC) $(CFLAGS) -c -o io_tipsy.o io_tipsy.c

io_tipsy_header.o: io_tipsy_header.h \
                    io_tipsy_header_def.h io_tipsy_def.h io_logging.h \
                    io_util.h \
                    io_tipsy_header.c
	$(CC) $(CFLAGS) -c -o io_tipsy_header.o io_tipsy_header.c

io_util.o: io_util.h io_file.h io_util.c
	$(CC) $(CFLAGS) -c -o io_util.o io_util.c

io_paramter.o: io_parameter.h \
               io_parameter_def.h io_util.h \
               io_parameter.c
	$(CC) $(CFLAGS) -c -o io_parameter.o io_paramter.h
