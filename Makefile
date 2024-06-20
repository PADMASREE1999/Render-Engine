CC = g++
CCFLAGS = -g -Wall

rd_view: libcs631.a rd_direct.o pnm_display.o pp_transform.o
	$(CC) -o rd_view $(CCFLAGS) libcs631.a rd_direct.o pnm_display.o pp_transform.o -lm -lX11

# Add whatever additional files and rules here, and also
# in the final linking rule above.

rd_display.o: rd_display.cc rd_display.h
	$(CC) $(CCFLAGS) -c rd_display.cc

rd_direct.o: rd_direct.h rd_direct.cc
	$(CC) $(CCFLAGS) -c rd_direct.cc -o rd_direct.o


pnm_display.o: pnm_display.cc pnm_display.h
	$(CC) $(CCFLAGS) -c pnm_display.cc

pp_transform.o: pp_transform.cc pp_transform.h
	$(CC) $(CCFLAGS) -c pp_transform.cc

clean:
	-rm *.o rd_view

