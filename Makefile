CC=			gcc
CXX=		g++
CFLAGS=		-g -Wall -O3
CXXFLAGS=	$(CFLAGS)
CPPFLAGS=
INCLUDES=
OBJS=
PROG=		srf
LIBS=		-lz

WFA_ROOT=WFA2-lib

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address
endif

.SUFFIXES:.c .cpp .o
.PHONY:all clean depend

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

.cpp.o:
		$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

srf:$(OBJS) srf.o
		$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

fmd-occ:rld0.o fmd-occ.o
		$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

clean:
		rm -fr gmon.out *.o a.out $(PROG) fmd-occ *~ *.a *.dSYM

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c *.cpp)

# DO NOT DELETE

fmd-occ.o: ketopt.h rld0.h
rld0.o: rld0.h
srf.o: khashl.h ketopt.h kseq.h ksort.h
