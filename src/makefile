dflags = -I/opt/local/include/ -L/opt/local/lib/
cflags = -Wno-deprecated -O3
lflags = -lgsl -lm

runsampler: runsampler.o model.o probdist.o
	g++ $(cflags) $(dflags) $(lflags) runsampler.o model.o probdist.o -o runsampler

runsampler.o: runsampler.c
	g++ $(cflags) $(dflags) -c runsampler.c

model.o: model.c
	g++ $(cflags) $(dflags) -c model.c

probdist.o: probdist.c 
	g++ $(cflags) $(dflags) -c probdist.c
