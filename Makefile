
INITIAL  = disk
HYDRO    = euler_sph
GEOMETRY = spherical
BOUNDARY = polar
RIEMANN  = hllc
TIMESTEP = rk2
OUTPUT   = h5out
RESTART  = h5in

UNAME = $(shell uname)
ifeq ($(UNAME),Linux)
H55 = /home/install/app/hdf5
endif
ifeq ($(UNAME),Darwin)
H55 = /opt/local
endif

CC = mpicc
FLAGS = -O3 -Wall -g

INC = -I$(H55)/include
LIB = -L$(H55)/lib -lm -lhdf5

OBJ = main.o readpar.o onestep.o exchange.o plm.o domain.o faces.o nozzle.o gravity.o misc.o mpisetup.o gridsetup.o $(RIEMANN).o $(TIMESTEP).o $(INITIAL).o $(HYDRO).o $(GEOMETRY).o $(BOUNDARY).o $(OUTPUT).o snapshot.o report.o $(RESTART).o

default: jet

%.o: %.c paul.h
	$(CC) $(FLAGS) $(INC) -c $<

$(RIEMANN).o: Riemann/$(RIEMANN).c paul.h
	$(CC) $(FLAGS) $(INC) -c Riemann/$(RIEMANN).c

$(TIMESTEP).o: Timestep/$(TIMESTEP).c paul.h
	$(CC) $(FLAGS) $(INC) -c Timestep/$(TIMESTEP).c

$(INITIAL).o : Initial/$(INITIAL).c paul.h
	$(CC) $(FLAGS) $(INC) -c Initial/$(INITIAL).c

$(HYDRO).o : Hydro/$(HYDRO).c paul.h
	$(CC) $(FLAGS) $(INC) -c Hydro/$(HYDRO).c

$(GEOMETRY).o : Geometry/$(GEOMETRY).c paul.h
	$(CC) $(FLAGS) $(INC) -c Geometry/$(GEOMETRY).c

$(BOUNDARY).o : Boundary/$(BOUNDARY).c paul.h
	$(CC) $(FLAGS) $(INC) -c Boundary/$(BOUNDARY).c

$(OUTPUT).o : Output/$(OUTPUT).c paul.h
	$(CC) $(FLAGS) $(INC) -c Output/$(OUTPUT).c

$(RESTART).o : Restart/$(RESTART).c paul.h
	$(CC) $(FLAGS) $(INC) -c Restart/$(RESTART).c

jet: $(OBJ) paul.h
	$(CC) $(FLAGS) $(LIB) -o jet $(OBJ)

clean:
	rm -f *.o jet
