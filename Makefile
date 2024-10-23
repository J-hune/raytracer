# Makefile pour un unique exécutable

CIBLE = main
SRCS = src/Camera.cpp main.cpp src/Trackball.cpp src/imageLoader.cpp src/Mesh.cpp
LIBS = -lglut -lGLU -lGL -lm -lpthread

INCDIR = .
LIBDIR = .
BINDIR = .

CC = g++
CPP = g++

CFLAGS = -Wall -O3
CXXFLAGS = -Wall -O3
CPPFLAGS = -I$(INCDIR)
LDFLAGS = -L/usr/X11R6/lib
LDLIBS = -L$(LIBDIR) $(LIBS)

OBJS = $(SRCS:.cpp=.o)

.PHONY: all clean veryclean install installdirs dep

$(CIBLE): $(OBJS)
	$(CC) $(OBJS) -o $(CIBLE) $(LDFLAGS) $(LDLIBS)

install: $(CIBLE)
	cp $(CIBLE) $(BINDIR)/

installdirs:
	test -d $(INCDIR) || mkdir -p $(INCDIR)
	test -d $(LIBDIR) || mkdir -p $(LIBDIR)
	test -d $(BINDIR) || mkdir -p $(BINDIR)

clean:
	rm -f *~ $(CIBLE) $(OBJS)

veryclean: clean
	rm -f $(BINDIR)/$(CIBLE)

dep:
	g++ $(CPPFLAGS) -MM $(SRCS)

all: dep $(CIBLE)

# liste des dépendances générée par 'make dep'
Camera.o: src/Camera.cpp src/Camera.h src/Vec3.h src/Trackball.h
main.o: main.cpp src/Vec3.h src/Camera.h src/Trackball.h
Trackball.o: src/Trackball.cpp src/Trackball.h
