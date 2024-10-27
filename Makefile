# Makefile for a single executable

TARGET = main
SRCS = src/Camera.cpp main.cpp src/Trackball.cpp src/imageLoader.cpp src/Mesh.cpp
LIBS = -lglut -lGLU -lGL -lm -lpthread

INCDIR = .
LIBDIR = .
BINDIR = build

CC = g++
CFLAGS = -Wall -O3 -std=c++20
CPPFLAGS = -I$(INCDIR)
LDFLAGS = -L/usr/X11R6/lib
LDLIBS = -L$(LIBDIR) $(LIBS)

OBJS = $(SRCS:%.cpp=$(BINDIR)/%.o)
DEPS = $(SRCS:%.cpp=$(BINDIR)/%.d)

.PHONY: all clean veryclean install installdirs dep

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(OBJS) -o $(TARGET) $(LDFLAGS) $(LDLIBS)

$(BINDIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CC) $(CPPFLAGS) $(CFLAGS) -MMD -c $< -o $@

install: $(TARGET)
	cp $(TARGET) $(BINDIR)/

installdirs:
	test -d $(INCDIR) || mkdir -p $(INCDIR)
	test -d $(LIBDIR) || mkdir -p $(LIBDIR)
	test -d $(BINDIR) || mkdir -p $(BINDIR)

clean:
	rm -f *~ $(TARGET) $(OBJS) $(DEPS)

veryclean: clean
	rm -f $(BINDIR)/$(TARGET)

dep:
	$(CC) $(CPPFLAGS) -MM $(SRCS) > .depend

-include $(DEPS)