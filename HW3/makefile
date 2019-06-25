# Example 8
EXE=HW3

# Main target
all: $(EXE)

#  MinGW
ifeq "$(OS)" "Windows_NT"
CFLG=-O3 -Wall
LIBS=-lglut32cu -lglu32 -lopengl32
CLEAN=del *.exe *.o *.a
else
#  OSX
ifeq "$(shell uname)" "Darwin"
CFLG=-O3 -Wall -Wno-deprecated-declarations
LIBS=-framework GLUT -framework OpenGL
#  Linux/Unix/Solaris
else
CFLG=-O3 -Wall -std=c++11
LIBS=-lglut -lGLU -lGL -lm
endif
#  OSX/Linux/Unix/Solaris
CLEAN=rm -f $(EXE) *.o *.a
endif

# Dependencies
HW3.o: HW3.cpp HW3.h
loadtexbmp.o: loadtexbmp.cpp HW3.h
fatal.o: fatal.cpp HW3.h
errcheck.o: errcheck.cpp HW3.h

#  Create archive
HW3.a:loadtexbmp.o fatal.o errcheck.o
	ar -rcs $@ $^

# Compile rules
.c.o:
	gcc -c $(CFLG) $<
.cpp.o:
	g++ -c $(CFLG) $<

#  Link
HW3:HW3.o HW3.a
	g++ -O3 -o $@ $^   $(LIBS)

#  Clean
clean:
	$(CLEAN)
