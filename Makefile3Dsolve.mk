CC			=/usr/local/bin/g++-10
# CFLAGS		=-c -Wall -O4 -fopenmp -std=c++17 -I/usr/local/include -I/.#-I/usr/local/opt/icu4c/include #-ffast-math -Ofast -ffinite-math-only
# LDFLAGS		=-fopenmp -std=c++17 -I/usr/local/include -I/.#-L/usr/local/opt/icu4c/lib #-Ofast
CFLAGS		=-c -Wall -O4 -fopenmp -std=c++17 -I../#-I/usr/local/opt/icu4c/include #-ffast-math -Ofast -ffinite-math-only
LDFLAGS		=-fopenmp -std=c++17 -I../#-L/usr/local/opt/icu4c/lib #-Ofast
SOURCES		=./testFMM3Dsolve.cpp
OBJECTS		=$(SOURCES:.cpp=.o)
EXECUTABLE	=./testFMM3Dsolve

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
		$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
		$(CC) $(CFLAGS) $(KERNEL) $(HOMOG) $< -o $@

clean:
	rm a.out testFMM3Dsolve *.o
