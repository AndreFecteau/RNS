# define some Makefile variables for the compiler and compiler flags
# to use Makefile variables later in the Makefile: $()
#
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
#
# for C++ define  CC = g++
CC = g++
# CFLAGS  = -Wall -Wextra -DNDEBUG -O3 -fopenmp -isystem /usr/include/eigen3/
CFLAGS  = -Wall -Wextra -DNDEBUG -O3 -fopenmp -isystem /usr/include/eigen3/

# CFLAGS  = -Wall -Wextra -DNDEBUG -O3 -isystem /usr/include/eigen3/

# typing 'make' will invoke the first target entry in the file
# (in this case the default target entry)
# you can name this target entry anything, but "default" or "all"
# are the most commonly used names by convention
#
# default: all

all:
	$(CC) $(CFLAGS) main.cpp -o RNS_Solver

# hello: main.o
# 	g++ main.o -o hello

# To create the executable file count we need the object files
# countwords.o, counter.o, and scanner.o:
#
# all:  main.o
# 	$(CC) $(CFLAGS) -o hello main.o
# # $(RM) all *.o *~

# To create the object file countwords.o, we need the source
# files countwords.c, scanner.h, and counter.h:
#
# main.o:  yeah.cpp
# 	$(CC) $(CFLAGS) -c yeah.cpp
#
# # To start over from scratch, type 'make clean'.  This
# # removes the executable file, as well as old .o object
# # files and *~ backup files:
# #
# clean:
# 	$(RM) all *.o *~
