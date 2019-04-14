

g4: g4.o database.o
	gcc -lpq -L/usr/lib -o g4 g4.o database.o -lpq -lm
	
g4.o: g4.c database.h
	gcc -I/usr/include/postgresql -c g4.c
	
database.o: database.c database.h
	gcc -I/usr/include/postgresql -c database.c
	
all: g4

clean: 
	rm *.o
