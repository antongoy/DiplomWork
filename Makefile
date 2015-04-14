all: two_triples one_triple three_triples

three_triples: ranlib.o three_triples.o com.o utils.o
	gcc -Wall -O3 utils.o com.o ranlib.o three_triples.o -o three_triples -llapacke -lm -llapack

two_triples: ranlib.o two_triples.o com.o utils.o
	gcc -Wall -O3 utils.o com.o ranlib.o two_triples.o -o two_triples -llapacke -lm -llapack

one_triple: ranlib.o one_triple.o com.o utils.o
	gcc -Wall -O3 utils.o com.o ranlib.o one_triple.o -o one_triple -llapacke -lm -llapack

ranlib.o: ranlib/ranlib.c
	gcc -O3 -c ranlib/ranlib.c 

com.o: ranlib/com.c
	gcc -O3 -c ranlib/com.c

$FLAGS = -g -ggdb -O0

three_triples.o: three_triples.c
	gcc -Wall -O3 -c three_triples.c

two_triples.o: two_triples.c
	gcc -Wall -O3 -c two_triples.c

one_triple.o: one_triple.c
	gcc -Wall -O3 -c one_triple.c

utils.o: utils.c utils.h
	gcc -Wall -O3 -c utils.c

clean:
	rm -f *.o


