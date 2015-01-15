aprox: main.o splines_2.o points.o aproksymator_fun_try.o 
	$(CC) -o aprox  main.o splines_2.o points.o aproksymator_fun_try.o -Wall -lm

aproksymator_na_bazie.o: makespl.h points.h 
	$(CC) -c aproksymator_fun_try.c

.PHONY: clean

clean:
	-rm *.o aprox 
