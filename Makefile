

# variables

CFLAGS = -W -Wall -m64 -O3
LFLAGS = -lm -lblas -llapack -lc
COMPILER = clang

# build rules

gaplagbb: main.o instance.o lagrelax.o knapsack.o bb.o minknap.o bundle.o
	$(COMPILER) $(CFLAGS) main.o instance.o lagrelax.o knapsack.o bb.o minknap.o bundle.o -o solve $(LFLAGS) 

main.o: main.c
	$(COMPILER) -c $(CFLAGS) main.c -o main.o 

bb.o: bb.c
	$(COMPILER) -c $(CFLAGS) bb.c -o bb.o

lagrelax.o: lagrelax.c
	$(COMPILER) -c $(CFLAGS) lagrelax.c -o lagrelax.o

knapsack.o: knapsack.c
	$(COMPILER) -c $(CFLAGS) knapsack.c -o knapsack.o 

instance.o: instance.c
	$(COMPILER) -c $(CFLAGS) instance.c -o instance.o 

minknap.o: minknap.c
	$(COMPILER) -c $(CFLAGS) minknap.c -o minknap.o 

bundle.o: bundle.c
	$(COMPILER) -c $(CFLAGS) bundle.c -o bundle.o 

# delete temporary files

clean :
	/bin/rm -rf *.o *~ *.class
	/bin/rm -rf *.mps *.ord *.sos *.lp *.sav *.net *.msg *.log


