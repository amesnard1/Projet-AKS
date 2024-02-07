run: aks
	./aks

aks: aks.c arithmetic.c polynomial.c
	gcc -O3 -lm -lgmp -o aks aks.c arithmetic.c polynomial.c

clean:
	rm aks

