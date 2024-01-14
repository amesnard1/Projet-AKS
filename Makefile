run: nogmp
	./nogmp

nogmp: nogmp.c
	gcc -lm -o nogmp nogmp.c

clean:
	rm nogmp

