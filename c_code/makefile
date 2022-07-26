obj = $(wildcard %.o)
%.o:	%.cpp 
	@gcc -c $< 
main:	main.o 
	@gcc *.o -o main	
run: main
	@./main
clean:
	@rm *.o
	@rm main