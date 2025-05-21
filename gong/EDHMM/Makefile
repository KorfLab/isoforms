CC = gcc

CFLAGS = -fsanitize=address -g
#CFLAGS = -g -Wall -Wextra	# normal with default debug; mostly useless for runtime bug

OBJS = main.o parser.o model.o

all: EDHMM

EDHMM: $(OBJS)
	$(CC) $(CFLAGS) -o EDHMM $(OBJS) -lm
	@chmod +x EDHMM

main.o: main.c parser.c model.c model.h
	@echo "Compiling main"
	$(CC) $(CFLAGS) -c main.c

parser.o: parser.c model.h
	@echo "Compiling function and data structure for computation"
	$(CC) $(CFLAGS) -c parser.c

model.o: model.c model.h
	@echo "Compiling model"
	$(CC) $(CFLAGS) -c model.c

clean:
	@echo "Removing everything but the source file"
	rm -f EDHMM *.o

# with default debug 
debug: EDHMM
	lldb ./EDHMM seq

# with default seq input
test: EDHMM
	./EDHMM seq
	
.PHONY: all clean debug