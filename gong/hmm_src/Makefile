CC = gcc
#CFLAGS = -fsanitize=address -g
CFLAGS = -g -Wall -Wextra	# normal with default debug; mostly useless for runtime bug

OBJS = main.o parser.o model.o output.o

all: edhmm2

edhmm2: $(OBJS)
	$(CC) $(CFLAGS) -o edhmm2 $(OBJS) -lm
	@chmod +x edhmm2

main.o: main.c parser.c model.c model.h output.c
	@echo "Compiling main"
	$(CC) $(CFLAGS) -c main.c

parser.o: parser.c model.h
	@echo "Compiling function and data structure for computation"
	$(CC) $(CFLAGS) -c parser.c

model.o: model.c model.h
	@echo "Compiling model"
	$(CC) $(CFLAGS) -c model.c

output.o: output.c model.h
	@echo "Compiling output session"
	$(CC) $(CFLAGS) -c output.c

clean:
	@echo "Removing everything but the source file"
	rm -f edhmm2 *.o

debug: edhmm2
	lldb ./edhmm2 --sequence seq --verbose

test: edhmm2
	@echo "Running test with debug information (verbose mode)"
	./edhmm2 --sequence seq \
		--don_emission ../../models/don.pwm \
		--acc_emission ../../models/acc.pwm \
		--exon_emission ../../models/exon.mm \
		--intron_emission ../../models/intron.mm \
		--ped_exon ../../models/exon.len \
		--ped_intron ../../models/intron.len \
		--verbose

hint: edhmm2
	./edhmm2 --sequence seq \
		--don_emission ../../models/don.pwm \
		--acc_emission ../../models/acc.pwm \
		--exon_emission ../../models/exon.mm \
		--intron_emission ../../models/intron.mm \
		--ped_exon ../../models/exon.len \
		--ped_intron ../../models/intron.len \
		--print_splice

test-short: edhmm2
	@echo "Running test (short flags)"
	./edhmm2 -s seq -v

hint-short: edhmm2  
	@echo "Running hint (short flags)"
	./edhmm2 -s seq -p

help: edhmm2
	./edhmm2 --help

.PHONY: all clean debug test hint test-short hint-short help