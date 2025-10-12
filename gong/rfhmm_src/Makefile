CC = gcc
CFLAGS = -fsanitize=address -g
#CFLAGS = -g -Wall -Wextra	# normal with default debug; mostly useless for runtime bug

OBJS = main.o parser.o model.o output.o viterbi.o randomf.o

all: hmm

hmm: $(OBJS)
	$(CC) $(CFLAGS) -o hmm $(OBJS) -lm
	@chmod +x hmm

main.o: main.c parser.c model.c model.h output.c randomf.h
	@echo "Compiling main"
	$(CC) $(CFLAGS) -c main.c

parser.o: parser.c model.h
	@echo "Compiling function and data structure for computation"
	$(CC) $(CFLAGS) -c parser.c

model.o: model.c model.h
	@echo "Compiling model"
	$(CC) $(CFLAGS) -c model.c

output.o: output.c model.h randomf.h
	@echo "Compiling output session"
	$(CC) $(CFLAGS) -c output.c

viterbi.o: viterbi.c model.h randomf.h
	@echo "Compiling Viterbi algorithm"
	$(CC) $(CFLAGS) -c viterbi.c

randomf.o: randomf.c randomf.h model.h
	@echo "Compiling random forest module"
	$(CC) $(CFLAGS) -c randomf.c

clean:
	@echo "Removing everything but the source file"
	rm -f hmm *.o *.json

debug: hmm
	lldb ./hmm --sequence seq --verbose

test: hmm
	@echo "Running test with debug information (verbose mode)"
	./hmm --sequence seq \
		--don_emission ../models/don.pwm \
		--acc_emission ../models/acc.pwm \
		--exon_emission ../models/exon.mm \
		--intron_emission ../models/intron.mm \
		--ped_exon ../models/exon.len \
		--ped_intron ../models/intron.len \
		--verbose

hint: hmm
	./hmm --sequence seq \
		--don_emission ../models/don.pwm \
		--acc_emission ../models/acc.pwm \
		--exon_emission ../models/exon.mm \
		--intron_emission ../models/intron.mm \
		--ped_exon ../models/exon.len \
		--ped_intron ../models/intron.len \
		--print_splice

# Random Forest Viterbi tests
rf: hmm
	@echo "Running Random Forest isoform generation"
	./hmm --sequence seq \
		--don_emission ../models/don.pwm \
		--acc_emission ../models/acc.pwm \
		--exon_emission ../models/exon.mm \
		--intron_emission ../models/intron.mm \
		--ped_exon ../models/exon.len \
		--ped_intron ../models/intron.len \
		--n_isoforms 100 \
		--m 0.5 \
		--stovit

rf-verbose: hmm
	@echo "Running Random Forest with verbose output"
	./hmm --sequence seq \
		--don_emission ../models/don.pwm \
		--acc_emission ../models/acc.pwm \
		--exon_emission ../models/exon.mm \
		--intron_emission ../models/intron.mm \
		--ped_exon ../models/exon.len \
		--ped_intron ../models/intron.len \
		--n_isoforms 100 \
		--m 0.5 \
		--stovit --verbose

rf-custom: hmm
	@echo "Running Random Forest with custom parameters"
	./hmm --sequence seq \
		--don_emission ../models/don.pwm \
		--acc_emission ../models/acc.pwm \
		--exon_emission ../models/exon.mm \
		--intron_emission ../models/intron.mm \
		--ped_exon ../models/exon.len \
		--ped_intron ../models/intron.len \
		--n_isoforms 10 \
		--m 0.5 \
		--stovit \
		--node_size 5 \
		--json

rf-restrict: hmm
	@echo "Running Random Forest with path restriction"
	./hmm --sequence seq \
		--don_emission ../models/don.pwm \
		--acc_emission ../models/acc.pwm \
		--exon_emission ../models/exon.mm \
		--intron_emission ../models/intron.mm \
		--ped_exon ../models/exon.len \
		--ped_intron ../models/intron.len \
		--n_isoforms 100 \
		--m 0.5 \
		--stovit

rf-short: hmm
	@echo "Running Random Forest (short flags)"
	./hmm -s seq -S -N 500 -m 3 -z 5

test-short: hmm
	@echo "Running test (short flags)"
	./hmm -s seq -v

hint-short: hmm
	@echo "Running hint (short flags)"
	./hmm -s seq -p

help: hmm
	./hmm --help

.PHONY: all clean debug test hint rf rf-verbose rf-custom rf-restrict rf-short test-short hint-short help