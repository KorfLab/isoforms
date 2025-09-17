CC = gcc
CFLAGS = -fsanitize=address -g
#CFLAGS = -g -Wall -Wextra	# normal with default debug; mostly useless for runtime bug

OBJS = main.o parser.o model.o output.o viterbi.o randomf.o

all: rfhmm

rfhmm: $(OBJS)
	$(CC) $(CFLAGS) -o rfhmm $(OBJS) -lm
	@chmod +x rfhmm

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
	rm -f rfhmm *.o *.json

debug: rfhmm
	lldb ./rfhmm --sequence seq --verbose

test: rfhmm
	@echo "Running test with debug information (verbose mode)"
	./rfhmm --sequence seq \
		--don_emission ../models/don.pwm \
		--acc_emission ../models/acc.pwm \
		--exon_emission ../models/exon.mm \
		--intron_emission ../models/intron.mm \
		--ped_exon ../models/exon.len \
		--ped_intron ../models/intron.len \
		--verbose

hint: rfhmm
	./rfhmm --sequence seq \
		--don_emission ../models/don.pwm \
		--acc_emission ../models/acc.pwm \
		--exon_emission ../models/exon.mm \
		--intron_emission ../models/intron.mm \
		--ped_exon ../models/exon.len \
		--ped_intron ../models/intron.len \
		--print_splice

# Random Forest Viterbi tests
rf: rfhmm
	@echo "Running Random Forest isoform generation"
	./rfhmm --sequence seq \
		--don_emission ../models/don.pwm \
		--acc_emission ../models/acc.pwm \
		--exon_emission ../models/exon.mm \
		--intron_emission ../models/intron.mm \
		--ped_exon ../models/exon.len \
		--ped_intron ../models/intron.len \
		--n_isoforms 100 \
		--m 0.5 \
		--stovit

rf-verbose: rfhmm
	@echo "Running Random Forest with verbose output"
	./rfhmm --sequence seq \
		--don_emission ../models/don.pwm \
		--acc_emission ../models/acc.pwm \
		--exon_emission ../models/exon.mm \
		--intron_emission ../models/intron.mm \
		--ped_exon ../models/exon.len \
		--ped_intron ../models/intron.len \
		--n_isoforms 100 \
		--m 0.5 \
		--stovit --verbose --r

rf-custom: rfhmm
	@echo "Running Random Forest with custom parameters"
	./rfhmm --sequence seq \
		--don_emission ../models/don.pwm \
		--acc_emission ../models/acc.pwm \
		--exon_emission ../models/exon.mm \
		--intron_emission ../models/intron.mm \
		--ped_exon ../models/exon.len \
		--ped_intron ../models/intron.len \
		--n_isoforms 100 \
		--m 0.5 \
		--stovit \
		--node_size 5

rf-restrict: rfhmm
	@echo "Running Random Forest with path restriction"
	./rfhmm --sequence seq \
		--don_emission ../models/don.pwm \
		--acc_emission ../models/acc.pwm \
		--exon_emission ../models/exon.mm \
		--intron_emission ../models/intron.mm \
		--ped_exon ../models/exon.len \
		--ped_intron ../models/intron.len \
		--n_isoforms 100 \
		--m 0.5 \
		--stovit \
		--restrict_path

rf-short: rfhmm
	@echo "Running Random Forest (short flags)"
	./rfhmm -s seq -S -N 500 -m 3 -z 5

test-short: rfhmm
	@echo "Running test (short flags)"
	./rfhmm -s seq -v

hint-short: rfhmm
	@echo "Running hint (short flags)"
	./rfhmm -s seq -p

help: rfhmm
	./rfhmm --help

.PHONY: all clean debug test hint rf rf-verbose rf-custom rf-restrict rf-short test-short hint-short help