CC = gcc
CFLAGS = -fsanitize=address -g -I include -I .
#CFLAGS = -g -Wall -Wextra -I include -I .

# Parser modules (in parser/)
PARSER_OBJS = parser/parser_seq.o parser/parser_pwm.o parser/parser_mm.o \
              parser/parser_len.o parser/parser_json.o parser/transition.o parser/memory.o

# Algorithm modules (in encoder/)
ALGO_OBJS = encoder/utility.o encoder/forward.o encoder/backward.o encoder/posterior.o encoder/allocate.o

# Decoder modules (in decoder/)
DECODER_OBJS = decoder/viterbi.o decoder/randomf.o

# Core modules (in src/)
CORE_OBJS = main.o output.o cJSON.o

OBJS = $(CORE_OBJS) $(PARSER_OBJS) $(ALGO_OBJS) $(DECODER_OBJS)

# Header files
HEADERS = include/constants.h include/log_math.h include/types.h \
          include/model_structs.h include/algorithm.h include/parser.h \
          include/output.h include/hmm.h model.h

all: hmm

hmm: $(OBJS)
	$(CC) $(CFLAGS) -o hmm $(OBJS) -lm
	@chmod +x hmm

# Core modules
main.o: main.c $(HEADERS) decoder/randomf.h
	@echo "Compiling main"
	$(CC) $(CFLAGS) -c main.c

output.o: output.c $(HEADERS) decoder/randomf.h
	@echo "Compiling output"
	$(CC) $(CFLAGS) -c output.c

# Decoder modules
decoder/viterbi.o: decoder/viterbi.c $(HEADERS) decoder/randomf.h
	@echo "Compiling decoder/viterbi"
	$(CC) $(CFLAGS) -c decoder/viterbi.c -o decoder/viterbi.o

decoder/randomf.o: decoder/randomf.c decoder/randomf.h $(HEADERS)
	@echo "Compiling decoder/randomf"
	$(CC) $(CFLAGS) -c decoder/randomf.c -o decoder/randomf.o

cJSON.o: cJSON.c cJSON.h
	@echo "Compiling cJSON"
	$(CC) $(CFLAGS) -c cJSON.c

# Parser modules
parser/parser_seq.o: parser/parser_seq.c $(HEADERS)
	@echo "Compiling parser/parser_seq"
	$(CC) $(CFLAGS) -c parser/parser_seq.c -o parser/parser_seq.o

parser/parser_pwm.o: parser/parser_pwm.c $(HEADERS)
	@echo "Compiling parser/parser_pwm"
	$(CC) $(CFLAGS) -c parser/parser_pwm.c -o parser/parser_pwm.o

parser/parser_mm.o: parser/parser_mm.c $(HEADERS)
	@echo "Compiling parser/parser_mm"
	$(CC) $(CFLAGS) -c parser/parser_mm.c -o parser/parser_mm.o

parser/parser_len.o: parser/parser_len.c $(HEADERS)
	@echo "Compiling parser/parser_len"
	$(CC) $(CFLAGS) -c parser/parser_len.c -o parser/parser_len.o

parser/parser_json.o: parser/parser_json.c $(HEADERS) cJSON.h
	@echo "Compiling parser/parser_json"
	$(CC) $(CFLAGS) -c parser/parser_json.c -o parser/parser_json.o

parser/transition.o: parser/transition.c $(HEADERS)
	@echo "Compiling parser/transition"
	$(CC) $(CFLAGS) -c parser/transition.c -o parser/transition.o

parser/memory.o: parser/memory.c $(HEADERS)
	@echo "Compiling parser/memory"
	$(CC) $(CFLAGS) -c parser/memory.c -o parser/memory.o

# Algorithm modules
encoder/utility.o: encoder/utility.c $(HEADERS)
	@echo "Compiling encoder/utility"
	$(CC) $(CFLAGS) -c encoder/utility.c -o encoder/utility.o

encoder/forward.o: encoder/forward.c $(HEADERS)
	@echo "Compiling encoder/forward"
	$(CC) $(CFLAGS) -c encoder/forward.c -o encoder/forward.o

encoder/backward.o: encoder/backward.c $(HEADERS)
	@echo "Compiling encoder/backward"
	$(CC) $(CFLAGS) -c encoder/backward.c -o encoder/backward.o

encoder/posterior.o: encoder/posterior.c $(HEADERS)
	@echo "Compiling encoder/posterior"
	$(CC) $(CFLAGS) -c encoder/posterior.c -o encoder/posterior.o

encoder/allocate.o: encoder/allocate.c $(HEADERS)
	@echo "Compiling encoder/allocate"
	$(CC) $(CFLAGS) -c encoder/allocate.c -o encoder/allocate.o

clean:
	@echo "Cleaning build files"
	rm -f hmm *.o *.json parser/*.o encoder/*.o decoder/*.o

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

json-model: hmm
	@echo "Running with JSON model file"
	./hmm --sequence seq --model ../models/worm.splicemodel --verbose

json-rf: hmm
	@echo "Running Random Forest with JSON model"
	./hmm --sequence seq --model ../models/worm.splicemodel --stovit --n_isoforms 100 --json

help: hmm
	./hmm --help

.PHONY: all clean debug test hint rf json-model json-rf help
