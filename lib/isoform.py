
import copy
import gzip
import itertools
import json
import math
import random
import sys

def anti(seq):
	comp = str.maketrans('ACGTRYMKWSBDHV', 'TGCAYRKMWSVHDB')
	anti = seq.translate(comp)[::-1]
	return anti

GCODE = {
	'AAA' : 'K',	'AAC' : 'N',	'AAG' : 'K',	'AAT' : 'N',
	'AAR' : 'K',	'AAY' : 'N',	'ACA' : 'T',	'ACC' : 'T',
	'ACG' : 'T',	'ACT' : 'T',	'ACR' : 'T',	'ACY' : 'T',
	'ACK' : 'T',	'ACM' : 'T',	'ACW' : 'T',	'ACS' : 'T',
	'ACB' : 'T',	'ACD' : 'T',	'ACH' : 'T',	'ACV' : 'T',
	'ACN' : 'T',	'AGA' : 'R',	'AGC' : 'S',	'AGG' : 'R',
	'AGT' : 'S',	'AGR' : 'R',	'AGY' : 'S',	'ATA' : 'I',
	'ATC' : 'I',	'ATG' : 'M',	'ATT' : 'I',	'ATY' : 'I',
	'ATM' : 'I',	'ATW' : 'I',	'ATH' : 'I',	'CAA' : 'Q',
	'CAC' : 'H',	'CAG' : 'Q',	'CAT' : 'H',	'CAR' : 'Q',
	'CAY' : 'H',	'CCA' : 'P',	'CCC' : 'P',	'CCG' : 'P',
	'CCT' : 'P',	'CCR' : 'P',	'CCY' : 'P',	'CCK' : 'P',
	'CCM' : 'P',	'CCW' : 'P',	'CCS' : 'P',	'CCB' : 'P',
	'CCD' : 'P',	'CCH' : 'P',	'CCV' : 'P',	'CCN' : 'P',
	'CGA' : 'R',	'CGC' : 'R',	'CGG' : 'R',	'CGT' : 'R',
	'CGR' : 'R',	'CGY' : 'R',	'CGK' : 'R',	'CGM' : 'R',
	'CGW' : 'R',	'CGS' : 'R',	'CGB' : 'R',	'CGD' : 'R',
	'CGH' : 'R',	'CGV' : 'R',	'CGN' : 'R',	'CTA' : 'L',
	'CTC' : 'L',	'CTG' : 'L',	'CTT' : 'L',	'CTR' : 'L',
	'CTY' : 'L',	'CTK' : 'L',	'CTM' : 'L',	'CTW' : 'L',
	'CTS' : 'L',	'CTB' : 'L',	'CTD' : 'L',	'CTH' : 'L',
	'CTV' : 'L',	'CTN' : 'L',	'GAA' : 'E',	'GAC' : 'D',
	'GAG' : 'E',	'GAT' : 'D',	'GAR' : 'E',	'GAY' : 'D',
	'GCA' : 'A',	'GCC' : 'A',	'GCG' : 'A',	'GCT' : 'A',
	'GCR' : 'A',	'GCY' : 'A',	'GCK' : 'A',	'GCM' : 'A',
	'GCW' : 'A',	'GCS' : 'A',	'GCB' : 'A',	'GCD' : 'A',
	'GCH' : 'A',	'GCV' : 'A',	'GCN' : 'A',	'GGA' : 'G',
	'GGC' : 'G',	'GGG' : 'G',	'GGT' : 'G',	'GGR' : 'G',
	'GGY' : 'G',	'GGK' : 'G',	'GGM' : 'G',	'GGW' : 'G',
	'GGS' : 'G',	'GGB' : 'G',	'GGD' : 'G',	'GGH' : 'G',
	'GGV' : 'G',	'GGN' : 'G',	'GTA' : 'V',	'GTC' : 'V',
	'GTG' : 'V',	'GTT' : 'V',	'GTR' : 'V',	'GTY' : 'V',
	'GTK' : 'V',	'GTM' : 'V',	'GTW' : 'V',	'GTS' : 'V',
	'GTB' : 'V',	'GTD' : 'V',	'GTH' : 'V',	'GTV' : 'V',
	'GTN' : 'V',	'TAA' : '*',	'TAC' : 'Y',	'TAG' : '*',
	'TAT' : 'Y',	'TAR' : '*',	'TAY' : 'Y',	'TCA' : 'S',
	'TCC' : 'S',	'TCG' : 'S',	'TCT' : 'S',	'TCR' : 'S',
	'TCY' : 'S',	'TCK' : 'S',	'TCM' : 'S',	'TCW' : 'S',
	'TCS' : 'S',	'TCB' : 'S',	'TCD' : 'S',	'TCH' : 'S',
	'TCV' : 'S',	'TCN' : 'S',	'TGA' : '*',	'TGC' : 'C',
	'TGG' : 'W',	'TGT' : 'C',	'TGY' : 'C',	'TTA' : 'L',
	'TTC' : 'F',	'TTG' : 'L',	'TTT' : 'F',	'TTR' : 'L',
	'TTY' : 'F',	'TRA' : '*',	'YTA' : 'L',	'YTG' : 'L',
	'YTR' : 'L',	'MGA' : 'R',	'MGG' : 'R',	'MGR' : 'R',
}

def translate_str(seq):
	pro = []
	for i in range(0, len(seq), 3):
		codon = seq[i:i+3]
		if codon in GCODE: pro.append(GCODE[codon])
		else: pro.append('X')
	return ''.join(pro)

def longest_orf(seq):
	orfs = []
	for f in range(3):
		pro = translate_str(seq[f:])
		start = 0
		while start < len(pro):
			stop = 0
			if pro[start] == 'M':
				for i, s in enumerate(pro[start+1:]):
					if s == '*':
						stop = i + start +1
						break
			if stop != 0:
				orfs.append( (pro[start:stop], start*3 + f) )
			start += 1
	orfs.sort(key=lambda t: len(t[0]), reverse=True)

	print(orfs)
	sys.exit('testing')

	if len(orfs) > 0: return orfs[0]
	else:             return (None, None)

#####################
## UTILITY SECTION ##
#####################

def randseq(n):
	seq = ''
	for i in range(n):
		seq += random.choice('ACGT')
	return seq

def get_filepointer(filename):
	fp = None
	if   filename.endswith('.gz'): fp = gzip.open(filename, 'rt')
	elif filename == '-':          fp = sys.stdin
	else:                          fp = open(filename)
	return fp

def read_fasta(filename):

	name = None
	seqs = []

	fp = get_filepointer(filename)

	for line in fp:
		line = line.rstrip()
		if line.startswith('>'):
			if len(seqs) > 0:
				seq = ''.join(seqs)
				yield(name, seq)
				name = line[1:]
				seqs = []
			else:
				name = line[1:]
		else:
			seqs.append(line)
	yield(name, ''.join(seqs))
	fp.close()

def prob2score(p):
	if p == 0: return -100
	return math.log2(p/0.25)

def entropy(ps):
	assert(math.isclose(sum(ps), 1))
	h = 0
	for p in ps:
		h -= p * math.log2(p)
	return h

def dkl(p, q):
	assert(math.isclose(sum(p), 1.0))
	assert(math.isclose(sum(q), 1.0))
	d = 0
	for pi, qi in zip(p, q):
		d += pi * math.log2(pi / qi)
	return d

def manhattan(p, q):
	assert(math.isclose(sum(p), 1.0, abs_tol=1e-6))
	assert(math.isclose(sum(q), 1.0, abs_tol=1e-6))
	d = 0
	for pi, qi in zip(p, q):
		d += abs(pi - qi)
	return d

def intersection(p, q): # aka dtc / 2
	assert(math.isclose(sum(p), 1.0, abs_tol=1e-6))
	assert(math.isclose(sum(q), 1.0, abs_tol=1e-6))
	d = 0
	for pi, qi in zip(p, q):
		d += min(pi, qi)
	return 1 - d

def d1(p, q):
	assert(math.isclose(sum(p), 1.0, abs_tol=1e-6))
	assert(math.isclose(sum(q), 1.0, abs_tol=1e-6))
	d = 0
	for pi, qi in zip(p, q):
		if abs(pi - qi) > d: d = abs(pi - qi)
	return d

#########################
## SPLICEMODEL SECTION ##
#########################

def read_splicemodel(file):
	with open(file) as fp: model = json.load(fp)
	return model

#################
## PWM SECTION ##
#################

def create_pwm(seqs):
	count = []
	for seq in seqs:
		for i, nt in enumerate(seq):
			if len(count) <= i:
				count.append({'A':0, 'C': 0, 'G': 0, 'T': 0})
			count[i][nt] += 1

	pwm = [{} for i in range(len(count))]
	for i in range(len(count)):
		for nt in count[i]:
			pwm[i][nt] = count[i][nt] / len(seqs)
	return pwm

def write_pwm(file, pwm):
	with open(file, 'w') as fp:
		fp.write(f'% PWM {file} {len(pwm)}\n')
		for pos in pwm:
			for nt in pos:
				fp.write(f'{pos[nt]:.6f} ')
			fp.write('\n')

def read_pwm(file):
	nts = ('A', 'C', 'G', 'T')
	pwm = []

	# read raw values
	with open(file) as fp:
		for line in fp:
			if line.startswith('%'): continue
			f = line.split()
			d = {}
			for nt, val in zip(nts, f):
				d[nt] = float(val)
			pwm.append(d)

	# convert to log-odds
	for i in range(len(pwm)):
		for nt in nts:
			pwm[i][nt] = prob2score(pwm[i][nt])

	return pwm

def score_pwm(pwm, seq, pos, memo=None):
	if memo is not None and pos in memo: return memo[pos]
	score = 0
	for i in range(len(pwm)):
		nt = seq[pos+i]
		score += pwm[i][nt]
	if memo is not None: memo[pos] = score
	return score

####################
## LENGTH SECTION ##
####################

def find_tail(val, x): ## this shouldn't be needed anymore...
	lo = 0
	hi = 1000

	while hi - lo > 1:
		m = (hi + lo) / 2;
		p = 1 / m
		f = (1 - p)**(x-1) * p
		if f < val: lo += (m - lo) / 2
		else:       hi -= (hi - m) / 2

	return m

def write_len(file, hist):
	with open(file, 'w') as fp:
		fp.write(f'% LEN {file} {len(hist)}\n')
		for val in hist:
			fp.write(f'{val:.6f}\n')

def read_len(file):
	model = []
	with open(file) as fp:
		for line in fp:
			if line.startswith('%'): continue
			line = line.rstrip()
			model.append(float(line))

	size = len(model)
	tail = find_tail(model[-1], size)
	expect = 1 / size;
	for i in range(len(model)):
		if model[i] == 0: model[i] = -100
		else:             model[i] = math.log2(model[i] / expect)

	return {'tail': tail, 'size':len(model), 'val': model}

def score_len(model, x, memo=None):
	if x >= model['size']:
		p = 1 / model['tail']
		q = (1-p)**(x-1) * p
		expect = 1 / model['size']
		s = math.log2(q/expect)
		return s
	else:
		return model['val'][x]

##########################
## MARKOV MODEL SECTION ##
##########################

def create_markov(seqs, order, beg, end):
	count = {}
	for seq in seqs:
		for i in range(beg+order, len(seq) - end):
			ctx = seq[i-order:i]
			nt = seq[i]
			if ctx not in count: count[ctx] = {'A':0, 'C':0, 'G':0, 'T':0}
			count[ctx][nt] += 1

	# these need to be probabilities
	mm = {}
	for kmer in count:
		mm[kmer] = {}
		total = 0
		for nt in count[kmer]: total += count[kmer][nt]
		for nt in count[kmer]: mm[kmer][nt] = count[kmer][nt] / total

	return mm

def write_markov(file, mm):
	with open(file, 'w') as fp:
		fp.write(f'% MM {file} {len(mm)*4}\n')
		for kmer in sorted(mm):
			for v in mm[kmer]:
				fp.write(f'{kmer}{v} {mm[kmer][v]:.6f}\n')
			fp.write('\n')

def read_markov(file):
	mm = {}
	k = None
	with open(file) as fp:
		for line in fp:
			if line.startswith('%'): continue
			f = line.split()
			if len(f) == 2:
				mm[f[0]] = prob2score(float(f[1]))
				if k == None: k = len(f[0])
	return {'k': k, 'mm': mm}

def score_markov(model, seq, beg, end, memo=None):
	if memo is not None and (beg,end) in memo: return memo[(beg,end)]
	score = 0
	k = model['k']
	mm = model['mm']

	for i in range(beg, end -k + 2):
		kmer = seq[i:i+k]
		score += mm[kmer]

	if memo is not None: memo[(beg,end)] = score
	return score


################################
## ISOFORM GENERATION SECTION ##
################################

def gtag_sites(seq, flank, minex):
	dons = []
	accs = []
	for i in range(flank + minex, len(seq) -flank -minex):
		if seq[i:i+2]   == 'GT': dons.append(i)
		if seq[i-1:i+1] == 'AG': accs.append(i)
	return dons, accs

def gff_sites(seq, gff, gtag=True):
	dond = {}
	accd = {}
	with open(gff) as fp:
		for line in fp:
			if line.startswith('#'): continue
			f = line.split()
			if len(f) < 8: continue
			if f[2] != 'intron': continue
			if f[6] != '+': continue
			beg = int(f[3]) -1
			end = int(f[4]) -1
			if gtag:
				if seq[beg:beg+2]   != 'GT': continue
				if seq[end-1:end+1] != 'AG': continue
			dond[beg] = True
			accd[end] = True

	dons = list(dond.keys())
	accs = list(accd.keys())
	dons.sort()
	accs.sort()

	return dons, accs

class Isoform:
	"""Class to represent a single isoform"""

	def __init__(self, seq, beg, end, dons, accs, model=None, weights=None,
			memo=None):
		assert(beg <= end)
		assert(len(dons) == len(accs))

		self.seq = seq
		self.beg = beg
		self.end = end
		self.dons = dons
		self.accs = accs
		self.exons = []
		self.introns = []
		self.model = model
		self.weights = weights
		self.memo = memo
		self.score = None
		self.prob = None   # set externally: via Locus
		self.txseq = ''
		self.rnaidx = []

		# translation specific tiggered by self.translate(pos)
		self.start_codon = None
		self.stop_codon = None
		self.cdsseq = None
		self.aaseq = None
		self.splice_after_stop = None # will be a distance
		self.utr3_length = None # will be a distance
		self.rnatype = None # coding, non-coding, non-stop, nmd-target

		# create exons, introns, and txseq
		if len(dons) == 0:
			self.exons.append((beg, end))
		else:
			for a, b in zip(dons, accs): self.introns.append((a, b))
			self.exons.append((beg, dons[0] -1))
			for i in range(1, len(dons)):
				a = accs[i-1] + 1
				b = dons[i] -1
				self.exons.append((a, b))
			self.exons.append((accs[-1] +1, end))
		for beg, end in self.exons: self.txseq += seq[beg:end+1]

		if model is not None: self.compute_score()

	def compute_score(self, reweight=None):

		if self.model is None:
			self.score = None
			return

		weights = reweight if reweight else self.weights

		seq = self.seq
		acc = self.model['acc']
		don = self.model['don']
		exs = self.model['exs']
		ins = self.model['ins']
		exl = self.model['exl']
		inl = self.model['inl']
		inf = self.model['inf']
		wacc = weights['wacc'] if weights is not None else 1
		wdon = weights['wdon'] if weights is not None else 1
		wexs = weights['wexs'] if weights is not None else 1
		wins = weights['wins'] if weights is not None else 1
		wexl = weights['wexl'] if weights is not None else 1
		winl = weights['winl'] if weights is not None else 1
		winf = weights['winf'] if weights is not None else 1
		macc = self.memo['acc'] if self.memo is not None else None
		mdon = self.memo['don'] if self.memo is not None else None
		mexs = self.memo['exs'] if self.memo is not None else None
		mins = self.memo['ins'] if self.memo is not None else None

		s = 0
		for i in self.accs:
			s += score_pwm(acc, seq, i -len(acc)+1, memo=macc) * wacc
		for i in self.dons:
			s += score_pwm(don, seq, i, memo=mdon) * wdon
		for b, e in self.exons:
			s += score_len(exl, e - b + 1) * wexl
		for b, e in self.introns:
			s += score_len(inl, e - b + 1) * winl
		for b, e in self.exons:
			s += score_markov(exs, seq, b, e, memo=mexs) * wexs
		for b, e in self.introns: s += score_markov(ins, seq,
			b + len(don), e - len(acc), memo=mins) * wins
		s += inf * len(self.introns) * winf

		self.score = s

	def _canonical_start_codon(self, x):
		pos = set()
		pos.add(self.dnaidx_to_rnaidx(x))
		pos.add(self.dnaidx_to_rnaidx(x+1))
		pos.add(self.dnaidx_to_rnaidx(x+2))
		if None in pos: return False

		aug = self.dnaidx_to_rnaidx(x)
		start_codon = self.txseq[aug:aug+3]
		if start_codon != 'ATG': return False

		return True

	def translate(self, dna_start, minejc=10, minutr=300):
		# find the start codon: canonical, alternate, or non-coding
		if self._canonical_start_codon(dna_start):
			self.start_codon = dna_start                # canonical
		else:
			x = self.txseq.find('ATG')
			if x == -1:
				self.rnatype = 'non-coding'             # non-coding
				return
			self.start_codon = self.rnaidx_to_dnaidx(x) # alternate

		# translate to find the stop codon
		self.stop_codon = None
		aas = []
		cds = []
		x = self.dnaidx_to_rnaidx(self.start_codon)
		while True:
			c0 = self.rnaidx_to_dnaidx(x)
			c1 = self.rnaidx_to_dnaidx(x+1)
			c2 = self.rnaidx_to_dnaidx(x+2)
			if c0 is None: break
			if c1 is None: break
			if c2 is None: break
			codon = self.seq[c0] + self.seq[c1] + self.seq[c2]
			cds.append(codon)
			aa = GCODE[codon]
			aas.append(aa)
			if codon == 'TAA' or codon == 'TAG' or codon == 'TGA':
				self.stop_codon = self.rnaidx_to_dnaidx(x + 2)
				break
			x += 3
		self.aaseq = ''.join(aas)
		self.cdsseq = ''.join(cds)

		# check for non-stop
		if self.stop_codon is None:
			self.rnatype = 'non-stop'
			return

		# find maximum distance of EJC/splice after stop codon
		self.splice_after_stop = -1
		for beg, end in self.introns:
			d = beg - self.stop_codon
			if d > self.splice_after_stop:
				self.splice_after_stop = d
				break
		self.utr3_length = len(self.seq) - self.stop_codon
		if self.splice_after_stop > minejc or self.utr3_length > minutr:
			self.rnatype = 'nmd-target'
			return

		# made it!
		self.rnatype = 'coding'

	def _compute_txcoords(self):
		rna = ['-'] * len(self.seq)
		for beg, end in self.exons:
			for i in range(beg, end+1): rna[i] = self.seq[i]
		tseq = ''.join(rna)
		txcoor = 0
		tpos = []
		for nt in tseq:
			if nt == '-':
				tpos.append(None)
			else:
				tpos.append(txcoor)
				txcoor += 1
		self.rnaidx = tpos

	def dnaidx_to_rnaidx(self, dnaidx):
		if len(self.rnaidx) == 0: self._compute_txcoords()
		return self.rnaidx[dnaidx]

	def rnaidx_to_dnaidx(self, rnaidx):
		if len(self.rnaidx) == 0: self._compute_txcoords()
		for i, pos in enumerate(self.rnaidx):
			if pos is not None and pos == rnaidx: return i
		return None

class Locus:
	"""Class to represent an alternatively spliced locus"""

	def __init__(self, desc, seq, model, constraints=None, weights=None,
			gff=None, limit=None, countonly=False, memoize=False):

		# sequence stuff
		self.desc = desc
		self.name = desc.split()[0]
		self.seq = seq

		# model, weights, and memoizers pass-through to transcript
		self.model = model
		self.weights = weights
		if memoize: self.memo = {'acc': {}, 'don': {}, 'exs': {}, 'ins': {}}
		else: self.memo = None

		# constraints
		if constraints is None:
			self.imin = 35
			self.emin = 25
			self.flank = 100
		else:
			self.imin = constraints['min_intron']
			self.emin = constraints['min_exon']
			self.flank = constraints['flank']

		# algorithm init
		if gff: self.dons, self.accs = gff_sites(seq, gff)
		else:   self.dons, self.accs = gtag_sites(seq, self.flank, self.emin)
		self.isoforms = []
		self.worst = None
		self.rejected = 0
		self.limit = limit
		if self.limit: self.resize = self.limit * 2
		else:          self.resize = None
		self.countonly = countonly
		self.isocount = 0

		# recursion
		introns = []
		for i in range(len(self.dons)):
			if self.countonly and self.limit and self.isocount >= self.limit: return
			self._build_isoforms(self.dons[i:], self.accs, introns)

		# finalization
		if self.limit:
			x = sorted(self.isoforms, key=lambda d: d.score, reverse=True)
			self.isoforms = x[:self.limit]
		self.calculate_isoform_probabilities()

	def calculate_isoform_probabilities(self):
		weight = []
		total = 0
		for tx in self.isoforms:
			w = 2 ** tx.score
			weight.append(w)
			total += w
		prob = [w / total for w in weight]
		for p, tx in zip(prob, self.isoforms): tx.prob = p
		x = sorted(self.isoforms, key=lambda d: d.score, reverse=True)
		self.isoforms = x

	def _build_isoforms(self, dons, accs, introns):
		if self.countonly and self.limit and self.isocount >= self.limit: return
		don = dons[0]
		for aix, acc in enumerate(accs):
			if acc - don + 1 < self.imin: continue
			intron = (don, acc)

			# legit isoform as is, save it
			if self.countonly and self.limit and self.isocount >= self.limit: return
			iso = copy.copy(introns)
			iso.append(intron)
			self._save_isoform(iso)

			# also extend it
			for dix, ndon in enumerate(dons):
				elen = ndon - acc -1
				if elen >= self.emin:
					ext = copy.copy(iso)
					self._build_isoforms(dons[dix:], accs[aix:], ext)

	def _save_isoform(self, introns):

		# counting only?
		if self.countonly:
			self.isocount += 1
			return

		# create transcript
		dsites = [x[0] for x in introns]
		asites = [x[1] for x in introns]
		tx = Isoform(self.seq, self.flank, len(self.seq) - self.flank -1,
			dsites, asites, model=self.model, weights=self.weights,
			memo=self.memo)

		# unlimited?
		if self.limit is None:
			self.isoforms.append(tx)
			return

		 # don't save low-scoring isoforms
		if self.worst is not None and tx.score < self.worst:
			self.rejected += 1
			return

		# store (sort and prune as necessary)
		self.isoforms.append(tx)
		if len(self.isoforms) > self.resize:
			x = sorted(self.isoforms, key=lambda d: d.score, reverse=True)
			before = len(x)
			self.isoforms = x[:self.limit]
			diff = before - len(self.isoforms)
			self.rejected += diff
			self.worst = self.isoforms[-1].score

	def write_gff(self, fp):
		"""write locus as GFF to stream"""

		print('# name:', self.name, file=fp)
		print('# length:', len(self.seq), file=fp)
		print('# donors:', len(self.dons), file=fp)
		print('# acceptors:', len(self.accs), file=fp)
		print('# isoforms:', len(self.isoforms), file=fp)
		print('# rejected:', self.rejected, file=fp)
		print(f'# maxprob: {self.isoforms[0].prob:.4g}', file=fp)
		print(f'# minprob {self.isoforms[-1].prob:.4g}', file=fp)
		print(f'# complexity: {complexity(self.isoforms):.3f}', file=fp)
		print(f'# model: {self.mod_name}')

		src = 'apc'
		cs = f'{self.name}\t{src}\t'
		b = self.flank + 1
		e = len(self.seq) - self.flank
		gene = f'Gene-{self.name}'
		print(f'{cs}gene\t{b}\t{e}\t.\t+\t.\tID={gene}\n', file=fp)
		for i, tx in enumerate(self.isoforms):
			b = tx.beg + 1
			e = tx.end + 1
			s = tx.prob
			tid = f'tx-{self.name}-{i+1}'
			print(f'{cs}mRNA\t{b}\t{e}\t{s:.4g}\t+\t.\tID={tid};Parent={gene}',
				file=fp)
			p = f'Parent={tid}'
			for b, e in tx.exons:
				print(f'{cs}exon\t{b+1}\t{e+1}\t{s:.4g}\t+\t.\t{p}', file=fp)
			for b, e in tx.introns:
				print(f'{cs}intron\t{b+1}\t{e+1}\t{s:.4g}\t+\t.\t{p}', file=fp)
			print(file=fp)


########################
## EXPRESSION SECTION ##
########################

def complexity(txs):
	prob = []
	total = 0
	for tx in txs:
		w = 2 ** tx.score
		prob.append(w)
		total += w
	for i in range(len(prob)):
		prob[i] /= total
	return entropy(prob)

def get_introns(gff):
	introns = {}
	total = 0
	with open(gff) as fp:
		for line in fp:
			if line.startswith('#'): continue
			f = line.split()
			if len(f) < 8: continue
			if f[2] != 'intron': continue
			if f[6] != '+': continue
			beg, end, score = f[3], f[4], f[5]
			beg = int(beg)
			end = int(end)
			score = 0 if score == '.' else float(score)
			if (beg,end) not in introns: introns[(beg,end)] = 0
			introns[(beg,end)] += score
			total += score

	# convert to histogram
	for k in introns: introns[k] /= total

	return introns

def expdiff(introns1, introns2):

	# non-mutating
	i1 = copy.deepcopy(introns1)
	i2 = copy.deepcopy(introns2)

	# ensure all introns are in both collections
	for k in i1:
		if k not in i2: i2[k] = 0
	for k in i2:
		if k not in i1: i1[k] = 0

	# compare distances
	p1 = []
	p2 = []
	details = []
	for k in i1:
		p1.append(i1[k])
		p2.append(i2[k])
		details.append((k, i1[k], i2[k]))

	distance = manhattan(p1, p2)
	return distance, details

