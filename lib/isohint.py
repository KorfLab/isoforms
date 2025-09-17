import subprocess
import os

import isoform

def countiso(dons, accs, min_intron, min_exon, limit=False):
    """count all possible isoform"""
    count   = 0

    def all_possible(dpos, apos, ipos, old):
        """recursion"""
        nonlocal count
        
        if limit:
            if count > limit: return
            
        test  = ipos == 0
        start = dpos      if test else apos
        end   = len(dons) if test else len(accs)
        
        for i in range(start, end):
            new = dons[i] if test else accs[i]
            if old != 0 and test:      
                if new - old - 1 < min_exon: continue
            else:                       
                if new - old + 1 < min_intron: continue
            if not test: count += 1
            if test:    all_possible(i + 1, apos, 1, dons[i])
            else:       all_possible(dpos, i + 1, 0, accs[i])
    
    # recursion
    all_possible(0, 0, 0, 0)
    return count

def hints(hmm, fasta, models_dir=None, method='gap', **kwargs):
    """Use HMM to cutoff whole gene for isoform prediction"""
    
    output     = run_hmm(hmm, fasta, models_dir)
    dons, accs = parse_hint(output)
    
    # methods
    cutoff_methods = {
        'gap':    gapstats,
        'smooth': lambda s: smoothed_gapstats(s, kwargs.get('k', 5)),
        'perc':   lambda s: percentile(s, kwargs.get('p', 25)),
        'topk':   lambda s: topk(s, kwargs.get('k', 10)),
        'thresh': lambda s: threshold(s, kwargs.get('t', 0.5)),
        'adapt':  lambda s: adaptive_cutoff(s, kwargs.get('min_sites', 5)),
    }
    
    if method not in cutoff_methods:
        raise ValueError(f"Unknown cutoff method: {method}")
    
    cutfn       = cutoff_methods[method]
    dons_cut    = cutfn(dons) if dons else []
    accs_cut    = cutfn(accs) if accs else []
    
    return dons_cut, accs_cut

############################################################
######################## Auxiliary #########################
############################################################

"""Pipeline Cmd"""

def run_hmm(hmm, fasta, models_dir=None):
    """Run HMM with hint settings"""

    if models_dir is None:
        exe_dir     = os.path.dirname(os.path.abspath(hmm))
        models_dir  = os.path.join(exe_dir, "..", "..", "models")
    
    models_dir = os.path.abspath(models_dir)
    
    if not os.path.exists(models_dir):
        raise ValueError(f"Models directory not found: {models_dir}")
    
    cmd = [
        hmm,
        "--sequence", fasta,
        "--don_emission", os.path.join(models_dir, "don.pwm"),
        "--acc_emission", os.path.join(models_dir, "acc.pwm"),
        "--exon_emission", os.path.join(models_dir, "exon.mm"),
        "--intron_emission", os.path.join(models_dir, "intron.mm"),
        "--ped_exon", os.path.join(models_dir, "exon.len"),
        "--ped_intron", os.path.join(models_dir, "intron.len"),
        "--print_splice"
    ]
    
    try:
        result = subprocess.run(cmd, check=True, text=True, capture_output=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error running HMM: {e}")
        print(f"Command: {' '.join(cmd)}")
        print(f"Stderr: {e.stderr}")
        raise

def run_geniso2(geniso, fasta, model, hints=False):
    """cmd to run geniso"""
    
    cmd   = [geniso, fasta, model]
    if hints:
        cmd.append('--hmm')
        cmd.append(str(hints))
    try:
        result = subprocess.run(cmd, check=True, text=True, capture_output=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Command failed with exit code {e.returncode}")
        print(f"STDERR: {e.stderr}")
        return f"Error in geniso2: {e.stderr}"

"""Parser"""

def parse_hmm(output, stovit=False):    
    dons = []
    accs = []
    switch = 0
    
    lines = output.strip().split('\n')
    for line in lines:
        if      line == "DONS": 
            switch = 0
            continue
        elif    line == "ACCS":
            switch = 1
            continue
        
        line = line.strip().split('\t')
        pos = int  (line[0])
        val = float(line[1])

        if   switch == 0: dons.append( (pos, val) )
        elif switch == 1: accs.append( (pos, val) )
    
    return dons, accs

def get_basis_sites(seq, hmm_dons, hmm_accs, flank=99, minex=25):
    """GT/AG sites that are NOT included in HMM output"""
    
    dons, accs  = isoform.gtag_sites(seq, flank, minex)    
    hmm_don_pos = set(pos for pos, score in hmm_dons) if hmm_dons else set()
    hmm_acc_pos = set(pos for pos, score in hmm_accs) if hmm_accs else set()
    basis_dons  = [pos for pos in dons if pos not in hmm_don_pos]
    basis_accs  = [pos for pos in accs if pos not in hmm_acc_pos]
    return basis_dons, basis_accs

def combine_with_basis(cutoff_sites, basis_sites):
    combined = set(cutoff_sites)
    combined.update(basis_sites)
    return sorted(list(combined))

def parse_rf(output):
    """Parse HMM for rf"""

    hints   = []
    typ     = None
    
    for line in output.strip().split('\n'):
        line = line.strip()
        
        if not line:    continue
            
        if line == "DONS":
            typ = "don"
            continue
        elif line == "ACCS":
            typ = "acc"
            continue
        
        if typ and '\t' in line or '     ' in line:
            parts = line.split()
            if len(parts) >= 2:
                pos     = int(parts[0])
                val     = float(parts[1])
                hints.append((pos, typ, val))
    
    return hints

"""Cutoff functions for hints"""

def gapstats(sites):
    if len(sites) < 2:
        return [s[0] for s in sites]
    
    sites   = sorted(sites, key=lambda x: x[1], reverse=True)
    max_gap = 0
    idx     = len(sites) - 1
    
    for i in range(len(sites) - 1):
        val  = sites[i][1]
        val2 = sites[i+1][1]
        gap  = val - val2  # log space
        
        if gap > max_gap:
            max_gap = gap 
            idx     = i
    
    return sorted([site[0] for site in sites[:idx+1]])

def smoothed_gapstats(sites, k=5):
    if len(sites) <= k:
        return [s[0] for s in sites]
    
    sites   = sorted(sites, key=lambda x: x[1], reverse=True)
    max_gap = 0
    idx     = len(sites) - 1
    
    for i in range(len(sites) - k):
        avg  = sum(s[1] for s in sites[i:i+k]) / k
        avg2 = sum(s[1] for s in sites[i+1:i+k+1]) / k
        gap  = avg - avg2
        
        if gap > max_gap:
            max_gap = gap
            idx     = i + k // 2
    
    return sorted([site[0] for site in sites[:idx+1]])

def percentile(sites, p=25):
    if not sites:
        return []
    
    sites = sorted(sites, key=lambda x: x[1], reverse=True)
    
    cutoff_idx = int(len(sites) * (100 - p) / 100)
    cutoff_idx = max(1, cutoff_idx)
    cutoff_val = sites[cutoff_idx - 1][1] if cutoff_idx <= len(sites) else sites[-1][1]
    
    filtered = [s for s in sites if s[1] >= cutoff_val]
    return sorted([s[0] for s in filtered])

def topk(sites, k=10):
    if not sites:
        return []
    
    sites = sorted(sites, key=lambda x: x[1], reverse=True)
    k = min(k, len(sites))
    
    return sorted([s[0] for s in sites[:k]])

def threshold(sites, t=-5.0):
    if not sites:
        return []
    
    filtered = [s for s in sites if s[1] >= t]
    return sorted([s[0] for s in filtered])

def adaptive_cutoff(sites, min_sites=5):
    if len(sites) <= min_sites:
        return [s[0] for s in sites]
    
    sites  = sorted(sites, key=lambda x: x[1], reverse=True)
    
    scores = [s[1] for s in sites]
    mean   = sum(scores) / len(scores)
    var    = sum((x - mean) ** 2 for x in scores) / len(scores)
    std    = var ** 0.5
    
    cutoff   = mean + std
    filtered = [s for s in sites if s[1] >= cutoff]
    
    if len(filtered) < min_sites:
        filtered = sites[:min_sites]
    
    return sorted([s[0] for s in filtered])

"""Debug Check"""

def mgtag_sites(seq, flank, minex, minin):
    dons = []
    accs = []
    
    for i in range(flank + minex, len(seq) - flank - minex):
        if seq[i:i+2] == 'GT':
            dons.append(i)
        if seq[i-1:i+1] == 'AG':
            accs.append(i)
    
    ndons = []
    naccs = []
    
    if dons and accs:
        first_donor   = dons[0]
        last_acceptor = accs[-1]

        naccs = [acc for acc in accs if acc >= first_donor + minin - 1]        
        ndons = [don for don in dons if don <= last_acceptor - minin + 1]
    
    return ndons, naccs