####################
## HMM Validation ##
####################

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