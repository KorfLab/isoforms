import random
import numpy as np

############################################################
##################### DevData Generator ####################
############################################################

''' simulator for output of the HMM model '''

def generate_dev_data(
        total_samples   = 20,
        don_ratio       = 0.6,
        value_range     = (0.0, 1.0),
        pos_range       = (0, 1000),
        seed            = None):

    if seed is not None:
        random.seed(seed)

    positions = random.sample(range(pos_range[0], pos_range[1] + 1), total_samples)
    random.shuffle(positions)

    num_dons = int(total_samples * don_ratio)
    dons     = positions[:num_dons]
    accs     = positions[num_dons:]

    pos2info = {
        p: (random.uniform(value_range[0], value_range[1]),
            'don' if p in dons else 'acc')
        for p in positions
    }

    return dons, accs, pos2info

############################################################
######################### Functions ########################
############################################################

def base2_to_int(bits) -> int:

    idx = 0
    for bit in bits:
        idx = (idx << 1) | bit
    return idx

def finding_outlier(doa, pos2info, z_threshold=None):
    '''
    doa     :   dons or accs
    pos2info:   dictionary for pos : val, typ
    '''
    assert(doa)

    mdoa = []

    for pos in doa:
        val, typ = pos2info[pos]
        mdoa.append((pos, typ, val))

    mdoa = sorted(mdoa, key=lambda x: x[2])
    vals  = np.asarray([val for _, _, val in mdoa], dtype=float)
    mean  = vals.mean()
    std   = vals.std(ddof=0)
    zscr  = (vals - mean)/std

    if not z_threshold: otest = np.abs(zscr) >= 3.0
    else:               otest = np.abs(zscr) >= z_threshold

    otl   = [info for info, tf in zip(mdoa, otest) if tf]
    doa   = [info for info, tf in zip(mdoa, otest) if not tf]

    return otl, doa

############################################################
##################### DecisionTree Class ###################
############################################################

'''
random forest parameter:
    mtry    : p/3
    nodesize: 30% of root node
    
    diff_split_coeff: current being inhibit
        if diff_coeff: expand the full binary tree
'''

class IsoformTree:

    def __init__(
        self,
        dons,
        accs,
        pos2info,
        min_sample_coeff    = 0.3,
        diff_split_coeff    = 0.9,
        gini_threshold      = None           # test each split for gini_impurity 
    ):
        
        self.dons               = dons
        self.accs               = accs
        self.pos2info           = pos2info
        self.min_samples_split  = int(len(self.dons+self.accs) * min_sample_coeff)
        self.gini_threshold     = gini_threshold
        self.diff_samples_split = self.min_samples_split * diff_split_coeff
        self.output             = dict()
        self.rules              = []

        self._btstrapping()

        dataset = [
            (pos,
             pos2info[pos][1],
             pos2info[pos][0])
             for pos in (self.dons + self.accs)
        ]

        self._recursion_tree(dataset, [])

    def _btstrapping(self):
        ''' Bootstrap original donor and acceptor sites '''
        self.dons = sorted(random.choices(self.dons, k=len(self.dons)))
        self.accs = sorted(random.choices(self.accs, k=len(self.accs)))

    def _fsubset(self, node, mtry):
        ''' Find random subset of current dataset '''
        return random.sample(node, mtry)

    def compute_mse(self, vals):
        ''' Compute the mean squared error '''
        vals = np.asarray(vals, dtype=float)

        return np.mean((vals - vals.mean())**2) if vals.size else 0.0

    def compute_gini(self, typs):
        ''' Compute the gini impurity '''
        count = {}
        total = len(typs)

        for typ in typs:
            count[typ] = count.get(typ, 0) + 1

        gini = 1.0
        for frq in count.values():
            p = frq / total
            gini -= p * p

        return gini

    def _gini_test(self, typs):
        ''' Compute the gini impurity test '''
        threshold   = self.gini_threshold if self.gini_threshold else 0.1
        gini        = self.compute_gini(typs)

        return gini < threshold

    def _find_threshold(self, right_split):
        ''' Find the split threshold for split original dataset '''
        right_split     = sorted(right_split, key=lambda x: x[2])
        _, _, threshold = right_split[0]
        return threshold

    def _split_node(self, threshold, node):
        ''' split the original node based on threshold '''
        node = sorted(node, key=lambda x: x[2])
        idx = 0
        for i, (*_, val) in enumerate(node):
            if val >= threshold:
                idx = i
                break
        
        left_node  = node[:idx]
        right_node = node[idx+1:]
        return left_node, right_node
        
    def _split(self, node):
        ''' 
            Compute the mean squared error gain for split criteria
            Additionally gini impurity is maintained for diversity of split quality
        '''
        subset      = self._fsubset(node, int(len(node)/3))
        parent_mse  = self.compute_mse([val for _, _, val in subset])
        subset      = sorted(subset, key=lambda x: x[2])

        all_psplit = {}

        for i in range(len(subset)-1):
            gain_mse    = parent_mse
            left_split  = subset[:i+1]
            right_split = subset[i+1:]
            left_mse    = self.compute_mse([val for _, _, val in left_split])
            right_mse   = self.compute_mse([val for _, _, val in right_split])
            gain_mse    = gain_mse - left_mse - right_mse
            threshold   = self._find_threshold(right_split)
            all_psplit.setdefault(gain_mse, []).append(threshold)

        max_gain = sorted(all_psplit.keys(), reverse=True)

        for mse in max_gain:
            thresholds = all_psplit[mse]

            if len(thresholds) == 1:
                ''' if single threshold '''
                threshold = thresholds[0]
                left_node, right_node = self._split_node(threshold, node)

                ltest = len(left_node)  < self.min_samples_split
                rtest = len(right_node) < self.min_samples_split
                if ltest and rtest: continue

                '''
                diff_test = abs(len(left_node)-len(right_node))
                if diff_test > self.diff_samples_split: continue
                '''
                
                ltest = self._gini_test([typ for _, typ, _ in left_node])
                rtest = self._gini_test([typ for _, typ, _ in right_node])
                if ltest or rtest:  continue

                self.rules.append(threshold)
                return

            else:
                ''' if mse collision happen '''
                split_ginis = []

                for threshold in thresholds:
                    left_node, right_node = self._split_node(threshold, node)
                    
                    ltest = len(left_node)  < self.min_samples_split
                    rtest = len(right_node) < self.min_samples_split
                    if ltest and rtest: continue

                    '''
                    diff_test = abs(len(left_split)-len(right_split))
                    if diff_test > self.diff_samples_split: continue
                    '''
                    
                    left_gini  = self.compute_gini([typ for _, typ, _ in left_node])
                    right_gini = self.compute_gini([typ for _, typ, _ in right_node])
                    split_ginis.append(left_gini)
                    split_ginis.append(right_gini)

                if not split_ginis or max(split_ginis) <= 0.1:  return False
                max_idx = split_ginis.index(max(split_ginis))
                
                if max_idx % 2 == 0:
                    threshold = thresholds[max_idx%2]
                    self.rules.append(threshold)
                    return
                else:
                    threshold = thresholds[(max_idx-1)%2]
                    self.rules.append(threshold)
                    return
                    
        return False

    def _store_output(self, node, path):
        ''' Store the leaf node of decision tree into base2 keys '''
        id = base2_to_int(path)
        self.output[id] = node

    def _recursion_tree(self, node, path):
        ''' Recursion until death '''

        if not node:
            self._store_output(node, path)
            return
        
        if len(node) < self.min_samples_split:
            self._store_output(node, path)
            return
        
        # split on subset
        split = self._split(node)

        if split is False:
            self._store_output(node, path)
            return
        
        # split original dataset
        split_threshold = self.rules[-1]
        left_node, right_node = self._split_node(split_threshold, node)

        self._recursion_tree(left_node,  [0] + path)
        self._recursion_tree(right_node, [1] + path)

############################################################
##################### Classifier Class #####################
############################################################

class IsoformClassifer:
    
    def __init__(
        self,
        data,              # (pos, typ, val) : (int, dons or accs, float)
        rules,
        output,
    ):
        self.data  = data
        self.rules  = rules
        self.output = output
        
        self.feature    = [val for _, _, val in self.data][0]
        self.prediction = []
        self._rules2base2(self.rules)
        
        idx        = base2_to_int(self.prediction)
        
        dons = [pos for typ, pos, _ in [self.output[idx] + self.input] if typ == 'dons']
        accs = [pos for typ, pos, _ in [self.output[idx] + self.input] if typ == 'accs']
        
        return dons, accs
    
    def _rules2base2(self, node):
        
        if not node:
            return
        
        root_split = node[0]
        
        if self.feature < root_split:
            next_node = [val for val in node if val < root_split]
            self.prediction.append(0)
        else:
            next_node = [val for val in node if val > root_split]
            self.prediction.append(1)
        
        self._rules2base2(next_node)