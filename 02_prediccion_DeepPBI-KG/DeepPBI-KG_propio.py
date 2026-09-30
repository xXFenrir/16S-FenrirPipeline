from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from itertools import cycle
import itertools
from Bio import SeqIO
from collections import OrderedDict
import numpy as np
import pandas as pd
import re
import os
import torch
#torch.cuda.current_device() 
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, classification_report, confusion_matrix
from sklearn.preprocessing import label_binarize
from sklearn.metrics import roc_curve, auc, mean_squared_error, precision_recall_curve
import argparse
import pickle
import shutil
import time
from concurrent.futures import ProcessPoolExecutor
N_WORKERS = max(1, (os.cpu_count() or 1) - 1)
start_time = time.time()
import warnings
warnings.filterwarnings("ignore")

#define model
class DNN(nn.Module):
    def __init__(self):
        super(DNN,self).__init__()
        self.encoder  =  nn.Sequential(
            nn.Linear(2676,1024),
            nn.ReLU(),
            nn.Linear(1024, 512),
            nn.ReLU(),
            nn.Linear(512, 256),
            nn.ReLU(),
            nn.Linear(256, 128),
            nn.ReLU(),
            nn.Linear(128, 64),
            nn.ReLU(),
            nn.Linear(64, 1),
            nn.Sigmoid()
        )

    def forward(self, x):
        x = self.encoder(x)
        return x

#Discard data with empty features
def check_blank_data(data):
    val = data.sum(axis = 1).tolist()
    ind = []
    for i in range(len(val)):
        if val[i] == 0:
            pass
        else:
            ind.append(i)
    data = data.iloc[ind,:]
    
    return data
    

def predict(phage_fea, host_fea, scaler, parameter):
    # read data
    phage = pd.read_csv(phage_fea, index_col = 0, header = 0)
    host= pd.read_csv(host_fea, index_col = 0, header = 0)

    phage = check_blank_data(phage)
    host = check_blank_data(host)

    print("loading data and predicting on CPU by batches...")
    phage_name = phage.index.tolist()
    host_name = host.index.tolist()
    n_phage = len(phage)
    n_host = len(host)

    # Cargar el scaler y el modelo en la CPU
    scaler = pickle.load(open(scaler, 'rb'))
    model = DNN()
    model.load_state_dict(torch.load((parameter), map_location = torch.device("cpu")))
    model.eval()
    torch.set_num_threads(N_WORKERS)

    phage_vals = phage.values
    host_vals = host.values

    # Group several phages per batch instead of one phage (= n_host rows) per
    # iteration: with n_phage in the thousands, per-iteration Python/torch call
    # overhead dominated runtime far more than the actual matrix math.
    # Preallocating a float array (instead of writing into an object-dtype
    # DataFrame row by row) also avoids pandas' slow per-row block updates.
    rows_target = 10000
    batch_size = max(1, rows_target // max(n_host, 1))
    out_arr = np.empty((n_phage, n_host), dtype = np.float32)

    with torch.no_grad():
        for start in range(0, n_phage, batch_size):
            end = min(start + batch_size, n_phage)
            k = end - start

            phage_feats = np.repeat(phage_vals[start:end], n_host, axis = 0)
            host_feats = np.tile(host_vals, (k, 1))

            fea = np.concatenate((phage_feats, host_feats), axis = 1)
            val_data = scaler.transform(fea)
            val_tensor = torch.from_numpy(val_data).float()

            result = model(val_tensor)
            out_arr[start:end, :] = result.detach().numpy().reshape(k, n_host)

    out = pd.DataFrame(out_arr, index = phage_name, columns = host_name)
    return out
    

def protein_features(protein_sequences):
    """
    This function calculates a number of basic properties for a list of protein sequences
    
    Input: list of protein sequences (as strings), length can also be 1
    Output: a dataframe of features
    """
    
    import numpy as np
    import pandas as pd
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    from collections import Counter
    
    # AA frequency and protein characteristics
    mol_weight = []; aromaticity = []; instability = []; flexibility = []; prot_length = []
    pI = []; helix_frac = []; turn_frac = []; sheet_frac = []
    frac_aliph = []; frac_unch_polar = []; frac_polar = []; frac_hydrophob = []; frac_pos = []; frac_sulfur = []
    frac_neg = []; frac_amide = []; frac_alcohol = []; C2 = []; H2 = []; O2 = []; N2 = []; S2 = []
    AA_dict = {'G': [], 'A': [], 'V': [], 'L': [], 'I': [], 'F': [], 'P': [], 'S': [], 'T': [], 'Y': [],
           'Q': [], 'N': [], 'E': [], 'D': [], 'W': [], 'H': [], 'R': [], 'K': [], 'M': [], 'C': []}
    
    # calculate physical_chemical_feature
    for items in protein_sequences:
        items=items.replace('X','').replace('U','').replace('B','').replace('Z','').replace('*','')
        CE = 'CHONS'
        Chemi_stats = {'A':{'C': 3, 'H': 7, 'O': 2, 'N': 1, 'S': 0},
                       'C':{'C': 3, 'H': 7, 'O': 2, 'N': 1, 'S': 1},
                       'D':{'C': 4, 'H': 7, 'O': 4, 'N': 1, 'S': 0},
                       'E':{'C': 5, 'H': 9, 'O': 4, 'N': 1, 'S': 0},
                       'F':{'C': 9, 'H': 11,'O': 2, 'N': 1, 'S': 0},
                       'G':{'C': 2, 'H': 5, 'O': 2, 'N': 1, 'S': 0},
                       'H':{'C': 6, 'H': 9, 'O': 2, 'N': 3, 'S': 0},
                       'I':{'C': 6, 'H': 13,'O': 2, 'N': 1, 'S': 0},
                       'K':{'C': 6, 'H': 14,'O': 2, 'N': 2, 'S': 0},
                       'L':{'C': 6, 'H': 13,'O': 2, 'N': 1, 'S': 0},
                       'M':{'C': 5, 'H': 11,'O': 2, 'N': 1, 'S': 1},
                       'N':{'C': 4, 'H': 8, 'O': 3, 'N': 2, 'S': 0},
                       'P':{'C': 5, 'H': 9, 'O': 2, 'N': 1, 'S': 0},
                       'Q':{'C': 5, 'H': 10,'O': 3, 'N': 2, 'S': 0},
                       'R':{'C': 6, 'H': 14,'O': 2, 'N': 4, 'S': 0},
                       'S':{'C': 3, 'H': 7, 'O': 3, 'N': 1, 'S': 0},
                       'T':{'C': 4, 'H': 9, 'O': 3, 'N': 1, 'S': 0},
                       'V':{'C': 5, 'H': 11,'O': 2, 'N': 1, 'S': 0},
                       'W':{'C': 11,'H': 12,'O': 2, 'N': 2, 'S': 0},
                       'Y':{'C': 9, 'H': 11,'O': 3, 'N': 1, 'S': 0}
                    }
    
        count = Counter(items)
        code = []
    
        for c in CE:
            abundance_c = 0
            for key in count:
                num_c = Chemi_stats[key][c]
                abundance_c += num_c * count[key]
                
            code.append(abundance_c)
            
        C, H, O, N, S = code
        C2.append(C)
        H2.append(H)
        O2.append(O)
        N2.append(N)
        S2.append(S)

    
    for item in protein_sequences:
        item=item.replace('X','').replace('U','').replace('B','').replace('Z','').replace('*','')
        # calculate various protein properties
        prot_length.append(len(item))
        frac_aliph.append((item.count('A')+item.count('G')+item.count('I')+item.count('L')+item.count('P')
                       +item.count('V'))/len(item))
        frac_unch_polar.append((item.count('S')+item.count('T')+item.count('N')+item.count('Q'))/len(item))
        frac_polar.append((item.count('Q')+item.count('N')+item.count('H')+item.count('S')+item.count('T')+item.count('Y')
                      +item.count('C')+item.count('M')+item.count('W'))/len(item))
        frac_hydrophob.append((item.count('A')+item.count('G')+item.count('I')+item.count('L')+item.count('P')
                        +item.count('V')+item.count('F'))/len(item))
        frac_pos.append((item.count('H')+item.count('K')+item.count('R'))/len(item))
        frac_sulfur.append((item.count('C')+item.count('M'))/len(item))
        frac_neg.append((item.count('D')+item.count('E'))/len(item))
        frac_amide.append((item.count('N')+item.count('Q'))/len(item))
        frac_alcohol.append((item.count('S')+item.count('T'))/len(item))
        protein_chars = ProteinAnalysis(item) 
        mol_weight.append(protein_chars.molecular_weight())
        aromaticity.append(protein_chars.aromaticity())
        instability.append(protein_chars.instability_index())
        flexibility.append(np.mean(protein_chars.flexibility()))
        pI.append(protein_chars.isoelectric_point())
        H, T, S = protein_chars.secondary_structure_fraction()
        helix_frac.append(H)
        turn_frac.append(T)
        sheet_frac.append(S)
    
        # calculate AA frequency
        for key in AA_dict.keys():
            AA_dict[key].append(item.count(key)/len(item))
            
    # make new dataframe & return
    features_protein = pd.DataFrame.from_dict(AA_dict)
    features_protein['protein_length'] = np.asarray(prot_length)
    features_protein['mol_weight'] = np.asarray(mol_weight)
    features_protein['aromaticity'] = np.asarray(aromaticity)
    features_protein['instability'] = np.asarray(instability)
    features_protein['flexibility'] = np.asarray(flexibility)
    features_protein['pI'] = np.asarray(pI)
    features_protein['frac_aliphatic'] = np.asarray(frac_aliph)
    features_protein['frac_uncharged_polar'] = np.asarray(frac_unch_polar)
    features_protein['frac_polar'] = np.asarray(frac_polar)
    features_protein['frac_hydrophobic'] = np.asarray(frac_hydrophob)
    features_protein['frac_positive'] = np.asarray(frac_pos)
    features_protein['frac_sulfur'] = np.asarray(frac_sulfur)
    features_protein['frac_negative'] = np.asarray(frac_neg)
    features_protein['frac_amide'] = np.asarray(frac_amide)
    features_protein['frac_alcohol'] = np.asarray(frac_alcohol)
    features_protein['AA_frac_helix'] = np.asarray(helix_frac)
    features_protein['AA_frac_turn'] = np.asarray(turn_frac)
    features_protein['AA_frac_sheet'] = np.asarray(sheet_frac)
    features_protein['C2'] = np.asarray(C2)
    features_protein['H2'] = np.asarray(H2)
    features_protein['O2'] = np.asarray(O2)
    features_protein['N2'] = np.asarray(N2)
    features_protein['S2'] = np.asarray(S2)
    
    return features_protein


# PROTEIN FEATURE: COMPOSITION
# --------------------------------------------------
def Count1(seq1, seq2):
	sum = 0
	for aa in seq1:
		sum = sum + seq2.count(aa)
	return sum

# C/T/D model
def CTDC(sequence):
	group1 = {
		'hydrophobicity_PRAM900101': 'RKEDQN',
		'hydrophobicity_ARGP820101': 'QSTNGDE',
		'hydrophobicity_ZIMJ680101': 'QNGSWTDERA',
		'hydrophobicity_PONP930101': 'KPDESNQT',
		'hydrophobicity_CASG920101': 'KDEQPSRNTG',
		'hydrophobicity_ENGD860101': 'RDKENQHYP',
		'hydrophobicity_FASG890101': 'KERSQD',
		'normwaalsvolume': 'GASTPDC',
		'polarity':        'LIFWCMVY',
		'polarizability':  'GASDT',
		'charge':          'KR',
		'secondarystruct': 'EALMQKRH',
		'solventaccess':   'ALFCGIVW'
	}
	group2 = {
		'hydrophobicity_PRAM900101': 'GASTPHY',
		'hydrophobicity_ARGP820101': 'RAHCKMV',
		'hydrophobicity_ZIMJ680101': 'HMCKV',
		'hydrophobicity_PONP930101': 'GRHA',
		'hydrophobicity_CASG920101': 'AHYMLV',
		'hydrophobicity_ENGD860101': 'SGTAW',
		'hydrophobicity_FASG890101': 'NTPG',
		'normwaalsvolume': 'NVEQIL',
		'polarity':        'PATGS',
		'polarizability':  'CPNVEQIL',
		'charge':          'ANCQGHILMFPSTWYV',
		'secondarystruct': 'VIYCWFT',
		'solventaccess':   'RKQEND'
	}

	property = ['hydrophobicity_PRAM900101', 'hydrophobicity_ARGP820101', 'hydrophobicity_ZIMJ680101', 'hydrophobicity_PONP930101',
	'hydrophobicity_CASG920101', 'hydrophobicity_ENGD860101', 'hydrophobicity_FASG890101', 'normwaalsvolume',
	'polarity', 'polarizability', 'charge', 'secondarystruct', 'solventaccess']

	for p in property:
		c1 = Count1(group1[p], sequence) / len(sequence)
		c2 = Count1(group2[p], sequence) / len(sequence)
		c3 = 1 - c1 - c2
		encoding = [c1, c2, c3]
        
	return encoding


# PROTEIN FEATURE: TRANSITION
# --------------------------------------------------
def CTDT(sequence):
    """
    Every number in the encoding tells us as a percentage over all AA pairs, how many transitions
    occured from group x to y or vice versa for a specific property p.
    """
    group1 = {'hydrophobicity_PRAM900101': 'RKEDQN','hydrophobicity_ARGP820101': 'QSTNGDE',
              'hydrophobicity_ZIMJ680101': 'QNGSWTDERA','hydrophobicity_PONP930101': 'KPDESNQT',
              'hydrophobicity_CASG920101': 'KDEQPSRNTG','hydrophobicity_ENGD860101': 'RDKENQHYP',
              'hydrophobicity_FASG890101': 'KERSQD','normwaalsvolume': 'GASTPDC','polarity': 'LIFWCMVY',
              'polarizability': 'GASDT','charge':'KR', 'secondarystruct': 'EALMQKRH','solventaccess': 'ALFCGIVW'}
    
    group2 = {'hydrophobicity_PRAM900101': 'GASTPHY','hydrophobicity_ARGP820101': 'RAHCKMV',
           'hydrophobicity_ZIMJ680101': 'HMCKV','hydrophobicity_PONP930101': 'GRHA',
           'hydrophobicity_CASG920101': 'AHYMLV','hydrophobicity_ENGD860101': 'SGTAW',
           'hydrophobicity_FASG890101': 'NTPG','normwaalsvolume': 'NVEQIL','polarity': 'PATGS',
           'polarizability': 'CPNVEQIL','charge': 'ANCQGHILMFPSTWYV', 
           'secondarystruct': 'VIYCWFT', 'solventaccess': 'RKQEND'}
    
    group3 = {'hydrophobicity_PRAM900101': 'CLVIMFW','hydrophobicity_ARGP820101': 'LYPFIW',
           'hydrophobicity_ZIMJ680101': 'LPFYI','hydrophobicity_PONP930101': 'YMFWLCVI',
           'hydrophobicity_CASG920101': 'FIWC','hydrophobicity_ENGD860101': 'CVLIMF',
           'hydrophobicity_FASG890101': 'AYHWVMFLIC','normwaalsvolume': 'MHKFRYW',
           'polarity': 'HQRKNED','polarizability': 'KMHFRYW','charge': 'DE',
           'secondarystruct': 'GNPSD','solventaccess': 'MSPTHY'}
    
    property = ['hydrophobicity_PRAM900101', 'hydrophobicity_ARGP820101', 'hydrophobicity_ZIMJ680101', 'hydrophobicity_PONP930101',
	'hydrophobicity_CASG920101', 'hydrophobicity_ENGD860101', 'hydrophobicity_FASG890101', 'normwaalsvolume',
	'polarity', 'polarizability', 'charge', 'secondarystruct', 'solventaccess']
    
    encoding = []
    aaPair = [sequence[j:j+2] for j in range(len(sequence)-1)]
    
    for p in property:
        c1221, c1331, c2332 = 0, 0, 0
        for pair in aaPair:
            if (pair[0] in group1[p] and pair[1] in group2[p]) or (pair[0] in group2[p] and pair[1] in group1[p]):
                c1221 += 1
                continue
            if (pair[0] in group1[p] and pair[1] in group3[p]) or (pair[0] in group3[p] and pair[1] in group1[p]):
                c1331 += 1
                continue
            if (pair[0] in group2[p] and pair[1] in group3[p]) or (pair[0] in group3[p] and pair[1] in group2[p]):
                c2332 += 1
        encoding.append(c1221/len(aaPair))
        encoding.append(c1331/len(aaPair))
        encoding.append(c2332/len(aaPair))
    
    return encoding


# PROTEIN FEATURE: Z-SCALE
# --------------------------------------------------
def zscale(sequence):
    zdict = {
		'A': [0.24,  -2.32,  0.60, -0.14,  1.30], # A
		'C': [0.84,  -1.67,  3.71,  0.18, -2.65], # C
		'D': [3.98,   0.93,  1.93, -2.46,  0.75], # D
		'E': [3.11,   0.26, -0.11, -0.34, -0.25], # E
		'F': [-4.22,  1.94,  1.06,  0.54, -0.62], # F
		'G': [2.05,  -4.06,  0.36, -0.82, -0.38], # G
		'H': [2.47,   1.95,  0.26,  3.90,  0.09], # H
		'I': [-3.89, -1.73, -1.71, -0.84,  0.26], # I
		'K': [2.29,   0.89, -2.49,  1.49,  0.31], # K
		'L': [-4.28, -1.30, -1.49, -0.72,  0.84], # L
		'M': [-2.85, -0.22,  0.47,  1.94, -0.98], # M
		'N': [3.05,   1.62,  1.04, -1.15,  1.61], # N
		'P': [-1.66,  0.27,  1.84,  0.70,  2.00], # P
		'Q': [1.75,   0.50, -1.44, -1.34,  0.66], # Q
		'R': [3.52,   2.50, -3.50,  1.99, -0.17], # R
		'S': [2.39,  -1.07,  1.15, -1.39,  0.67], # S
		'T': [0.75,  -2.18, -1.12, -1.46, -0.40], # T
		'V': [-2.59, -2.64, -1.54, -0.85, -0.02], # V
		'W': [-4.36,  3.94,  0.59,  3.44, -1.59], # W
		'Y': [-2.54,  2.44,  0.43,  0.04, -1.47], # Y
		'-': [0.00,   0.00,  0.00,  0.00,  0.00], # -
	}

    sequence=sequence.replace('X','-').replace('U','-').replace('B','-').replace('Z','-').replace('*','-')
    z1, z2, z3, z4, z5 = 0, 0, 0, 0, 0
    for aa in sequence:
        z1 += zdict[aa][0]
        z2 += zdict[aa][1]
        z3 += zdict[aa][2]
        z4 += zdict[aa][3]
        z5 += zdict[aa][4]
    encoding = [z1/len(sequence), z2/len(sequence), z3/len(sequence), z4/len(sequence), z5/len(sequence)]
    
    return encoding

# DNA feature
def dna_features(dna_sequences):
    """
    This function calculates a variety of properties from a DNA sequence.
    
    Input: a list of DNA sequence (can also be length of 1)
    Output: a dataframe of features
    """
    
    import numpy as np
    import pandas as pd
    from Bio.SeqUtils import GC, CodonUsage

    A_freq = []; T_freq = []; C_freq = []; G_freq = []; GC_content = []
    codontable = {'ATA':[], 'ATC':[], 'ATT':[], 'ATG':[], 'ACA':[], 'ACC':[], 'ACG':[], 'ACT':[],
    'AAC':[], 'AAT':[], 'AAA':[], 'AAG':[], 'AGC':[], 'AGT':[], 'AGA':[], 'AGG':[],
    'CTA':[], 'CTC':[], 'CTG':[], 'CTT':[], 'CCA':[], 'CCC':[], 'CCG':[], 'CCT':[],
    'CAC':[], 'CAT':[], 'CAA':[], 'CAG':[], 'CGA':[], 'CGC':[], 'CGG':[], 'CGT':[],
    'GTA':[], 'GTC':[], 'GTG':[], 'GTT':[], 'GCA':[], 'GCC':[], 'GCG':[], 'GCT':[],
    'GAC':[], 'GAT':[], 'GAA':[], 'GAG':[], 'GGA':[], 'GGC':[], 'GGG':[], 'GGT':[],
    'TCA':[], 'TCC':[], 'TCG':[], 'TCT':[], 'TTC':[], 'TTT':[], 'TTA':[], 'TTG':[],
    'TAC':[], 'TAT':[], 'TAA':[], 'TAG':[], 'TGC':[], 'TGT':[], 'TGA':[], 'TGG':[]}
    
    for item in dna_sequences:
        # nucleotide frequencies
        A_freq.append(item.count('A')/len(item))
        T_freq.append(item.count('T')/len(item))
        C_freq.append(item.count('C')/len(item))
        G_freq.append(item.count('G')/len(item))
    
        # GC content
        GC_content.append(GC(item))
    
        # codon frequency: count codons, normalize counts, add to dict
        codons = [item[i:i+3] for i in range(0, len(item), 3)]
        l = []
        for key in codontable.keys():
            l.append(codons.count(key))
        l_norm = [float(i)/sum(l) for i in l]
        
        for j, key in enumerate(codontable.keys()):
            codontable[key].append(l_norm[j])
     
    # codon usage bias (_b)
    synonym_codons = CodonUsage.SynonymousCodons
    codontable2 = {'ATA_b':[], 'ATC_b':[], 'ATT_b':[], 'ATG_b':[], 'ACA_b':[], 'ACC_b':[], 'ACG_b':[], 'ACT_b':[],
    'AAC_b':[], 'AAT_b':[], 'AAA_b':[], 'AAG_b':[], 'AGC_b':[], 'AGT_b':[], 'AGA_b':[], 'AGG_b':[],
    'CTA_b':[], 'CTC_b':[], 'CTG_b':[], 'CTT_b':[], 'CCA_b':[], 'CCC_b':[], 'CCG_b':[], 'CCT_b':[],
    'CAC_b':[], 'CAT_b':[], 'CAA_b':[], 'CAG_b':[], 'CGA_b':[], 'CGC_b':[], 'CGG_b':[], 'CGT_b':[],
    'GTA_b':[], 'GTC_b':[], 'GTG_b':[], 'GTT_b':[], 'GCA_b':[], 'GCC_b':[], 'GCG_b':[], 'GCT_b':[],
    'GAC_b':[], 'GAT_b':[], 'GAA_b':[], 'GAG_b':[], 'GGA_b':[], 'GGC_b':[], 'GGG_b':[], 'GGT_b':[],
    'TCA_b':[], 'TCC_b':[], 'TCG_b':[], 'TCT_b':[], 'TTC_b':[], 'TTT_b':[], 'TTA_b':[], 'TTG_b':[],
    'TAC_b':[], 'TAT_b':[], 'TAA_b':[], 'TAG_b':[], 'TGC_b':[], 'TGT_b':[], 'TGA_b':[], 'TGG_b':[]}

    for item1 in dna_sequences:
        codons = [item1[l:l+3] for l in range(0, len(item1), 3)]
        codon_counts = []
    
        # count codons corresponding to codontable (not codontable2 because keynames changed!)
        for key in codontable.keys():
            codon_counts.append(codons.count(key))
        
        # count total for synonymous codons, divide each synonym codon count by total
        for key_syn in synonym_codons.keys():
            total = 0
            for item2 in synonym_codons[key_syn]:
                total += codons.count(item2)
            for j, key_table in enumerate(codontable.keys()):
                if (key_table in synonym_codons[key_syn]) & (total != 0):
                    codon_counts[j] /= total
                
        # add corrected counts to codontable2 (also corresponds to codontable which was used to count codons)
        for k, key_table in enumerate(codontable2.keys()):
            codontable2[key_table].append(codon_counts[k])
            
    # make new dataframes & standardize
    features_codonbias = pd.DataFrame.from_dict(codontable2)
    features_dna = pd.DataFrame.from_dict(codontable)
    features_dna['A_freq'] = np.asarray(A_freq)
    features_dna['T_freq'] = np.asarray(T_freq)
    features_dna['C_freq'] = np.asarray(C_freq)
    features_dna['G_freq'] = np.asarray(G_freq)
    features_dna['GC'] = np.asarray(GC_content)
    
    # concatenate dataframes & return
    features = pd.concat([features_dna, features_codonbias], axis=1)
    return features



def _dna_process_one(args):
    """Compute the 798-length dna feature vector for a single cds file. Runs in a worker process."""
    ls, path, ref, pattern = args
    res_dir = path + os.sep + ls
    records = [r for r in SeqIO.parse(res_dir, "fasta")]

    ID = []
    seq = []
    kg_num = 0
    for i in records:
        if pattern == 'key_gene':
            tags = re.findall(r'\[(.+?)\]', i.description)
            if tags[0][:4] == 'gene' and tags[0] in ref:
                ID.append(i.id)
                seq.append(''.join(list(i.seq)))
                kg_num += 1
            elif tags[1] != 'protein=hypothetical protein' and tags[1] in ref:
                ID.append(i.id)
                seq.append(''.join(list(i.seq)))
                kg_num += 1
        else:
            ID.append(i.id)
            seq.append(''.join(list(i.seq)))

    features = dna_features(seq)
    features.index = ID

    key = ls[4:-6]
    if features.empty:
        return (key, None, kg_num if pattern == 'key_gene' else None)

    KP = features
    KP_diff = pd.DataFrame(index = ['mean', 'max', 'min', 'std', 'median', 'var'], columns = KP.columns[1:])
    KP_diff = KP_diff.replace(np.nan, 0)
    for i in KP.columns:
        KP_diff.loc[:, i] = [np.mean(KP[i]), max(KP[i]), min(KP[i]), np.std(KP[i]), np.median(KP[i]), np.var(KP[i])]
    KP_diff = pd.concat([KP_diff.iloc[0,:], KP_diff.iloc[1,:], KP_diff.iloc[2,:], KP_diff.iloc[3,:], KP_diff.iloc[4,:], KP_diff.iloc[5,:]])

    return (key, np.asarray(KP_diff), kg_num if pattern == 'key_gene' else None)


def dna_process(path, key_gene, pattern):

    walkfile = walkFile(path)
    if pattern == 'key_gene':
        ref = pd.read_csv(key_gene, header = 0)
        ref = ref.iloc[:,0].tolist()
    else:
        ref = None
    # create feature table
    CDD_diff = list(range(798))
    ID_diff = []
    for ls in walkfile:
        ID_diff.append(ls[4:-6])

    #feature process (parallel across cds files)
    key_gene_num = {}
    results = {}
    cds_files = [ls for ls in walkfile if ls[:3] == 'cds']
    tasks = [(ls, path, ref, pattern) for ls in cds_files]

    with ProcessPoolExecutor(max_workers = N_WORKERS) as executor:
        for key, values, kg_num in executor.map(_dna_process_one, tasks):
            if pattern == 'key_gene':
                key_gene_num[key] = kg_num
            if values is not None:
                results[key] = values

    # build the table in one vectorized step instead of assigning row-by-row (which is O(n^2) in pandas)
    KP_df = pd.DataFrame.from_dict(results, orient = 'index', columns = CDD_diff)
    KP_df = KP_df.reindex(ID_diff)
    KP_df = KP_df.fillna(0)

    return KP_df, key_gene_num


def _protein_process_one(args):
    """Compute the 540-length protein feature vector for a single cds file. Runs in a worker process."""
    ls, path, ref, pattern = args
    res_dir = path + os.sep + ls
    records = [r for r in SeqIO.parse(res_dir, "fasta")]

    ID = []
    seq = []
    for i in records:
        if pattern == 'key_gene':
            tags = re.findall(r'\[(.+?)\]', i.description)
            if tags[0][:4] == 'gene' and tags[0] in ref:
                ID.append(i.id)
                seq.append(''.join(list(i.seq)))
            elif tags[1] != 'protein=hypothetical protein' and tags[1] in ref:
                ID.append(i.id)
                seq.append(''.join(list(i.seq)))
        else:
            ID.append(i.id)
            seq.append(''.join(list(i.seq)))

    # protein features: CTD & Z-scale & protein_features
    extra_feats = np.zeros((len(seq), 47))
    for i, item in enumerate(seq):
        feature_lst = []
        feature_lst += CTDC(item)
        feature_lst += CTDT(item)
        feature_lst += zscale(item)
        extra_feats[i,:] = feature_lst

    extra_feats_df = pd.DataFrame(extra_feats, columns=['CTDC1', 'CTDC2', 'CTDC3', 'CTDT1', 'CTDT2', 'CTDT3',
                             'CTDT4', 'CTDT5', 'CTDT6', 'CTDT7', 'CTDT8', 'CTDT9', 'CTDT10', 'CTDT11', 'CTDT12', 'CTDT13',
                             'CTDT14', 'CTDT15', 'CTDT16', 'CTDT17', 'CTDT18', 'CTDT19', 'CTDT20', 'CTDT21', 'CTDT22',
                             'CTDT23', 'CTDT24', 'CTDT25', 'CTDT26', 'CTDT27', 'CTDT28', 'CTDT29', 'CTDT30', 'CTDT31',
                             'CTDT32', 'CTDT33', 'CTDT34', 'CTDT35', 'CTDT36', 'CTDT37', 'CTDT38', 'CTDT39', 'Z1', 'Z2',
                             'Z3', 'Z4', 'Z5'])

    protein_feats = protein_features(seq)
    features = pd.concat([protein_feats, extra_feats_df], axis=1)
    features.index = ID

    key = ls[8:-6]
    if features.empty:
        return (key, None)

    KP = features
    KP_diff = pd.DataFrame(index = ['mean', 'max', 'min', 'std', 'median', 'var'], columns = KP.columns[1:])
    KP_diff = KP_diff.replace(np.nan, 0)
    for i in KP.columns:
        KP_diff.loc[:, i] = [np.mean(KP[i]), max(KP[i]), min(KP[i]), np.std(KP[i]), np.median(KP[i]), np.var(KP[i])]
    KP_diff = pd.concat([KP_diff.iloc[0,:], KP_diff.iloc[1,:], KP_diff.iloc[2,:], KP_diff.iloc[3,:], KP_diff.iloc[4,:], KP_diff.iloc[5,:]])

    return (key, np.asarray(KP_diff))


def protein_process(path, key_gene, pattern):

    walkfile = walkFile(path)
    if pattern == 'key_gene':
        ref = pd.read_csv(key_gene, header = 0)
        ref = ref.iloc[:,0].tolist()
    else:
        ref = None
    # create feature table
    CDD_diff = list(range(540))
    ID_diff = []
    for ls in walkfile:
        ID_diff.append(ls[8:-6])

    #feature process (parallel across cds files)
    results = {}
    cds_files = [ls for ls in walkfile if ls[:3] == 'cds']
    tasks = [(ls, path, ref, pattern) for ls in cds_files]

    with ProcessPoolExecutor(max_workers = N_WORKERS) as executor:
        for key, values in executor.map(_protein_process_one, tasks):
            if values is not None:
                results[key] = values

    # build the table in one vectorized step instead of assigning row-by-row (which is O(n^2) in pandas)
    KP_df = pd.DataFrame.from_dict(results, orient = 'index', columns = CDD_diff)
    KP_df = KP_df.reindex(ID_diff)
    KP_df = KP_df.fillna(0)

    return KP_df

def fg(dna_file, protein_file, Gene, Pattern, species, res_path):
    phage_dna, kg_num = dna_process(dna_file, Gene, Pattern)
    phage_protein = protein_process(protein_file, Gene, Pattern)
    phage = pd.concat([phage_dna, phage_protein], axis = 1)
    index = []
    for i in range(6):
        index.extend([132 + i * 133] + list(range(0 + i * 133, 132 + i * 133)))
    for i in range(6):
        index.extend([89 + i * 90 + 798] + list(range(0 + i * 90 + 798, 89 + i * 90 + 798)))
    phage.iloc[:,index].to_csv(res_path + os.sep + species + '_' + Pattern + '_feature.csv', index = True)
    return kg_num

    
def format_fasta(ana, seq, num):
    """
    Format the text in fasta format
    :param ana: annotation information
    :param seq: sequence
    :param num: The number of characters when the sequence is wrapped
    :return: fasta format text
    """
    format_seq = ""
    for i, char in enumerate(seq):
        format_seq += char
        if (i + 1) % num == 0:
            format_seq += "\n"
    return ana + format_seq + "\n"


def get_complete(gb_file):
    """
    Extract the cds sequence and its complete sequence from the genbank file
    :param gb_file: genbank file path
    :param f_cds: Whether to obtain only one CDS sequence
    :return: fasta CDS sequence in fasta format, complete sequence in fasta format
    """
    # Extract the complete sequence and format it in fasta
    gb_seq = SeqIO.parse(open(gb_file), "genbank")
    complete_total = []
    complete_seq_total = []
    for seq_record in gb_seq:
        complete_seq = str(seq_record.seq)
        complete_seq_total.append(complete_seq)
        complete_ana = ">" + seq_record.id + ":" + " " + seq_record.description + "\n"
        complete_fasta = format_fasta(complete_ana, complete_seq, 70)
        complete_total.append(complete_fasta)
        
    return complete_total, complete_seq_total    


def get_cds_dna(gb_file, complete_seq_total):

    gb_seq = SeqIO.parse(open(gb_file), "genbank")
    
    # The CDS dna sequence was extracted and formatted as fasta
    cds_num = 1
    cds_fasta = ""
    for i, seq_record in enumerate(gb_seq):
        for ele in seq_record.features:
            if ele.type == "CDS":
                cds_seq = ""
                if list(ele.qualifiers)[0] == 'gene':
                    cds_ana = ">" + seq_record.id + "_cds_"  + "_" + str(cds_num) + \
                              " [gene=" + ele.qualifiers['gene'][0] + "]" + "[inference=" + ele.qualifiers['inference'][1] + "]" + \
                              " [locus_tag=" + ele.qualifiers['locus_tag'][0] + "]" + " [protein=" + ele.qualifiers['product'][0] + "]" + \
                              " [gbkey=CDS]\n"
                else:
                    cds_ana = ">" + seq_record.id + "_cds_"  + "_" + str(cds_num) + \
                              " [locus_tag=" + ele.qualifiers['locus_tag'][0] + "]" + " [protein=" + ele.qualifiers['product'][0] + "]" + \
                              " [gbkey=CDS]\n"
                cds_num += 1
                for ele1 in ele.location.parts:
                    cds_seq += complete_seq_total[i][ele1.start:ele1.end]
                cds_fasta += format_fasta(cds_ana, cds_seq, 70)
    return cds_fasta


def get_cds_protein(gb_file, f_cds):
    """
    Extract the cds sequence and its complete sequence from the genbank file
    :param gb_file: genbank file path
    :param f_cds: Whether to obtain only one CDS sequence
    :return: CDS sequence in fasta format, complete sequence in fasta format
    """
    gb_seq = SeqIO.parse(open(gb_file), "genbank")

    # The CDS protein sequence was extracted and the format was fasta
    cds_num = 1
    cds_fasta = ""
    for i, seq_record in enumerate(gb_seq):
        for ele in seq_record.features:
            if ele.type == "CDS":
                cds_seq = ""
                if list(ele.qualifiers)[0] == 'gene':
                    cds_ana = ">" + seq_record.id + "_cds_"  + "_" + str(cds_num) + \
                            " [gene=" + ele.qualifiers['gene'][0] + "]" + "[inference=" + ele.qualifiers['inference'][1] + "]" + \
                            " [locus_tag=" + ele.qualifiers['locus_tag'][0] + "]" + " [protein=" + ele.qualifiers['product'][0] + "]" + \
                            " [gbkey=CDS]\n"
                else:
                    cds_ana = ">" + seq_record.id + "_cds_"  + "_" + str(cds_num) + \
                            " [locus_tag=" + ele.qualifiers['locus_tag'][0] + "]" + " [protein=" + ele.qualifiers['product'][0] + "]" + \
                            " [gbkey=CDS]\n"
                cds_num += 1
                cds_fasta += format_fasta(cds_ana, ele.qualifiers['translation'][0], 70)
                if (f_cds):
                    break
    return cds_fasta

def process_blank(path, walkfile):
    #Handle bp Spaces in gbk comment files
    for ls in walkfile:
        res_dir = path + os.sep + ls
        for file in os.listdir(res_dir):
            if file[-3:] == 'gbk' or file[-3:] == 'gbf':
                fo = open(res_dir + os.sep + file, 'r+')
                flist = fo.readlines()
                for k, line in enumerate(flist):
                    if line[:5] == 'LOCUS' and 'length' in line :
                        line = line[:line.find('length')] + line[line.find('length'):line.find('cov')].replace(' ', '') + line[line.find('cov'):]
                        res = line.split('_')
                        num = res[res.index('length') + 1]

                        index_list = []
                        index = line.find(num)
                        while index != -1:
                            index_list.append(index)
                            index = line.find(num, index + 1)
                        
                        line = line[:index_list[0]] + line[index_list[0]:index_list[-1]].replace(" ", "") + ' ' + line[index_list[-1]:]
                        flist[k] = line

                fo = open(res_dir + os.sep + file, 'w+')
                fo.writelines(flist)
                fo.close()


def cds_dna_protein_generate(inter, path, walkfile, obj):
    #Process the phage annotation file to generate the phage_dna file and phage_protein file
    if obj == 'phage':
        dna = inter + os.sep + 'phage_dna'
        protein = inter + os.sep + 'phage_protein'
        os.mkdir(dna)
        os.mkdir(protein)
    else:
        dna = inter + os.sep + 'host_dna'
        protein = inter + os.sep + 'host_protein'
        os.mkdir(dna)
        os.mkdir(protein)
    for ls in walkfile: 
        #cds_dna_file = opt.dna_output + "cds_" + ls[7:] + ".fasta"
        #cds_protein_file = opt.protein_output + "cds_pro_" + ls[7:] + ".fasta"
        cds_dna_file = dna + os.sep + "cds_" + ls[7:] + ".fasta"
        cds_protein_file = protein + os.sep + "cds_pro_" + ls[7:] + ".fasta"
        # genbank file path
        res_dir = path + os.sep + ls
        cds_dna_file_obj = open(cds_dna_file, "w")
        cds_protein_file_obj = open(cds_protein_file, 'w')
        for file in os.listdir(res_dir):
            if file[-3:] == 'gbk' or file[-3:] == 'gbf':
                complete_fasta, complete_seq_fasta = get_complete(res_dir + os.sep + file)
                cds_dna_fasta = get_cds_dna(res_dir + os.sep + file, complete_seq_fasta)
                cds_dna_file_obj.write(cds_dna_fasta)   
                
                cds_protein_fasta = get_cds_protein(res_dir + os.sep + file, False)
                cds_protein_file_obj.write(cds_protein_fasta)                
                
        cds_dna_file_obj.close()
        cds_protein_file_obj.close()

        
def walkDir(file):
    phage_name = []
    for root, dirs, files in os.walk(file):

        # root 	The path of the currently accessed folder
        # dirs 	list of subdirectories in this folder
        # files 	list of files under this folder
        # traversal file
        for f in dirs:
            phage_name.append(os.path.join(f))
    return phage_name

def walkFile(file):
    phage_name = []
    for root, dirs, files in os.walk(file):

        # root 	The path of the currently accessed folder
        # dirs 	list of subdirectories in this folder
        # files 	list of files under this folder
        # traversal file
        for f in files:
            phage_name.append(os.path.join(f))
    return phage_name  


def walkFile_all_path(file):
    phage_name = []
    for root, dirs, files in os.walk(file):

        # root 	The path of the currently accessed folder
        # dirs 	list of subdirectories in this folder
        # files 	list of files under this folder
        # traversal file
        for f in files:
            phage_name.append(file + os.sep + os.path.join(f))
    return phage_name


def get_align_and_interaction_infor(phage_align, phage_raw_data, host_align, host_raw_data, inter, result, pa_fasta_path, ha_fasta_path,
                                   pi_infor_path, hi_infor_path):
    
    # store html link text content folder
    # (re)create these dirs so re-running after a failed attempt doesn't crash on FileExistsError
    # or mix stale files from a previous partial run with the current one
    for d in (pa_fasta_path, ha_fasta_path, pi_infor_path, hi_infor_path):
        if os.path.exists(d):
            shutil.rmtree(d)
        os.mkdir(d)

    # add the extra infor to html / add the extra column for result table
    result['Phage_optimum_align'] = np.nan
    result['Bacterium_optimum_align'] = np.nan
    result['Phage_interaction_infor_in_raw_data'] = np.nan
    result['Bacterium_interaction_infor_in_raw_data'] = np.nan

    # fill the result extra column and copy align target fasta/interaction infor to output folder
    phage_raw_data = walkFile_all_path(phage_raw_data)
    host_raw_data = walkFile_all_path(host_raw_data)
    dic_p_rd = {}
    dic_h_rd = {}
    # Map by sequence accession (parsed from the fasta headers) rather than by filename stem:
    # phage raw files hold a single record whose accession equals the filename, but bacterium
    # raw files are multi-contig draft assemblies (filename is a genus label), so blast hits
    # come back as a contig accession that never matches the filename-based key.
    for path in phage_raw_data:
        for record in SeqIO.parse(path, "fasta"):
            dic_p_rd[record.id.split('.')[0]] = path
    for path in host_raw_data:
        for record in SeqIO.parse(path, "fasta"):
            dic_h_rd[record.id.split('.')[0]] = path

    dic_p = {}
    for i in range(len(phage_align)):
        if os.stat(phage_align[i]).st_size == 0:
            dic_p[phage_align[i].split('/')[-1][:-4]] = ['not hit']
        else:
            f = open(phage_align[i], 'r')
            line = f.readlines()
            f.close()
            sp = line[0].split('\t')
            hit_id = sp[1].split('.')[0]
            if hit_id not in dic_p_rd:
                print("Warning: no raw phage fasta found for hit id " + hit_id + ", marking as not hit")
                dic_p[phage_align[i].split('/')[-1][:-4]] = ['not hit']
                continue
            dic_p[phage_align[i].split('/')[-1][:-4]] = [hit_id, sp[-2]]
            shutil.copyfile(dic_p_rd[hit_id], pa_fasta_path + os.sep + hit_id + '.fasta')
            inter[inter['phage'] == hit_id].iloc[:, [1,3]].to_csv(pi_infor_path + os.sep + hit_id + '.log', index = False)
            #dic_p[phage_align[i].split('/')[-1][:-4]]


    dic_h = {}
    for i in range(len(host_align)):
        if os.stat(host_align[i]).st_size == 0:
            dic_h[host_align[i].split('/')[-1][:-4]] = ['not hit']
        else:
            f = open(host_align[i], 'r')
            line = f.readlines()
            f.close()
            st = line[0].split('\t')
            hit_id = st[1].split('.')[0]
            if hit_id not in dic_h_rd:
                print("Warning: no raw bacterium fasta found for hit id " + hit_id + ", marking as not hit")
                dic_h[host_align[i].split('/')[-1][:-4]] = ['not hit']
                continue
            dic_h[host_align[i].split('/')[-1][:-4]] = [hit_id, st[-2]]
            shutil.copyfile(dic_h_rd[hit_id], ha_fasta_path + os.sep + hit_id + '.fasta')
            inter[inter['host'] == hit_id].iloc[:, [0,2]].to_csv(hi_infor_path + os.sep + hit_id + '.log', index = False)
            #dic_h[host_align[i].split('/')[-1][:-4]]
        
    for i in range(len(result)):
        if len(dic_p[result.iloc[i,0]]) == 1:
            result.iloc[i, 6] = dic_p[result.iloc[i,0]][0]
        else:
            result.iloc[i, 6] = dic_p[result.iloc[i,0]][0] + '(' + dic_p[result.iloc[i,0]][1] + ')' + ' ' + 'phage_target_align' + os.sep + dic_p[result.iloc[i,0]][0] + '.fasta'
            result.iloc[i, 8] = dic_p[result.iloc[i,0]][0] + '_interact_bacteria' + ' ' + 'phage_target_interaction' + os.sep + dic_p[result.iloc[i,0]][0] + '.log'
        
        if len(dic_h[result.iloc[i,1]]) == 1:
            result.iloc[i, 7] = dic_h[result.iloc[i,1]][0]
        else:
            result.iloc[i, 7] = dic_h[result.iloc[i,1]][0] + '(' + dic_h[result.iloc[i,1]][1] + ')' + ' ' + 'host_target_align' + os.sep + dic_h[result.iloc[i,1]][0] + '.fasta'
            result.iloc[i, 9] = dic_h[result.iloc[i,1]][0] + '_interact_phage' + ' ' + 'host_target_interaction' + os.sep + dic_h[result.iloc[i,1]][0] + '.log'
            
    return result
            
            
#html function
style_applied = '''
        body{
			font-family: verdana,arial,sans-serif;
			font-size:11px;
		}
		table.gridtable {
			color: #333333;
			border-width: 1px;
			border-color: #666666;
			border-collapse: collapse;
            font-size:11px;
		}
		table.gridtable th {
			border-width: 1px;
			padding: 8px;
			border-style: solid;
			border-color: #666666;
			background-color: #DDEBF7;
		}
		table.gridtable td {
			border-width: 1px;
			padding: 8px;
			border-style: solid;
			border-color: #666666;
			background-color: #ffffff;
			text-align:center;
		}
		table.gridtable td.key_gene_output {
            font-weight:bold;
			color:#ED5F5F;
		}
		table.gridtable td.wgs_output {
			font-weight:bold;
			color:green;
		}
		li {
			margin-top:5px;
		}
        div{
            font-family: arial;
            font-size:20px;
            color:#0B610B;
            margin-top:15px;
        }
    '''

_TABLE_HEAD_HTML = (
    '<tr><th style="background-color:white"></th><th>Phage</th><th>Bacterium</th>'
    '<th>Key_gene_output</th><th>Wgs_output</th><th>Phage_key_gene_num</th><th>Host_key_gene_num</th>'
    '<th>Phage_optimum_align</th><th>Bacterium_optimum_align</th>'
    '<th>Phage_interaction_infor_in_raw_data</th><th>Bacterium_interaction_infor_in_raw_data</th></tr>'
)

_ROW_TEMPLATE = (
    '<tr><td>{0}</td><td>{1}</td><td>{2}</td>'
    '<td class="key_gene_output">{3}</td><td class="wgs_output">{4}</td>'
    '<td>{5}</td><td>{6}</td>'
    '<td><a href="{8}">{7}</a></td><td><a href="{10}">{9}</a></td>'
    '<td><a href="{12}">{11}</a></td><td><a href="{14}">{13}</a></td></tr>'
)

def generate_html_report(res, path):
    # Building one dominate tag object per table cell (~10 per row) does not
    # scale to hundreds of thousands of rows: the object graph and its final
    # render() traversal dominated the runtime for large result sets. Writing
    # the markup directly as strings, in chunks, keeps the same visual output
    # while avoiding that per-cell object overhead and bounding peak memory.
    col = [res.iloc[:, i].tolist() for i in range(10)]

    with open(path + os.sep + 'result.html', 'w') as f:
        f.write('<html><head><style type="text/css">' + style_applied + '</style></head><body>')
        f.write('<div id="hello"><b><font>Hi Users,</font></b>'
                '<b><font>this is the prediction result of DeepPBI-KG</font></b></div>')
        f.write('<div id="test case result"><table class="gridtable"><tbody>')
        f.write(_TABLE_HEAD_HTML)

        buf = []
        for i in range(len(res)):
            p6, p7, p8, p9 = col[6][i], col[7][i], col[8][i], col[9][i]
            if p6 == 'not hit' and p7 == 'not hit':
                plink_text, plink_url, hlink_text, hlink_url = 'not hit', ' ', 'not hit', ' '
                pilink_text, pilink_url, hilink_text, hilink_url = ' ', ' ', ' ', ' '
            elif p6 == 'not hit' and p7 != 'not hit':
                plink_text, plink_url = 'not hit', ' '
                hlink_text, hlink_url = p7.split(' ')[0], p7.split(' ')[1]
                pilink_text, pilink_url = ' ', ' '
                hilink_text, hilink_url = p9.split(' ')[0], p9.split(' ')[1]
            elif p6 != 'not hit' and p7 == 'not hit':
                plink_text, plink_url = p6.split(' ')[0], p6.split(' ')[1]
                hlink_text, hlink_url = 'not hit', ' '
                pilink_text, pilink_url = p8.split(' ')[0], p8.split(' ')[1]
                hilink_text, hilink_url = ' ', ' '
            else:
                plink_text, plink_url = p6.split(' ')[0], p6.split(' ')[1]
                hlink_text, hlink_url = p7.split(' ')[0], p7.split(' ')[1]
                pilink_text, pilink_url = p8.split(' ')[0], p8.split(' ')[1]
                hilink_text, hilink_url = p9.split(' ')[0], p9.split(' ')[1]

            buf.append(_ROW_TEMPLATE.format(
                i + 1, col[0][i], col[1][i], col[2][i], col[3][i], col[4][i], col[5][i],
                plink_text, plink_url, hlink_text, hlink_url,
                pilink_text, pilink_url, hilink_text, hilink_url,
            ))
            if len(buf) >= 5000:
                f.write(''.join(buf))
                buf = []
        if buf:
            f.write(''.join(buf))

        f.write('</tbody></table></div>')
        f.write('<br><p>** Feel free to connect 20210700125@fudan.edu.cn if you have any question. **</p>')
        f.write('</body></html>')


if __name__ == '__main__':
    # The parameters are defined and encapsulated
    parser = argparse.ArgumentParser()
    parser.add_argument('--phage_annotation', type=str, help = 'input phage annotation file path')
    parser.add_argument('--bacterium_annotation', type=str, help = 'input bacterium annotation file path')
    parser.add_argument('--output', type=str, help='predict result output file path')
    parser.add_argument('--model', type=str, help='model relevant parameter file path')
    parser.add_argument('--template', type=str, help='Temporary folder for intermediate results')
    parser.add_argument('--phage_raw_data', type=str, help = 'Our collected raw phage fasta dataset')
    parser.add_argument('--bacterium_raw_data', type=str, help = 'Our collected raw bacterium fasta dataset')
    parser.add_argument('--phage_align_res', type=str, help = 'the folder of phage alignment result from the first step')
    parser.add_argument('--bacterium_align_res', type=str, help = 'the folder of bacterium alignment result from the first step')
    opt = parser.parse_args()

    # input the annotation file path for phage and bacterium
    phage_walkfile = walkDir(opt.phage_annotation)
    host_walkfile = walkDir(opt.bacterium_annotation)
    
    fea_path = opt.template + os.sep + 'feature_file'
    
    # Condicional para saltar la extraccion si la carpeta ya existe
    if os.path.exists(fea_path):
        print('Carpeta template detectada. Saltando generacion de caracteristicas...')
        from collections import defaultdict
        phage_kg_num = defaultdict(int)
        host_kg_num = defaultdict(int)
    else:
        #Generate dna protein sequence file of phage and bacterium
        cds_dna_protein_generate(opt.template, opt.phage_annotation, phage_walkfile, 'phage')
        cds_dna_protein_generate(opt.template, opt.bacterium_annotation, host_walkfile, 'host')
        print('phage and bacterium cds dna(protein) fasta file process complete!')
        
        # Features are generated from the dna protein sequence file
        os.mkdir(fea_path)
        phage_kg_num = fg(opt.template + os.sep + 'phage_dna', opt.template + os.sep + 'phage_protein', opt.model + os.sep + 'phage_key_gene_0.0001.csv', 'key_gene', 'phage', fea_path)
        fg(opt.template + os.sep + 'phage_dna', opt.template + os.sep + 'phage_protein', opt.model + os.sep + 'phage_key_gene_0.0001.csv', 'wgs', 'phage', fea_path)
        host_kg_num = fg(opt.template + os.sep + 'host_dna', opt.template + os.sep + 'host_protein', opt.model + os.sep + 'host_key_gene_0.00004.csv', 'key_gene', 'bacterium', fea_path)
        fg(opt.template + os.sep + 'host_dna', opt.template + os.sep + 'host_protein', opt.model + os.sep + 'host_key_gene_0.00004.csv', 'wgs', 'bacterium', fea_path)
        print('phage and bacterium feature generation complete!')
    
    #input feature, use the model to predict   
    kg_out = predict(opt.template + os.sep + 'feature_file/phage_key_gene_feature.csv', opt.template + os.sep + 'feature_file/bacterium_key_gene_feature.csv', opt.model + os.sep + 'key_gene_scaler.pkl', opt.model + os.sep + 'key_gene_model.pth')
    wgs_out = predict(opt.template + os.sep + 'feature_file/phage_wgs_feature.csv', opt.template + os.sep + 'feature_file/bacterium_wgs_feature.csv', opt.model + os.sep + 'wgs_scaler.pkl', opt.model + os.sep + 'wgs_model.pth')
    
    phage = wgs_out.index.tolist()
    host = wgs_out.columns.tolist()

    # build the phage x host cross product and fill in outputs via vectorized index alignment
    # instead of an O(n^2) boolean-mask lookup per pair (which never finished for real dataset sizes)
    result = pd.DataFrame(list(itertools.product(phage, host)), columns = ['phage', 'bacterium'])
    result = result.set_index(['phage', 'bacterium'])

    kg_long = kg_out.stack()
    kg_long.index.names = ['phage', 'bacterium']
    wgs_long = wgs_out.stack()
    wgs_long.index.names = ['phage', 'bacterium']

    result['key_gene_output'] = kg_long
    result['wgs_output'] = wgs_long
    result = result.reset_index()

    result['phage_key_gene_num'] = result['phage'].map(phage_kg_num).fillna(0).astype(int)
    result['host_key_gene_num'] = result['bacterium'].map(host_kg_num).fillna(0).astype(int)
    result = result[['phage', 'bacterium', 'key_gene_output', 'wgs_output', 'phage_key_gene_num', 'host_key_gene_num']]
    result.to_csv(opt.output + os.sep + 'result.csv', index = False)
    
    # process phage(bacteria) alignment result and interaction infor
    phage_align = walkFile_all_path(opt.phage_align_res)
    host_align = walkFile_all_path(opt.bacterium_align_res)
    inter = pd.read_csv(opt.model + os.sep + 'interaction_ref_infor.csv', header = 0)
    result = get_align_and_interaction_infor(phage_align, opt.phage_raw_data, host_align, opt.bacterium_raw_data, inter, result, 
                                             opt.output + os.sep + 'phage_target_align', opt.output + os.sep + 'host_target_align', 
                                             opt.output + os.sep + 'phage_target_interaction', opt.output + os.sep + 'host_target_interaction')
    
    # generate html result
    generate_html_report(result, opt.output)
    print('predict complete!')
