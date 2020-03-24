#!/usr/bin/env python3.7


import numpy as np
from Bio import SeqIO
import math
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.signal import argrelextrema
from scipy.stats import chi2

def seq_extract(in_file): 
    """Extrct the sequence from the fasta file
    Args: 
        input file [str]
    REturns: 
        Sequence [1d list]
    """
    record = SeqIO.parse(in_file, "fasta")
    seqs = []
    for i in record: 
        seqs.append(str(i.seq))
    return seqs
class Analysis: 
    """
    Class to analyze the multiple sequence alignment. 
    Implements the Shannon Entropy, chi-squared and 
    per-column variation analysis algorithms. 
    """        
    def __init__(self, seq): 
        """
        Initialize class Analysis. 
        
        Args: 
            Sequence list [list of lists]
            
        Returns: 
            Conservation score 
            Chi-squared table 
            etc etc
        """
        def seq2np(seq): 
            return np.asarray(seq, dtype=object)
        
        np_seq = seq2np(seq)
        self.seq = seq
        self.np_seq = np_seq   
    

    def conservation_score(self): 
        """
        Calculate the Shannon Entropy vertically
        for each position in the amino acid msa sequence.
        
        Args: 
            Numpy sequence [nd array]
        Returns: 
            Vertical-orientation conservation score [nd array]    
        """

        def _shannon(array): 
            array = str(array)
            entropy = 0 
            for x in array: 
                char_count = np.char.count(array, sub=x)

                if x == "-": 
                    p_x = float(char_count/len(array)*2) 
                    
                else: p_x = float(char_count/len(array))
                if p_x > 0: 
                    entropy += (- p_x*math.log(p_x, 2)/20)
            return entropy
        return np.apply_along_axis(_shannon, 0, self.np_seq)      

    def normalize_data(self, ent_list): 
        """Takes the entropy array and normalizes the data. 
        
        Args: 
            Shannon Entropy array [numpy ndarray]
        Returns: 
            Normalized numpy array [numpy ndarray]
        """
        return (2.*(ent_list - \
            np.min(ent_list))/np.ptp(ent_list)-1)       
        
    def apply_chisquare(self, norm_list): 
        return chi2(norm_list)
    
    def find_local_minima(self, data): 
        """
        Find the local minima of the given data
        and then only return a list of values that is 
        less than the mean of the data set. 
        
        Args: 
            Normalized np array [nd array]
        Returns: 
            Local minma less than mean [list]
        """
        
        local_minima = argrelextrema(data, np.less)
        data_mean = np.mean(data)
        polished_minima = []
        for i in local_minima: 
            for j in i: 
                if data[j] < data_mean: 
                    polished_minima.append(j)
        
        return polished_minima

    def find_motif(self, data, minima): 
        """
        In a given stretch of 4 word size,, find possible motifs
        by seeing if all the values are less than a threshold, which 
        is 1/4th of the mean of the normalized data. 
        
        Args; 
            Normalized data [nd array]
            Local minima [nd array]
            
        Returns: 
            List of possible motifs [list]
        """
        
        pos_motif = []
        pos = []
        threshold = float(np.mean(data)/2)
        for i in minima: 
            if i - 4 > 0: 
                motif_stretch = data[i : i + 2]
                if np.all(motif_stretch < threshold): 
                    # Take the first value. 
                    pos_motif.append(motif_stretch[0])
                    pos.append(i)
                       
        return (pos_motif, pos)
        
seq = seq_extract("/home/nadzhou/Desktop/clustal.fasta")
seq = [[x for x in y] for y in seq]

c = Analysis(seq)

c_ent = c.conservation_score()
norm_data = c.normalize_data(c_ent)
norm_data_len = [x for x in range(len(norm_data))]
minima = c.find_local_minima(norm_data)

pos_motif, pos = c.find_motif(norm_data, minima)

print(pos)
sns.lineplot(x=norm_data_len, y=norm_data)
plt.show()
