#!/usr/bin/env python3.7

import numpy as np
from Bio import SeqIO
import math
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.signal import argrelextrema
from pathlib import Path

def seq_extract(in_file): 
    """Extrct the sequence from the fasta file
    
    Args: 
        in_file [str]: Input file address
        
    Returns: 
        seqs [1d list]: Protein sequences lisst
        
    """
    record = SeqIO.parse(in_file, "fasta")
    seqs = [i.seq for i in record]
    return seqs

class Analysis: 
    """Class will calculate conservation, visualize, and write PyMol script.
    
    Args: 
        seq [1d list]: Sequences of amino acids from fasta file 
    
    """  
    def seq2np(self): 
        """"Turn the sequence into numpy S1 array for calculations later. 
        
        Args: 
            seq [2d list]: List of lists that contain sequences 
        
        Returns: 
            np array [2d np array]: Np array that turns the chars into bytes
            
        """
        
        return np.asarray(self.seq, dtype='S1')
              
    def __init__(self, seq, pdb_id): 
        """Initialize class Analysis. 
        
        Convert the 2d list to np array and make accessible
        globally for other other functions. 
        
        Args: 
            seq [list of lists]: 2d list of sequences
            pdb_id [str]: PDB ID that will be written in the 
                    final PyMol script file.  
            
        """   
        self.seq = seq
        self.pdb_id = pdb_id
        
    def _shannon(self, array): 
        """Calculate Shannon Entropy vertically via loop. 
        
        Args: 
            array [nd array]: 1d array of sequences from all the species
        
        Returns: 
            entropy [nd float array]: Per column value for each position vertically
        
        """
        
        array = str(array)
        entropy = 0 
        
        for x in array: 
            char_count = np.char.count(array, sub=x)

            if x == "-": 
                # Penalize columns with all - values
                p_x = float(char_count/len(array)*2) 
                
            else: p_x = float(char_count/len(array))
            
            if p_x > 0: 
                entropy += (- p_x*math.log(p_x, 2)/20)
        return entropy
        
    def conservation_score(self, np_seq): 
        """Calculate the Shannon Entropy vertically
        for each position in the amino acid msa sequence.
        
        Args: 
            np_seq [Numpy nd array]: Np array of sequences
            
        Returns: 
            np apply array [nd float array]: Calculate conservation 
            scores vertically into a float nd array   
        """
        
        return np.apply_along_axis(self._shannon, 0, np_seq)      

    def normalize_data(self, ent_list): 
        """Takes the entropy array and normalizes the data. 
        
        Args: 
            ent_list [Nd array]: Entropy float array
            
        Returns: 
            Normalized list [nd array]: Values between -1 and 1
            
        """
        
        return (2.*(ent_list - \
            np.min(ent_list))/np.ptp(ent_list)-1)       
        
    
    def find_local_minima(self, data): 
        """Find local minima that are lower than mean of data
        
        Args: 
            data [nd float array]: Normalized conservation data
            
        Returns: 
            polished_minima [list]: List that contains minima 
            lower than the mean.
            
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
        """Given the minima and data, look for consecutive
        set (word size of 4) of values that are 1/4th of the mean. 
        
        Args; 
            data [nd float array]: Normalized conservation data
            minima [1d list]: List of minima lower than the mean 
                        
        Returns: 
            pos_motif [1d list]: Possible motifs that have 4 consecutive 
                lower scores than 1/4th of the mean, only the 1st index. 
                
            pos [1d list]: Positions of these motifs.        
             
        """
        
        pos_motif = []
        pos = []
        threshold = float(np.mean(data)/2)
        for i in minima: 
            if i - 4 > 0: 
                motif_stretch = data[i - 1 : i + 1]
                if np.all(motif_stretch < threshold): 
                    # Take the first value. 
                    pos_motif.append(motif_stretch[0])
                    pos.append(i)
                       
        return (pos_motif, pos)
    
    def pymol_script_writer(self, out_file, pos): 
        """Motifs that are found are then written a txt file.
        Thhis is then run on PyMol by typing the following on terminal: 
            @pymol_script.txt 
        
        Args: 
            out_file [str]: Address for where the file should be written
            pos [1d list]: Posiitons of the motifs as list
            
        """
        
        path = Path(out_file)
        with open(path, "w") as file: 
            file.write(f"fetch {self.pdb_id}\n\n")
            for i in range(len(pos)): 
                file.write(f"create mot{pos[i]}, resi {pos[i]}-{pos[i]+4} \n")
                
            file.write("\nhide_all\n")
            file.write("\n")
            for i in range(len(pos)): 
                file.write(f"show cartoon, resi{pos[i]}\n")
            file.write("\n")
            for i in range(len(pos)): 
                file.write(f"color red, mot{pos[i]}\n")
        

def main(): 
    seq = seq_extract("/home/nadzhou/Desktop/clustal.fasta")
    seq = [[x for x in y] for y in seq]

    c = Analysis(seq)

    c_ent = c.conservation_score()
    norm_data = c.normalize_data(c_ent)
    norm_data_len = [x for x in range(len(norm_data))]
    minima = c.find_local_minima(norm_data)

    print(norm_data)
    pos_motif, pos = c.find_motif(norm_data, minima)
    
if __name__ == "__main__": 
    main()

