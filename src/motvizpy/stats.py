#!/usr/bin/env python3.7

import numpy as np
import pandas as pd
from Bio import SeqIO

from scipy.signal import argrelextrema
from pathlib import Path
from collections import Counter

from typing import List, Dict, Sequence


class Analysis: 
    """Class will calculate conservation, visualize, and write PyMol script.
    """  

    def seq2np(self) -> Sequence:     
        return np.asarray(self.seq, dtype='S1')
              
    def __init__(self, seq): 
        self.seq = seq
        self.np_seq = self.seq2np()
        
    def _shannon(self, array: List) -> Sequence: 
        """Apply the Shannon Entropy for each vertical amino acid column
        """
        
        aa_count = Counter(array)
        sum_aa_s = sum(aa_count.values())
        
        pA = 1
        for aa, aa_freq in aa_count.items(): 
            pA *= (aa_freq / sum_aa_s)
        
        return -(np.sum(pA*np.log2(pA)))
        
        
    def conservation_score(self) -> Sequence:       
        """Calculate the Shannon Entropy score. 
        """
        return np.apply_along_axis(self._shannon, 0, self.np_seq)      


    def normalize_data(self, ent_list: List[float]) -> List[float]: 
        """Normalize the data between 1 and -1
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

        for local_minimum in local_minima: 
            for value in local_minimum: 

                if data[value] < data_mean: 
                    polished_minima.append(value)
        
        return polished_minima


    def find_motif(self, data, minima, threshold): 
        """Find possible motifs in the given dataset            
        """
        conserved_posits = {}

        for min_val in minima:         

            for stretch in range(6, 30): 

                motif_stretch = data[min_val : min_val + stretch]

                if np.all(motif_stretch < threshold): 
                    conserved_posits[min_val] = min_val + stretch

                print(f"Maximum stretch of amino acids: {stretch}")

                    
        return conserved_posits
    
        
    def pymol_script_writer(self, out_file: str, pos: List): 
        """Write PyMol script for the given possibly important motifs   
        """
        path = Path(out_file)

        with open(path, "w") as file: 
            file.write("fetch 1fmw\n\n")

            for motif_start, motif_end in pos.items():  
                file.write(f"create mot{motif_start}, resi {motif_start}-{motif_end} \n")
                
            file.write("\nhide_all\n")
            file.write("\n")

            for motif_start, motif_end in pos.items():  
                file.write(f"show cartoon, resi{motif_start}\n")
            file.write("\n")
            
            for motif_start, motif_end in pos.items():  
                file.write(f"color red, mot{motif_start}\n")

        

    def find_conservation(self): 
        return np.apply_along_axis(self._find_conservation_depth, 0, self.np_seq)


    def _find_conservation_depth(self, array): 
        unique_aas, aa_counts = np.unique(array, return_counts=True)

        aa_frequency_dict = dict(zip(unique_aas, aa_counts))

        most_conserved_aa = np.max(aa_counts)

        for aa, aa_count in aa_frequency_dict.items(): 
            if aa_count == most_conserved_aa: 

                return aa_count / sum(aa_counts)


def plotter(norm_data): 
    norm_data_len = np.arange(1, len(norm_data) + 1)

    fig = sns.lineplot(norm_data_len, norm_data)
    plt.title("Conservation score per amino acid position")
    plt.xlabel("Amino acid position")
    plt.ylabel("Normalized conservation score")
    return fig


def main(): 
    orig_file = Path("/home/nadzhou/Desktop/6x2c/6x2c.fasta")
    seq_dict = seq_extract("/home/nadzhou/Desktop/6x2c/aligned_seq.fasta", "fasta")
    seq = [[x for x in y] for y in seq_dict.values()]
    c = Analysis(seq, "1xef")

    c_ent = c.conservation_score(c.seq2np())


    norm_data = c.normalize_data(c_ent)


    #norm_data = c.moving_average(norm_data)

    norm_data_len = [i for i,_ in enumerate(norm_data)]
    # minima = c.find_local_minima(norm_data)

    # cons_data = {num : data for num, data in zip(norm_data_len, norm_data)}

    # c.csv_writer("/home/nadzhou/Desktop/results.csv", cons_data)
    
    
    fig = plotter(norm_data)

    plt.show()



if __name__ == "__main__": 
    main()




    # def csv_writer(self, file, cons_data): 
    #     """Write CSV file 
            
    #     """
    #     print(cons_data)
        
    #     df = pd.DataFrame.from_dict(cons_data, orient="index")
    #     df.index.name = "Amino acid position"
    #     df.columns = ["Conservation score"]
    #     df.to_csv(file)
