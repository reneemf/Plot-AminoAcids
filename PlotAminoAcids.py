MWData = {'A': 89.0941, 'C': 121.1541, 'D': 133.1039, 'E': 147.1308, 'F': 165.1919, 'G': 75.0672, 'H': 155.1564,
'I': 131.1747, 'K': 146.1894, 'L': 131.1747, 'M': 149.2079, 'N': 132.1191, 'P': 115.1320, 'Q': 146.1460, 'R': 174.2028,
'S': 105.0935, 'T': 119.1204, 'V': 117.1479, 'W': 204.2285, 'Y': 181.1913, 'water': 18.0153, 'X': 75.0672}

fasta_file = raw_input('Enter a file path associated with a FASTA-formatted protein sequence: ')
#test path: /Users/reneefonseca/Desktop/sequence.fasta
with open(fasta_file, 'r') as input_file:
    protein_fasta = input_file.read().replace('\n', '')
#print protein_fasta
#print type(protein_fasta)
protein_seq = protein_fasta.split("]")[1]
#print protein_seq
water_count = len(protein_seq) - 1
#print water_count
    
class Sequence:
    def __init__(self, sequence):
        self.sequence=sequence
            
class Protein(Sequence):
    def get_mw(self):
        amino_mw = 0
        for n in self.sequence:
            if n in MWData:
                amino_mw += MWData[n]
        water = water_count * MWData["water"] 
        amino_mw = amino_mw - water   
        return amino_mw
    def count_aa(self):
        aa_dict = {}
        for n in self.sequence:
            if n in MWData:
                aa_dict[n] = 0
        for n in self.sequence:
            if n in MWData:
                aa_dict[n] += 1
        return aa_dict

protein_obj=Protein(protein_seq)
print "Molecular Mass:",protein_obj.get_mw()
#print protein_obj.count_aa()
aa_prop = protein_obj.count_aa()
print "Amino acid counts:",aa_prop

#Convert Amino Acid Dictionary counts into proportions 
def prop_dict(aa_count):
    total = 0
    for i in aa_count:
        total = total + aa_count[i]
    for j in aa_count:
        aa_count[j] = (float)(aa_count[j])/total
    return aa_count
#print prop_dict(aa_prop)
plot_aa = prop_dict(aa_prop)
print "Amino acid proportions:",plot_aa

#Sort Amino Acid Dictionary into ordered list of tuples
sorted(plot_aa.values())
sorted(plot_aa, key=plot_aa.get)
plot_ordered = sorted(plot_aa.items(), key=lambda x:x[1], reverse=True)
#print plot_ordered

#Bar plot of Amino Acid proportions
import pandas as pd
import matplotlib.pyplot as plt

df = pd.DataFrame(plot_ordered, columns=['Amino Acids','Proportion'])
df.plot(x='Amino Acids', kind='bar', legend=False)
plt.ylabel("Proportion of Each A.A. in the Sequence")
plt.savefig('AAproportions')
plt.show()
