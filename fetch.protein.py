import requests

def fetch_sequence(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    response = requests.get(url)
    return response.text

data = fetch_sequence("P04637")  # TP53
print(data[:200])
from Bio import SeqIO
from io import StringIO

def parse_fasta(fasta_data):
    record = SeqIO.read(StringIO(fasta_data), "fasta")
    return str(record.seq)

sequence = parse_fasta(data)
print("Sequence length:", len(sequence))
from collections import Counter

def aa_composition(seq):
    count = Counter(seq)
    total = len(seq)
    return {aa: round(count[aa]/total, 3) for aa in count}

composition = aa_composition(sequence)
print(composition)
from Bio.SeqUtils import molecular_weight

mw = molecular_weight(sequence, seq_type="protein")
print("Molecular Weight:", mw)
hydrophobic = "AILMFWYV"
score = sum([1 for aa in sequence if aa in hydrophobic]) / len(sequence)

print("Hydrophobicity Score:", round(score, 3))
import matplotlib.pyplot as plt

plt.bar(composition.keys(), composition.values())
plt.title("Amino Acid Composition")
plt.xticks(rotation=90)
plt.show()
# import matplotlib.pyplot as plt
# plt.bar(...)
# plt.show()
import matplotlib.pyplot as plt

plt.bar(composition.keys(), composition.values())
plt.title("Amino Acid Composition")
plt.xticks(rotation=90)
plt.show()


