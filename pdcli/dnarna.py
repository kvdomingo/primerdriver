convert = {
    "U": {
        "U": {
            "U": "F",
            "C": "F",
            "A": "L",
            "G": "L"
        },
        "C": {
            "U": "S",
            "C": "S",
            "A": "S",
            "G": "S"
        },
        "A": {
            "U": "Y",
            "C": "Y",
            "A": "*",
            "G": "*"
        },
        "G": {
            "U": "C",
            "C": "C",
            "A": "*",
            "G": "W"
        }
    },
    "C": {
        "U": {
            "U": "L",
            "C": "L",
            "A": "L",
            "G": "L"
        },
        "C": {
            "U": "P",
            "C": "P",
            "A": "P",
            "G": "P"
        },
        "A": {
            "U": "H",
            "C": "H",
            "A": "Q",
            "G": "Q"
        },
        "G": {
            "U": "R",
            "C": "R",
            "A": "R",
            "G": "R"
        }
    },
    "A": {
        "U": {
            "U": "I",
            "C": "I",
            "A": "I",
            "G": "M"
        },
        "C": {
            "U": "T",
            "C": "T",
            "A": "T",
            "G": "T"
        },
        "A": {
            "U": "N",
            "C": "N",
            "A": "K",
            "G": "K"
        },
        "G": {
            "U": "S",
            "C": "S",
            "A": "R",
            "G": "R"
        }
    },
    "G": {
        "U": {
            "U": "V",
            "C": "V",
            "A": "V",
            "G": "V"
        },
        "C": {
            "U": "A",
            "C": "A",
            "A": "A",
            "G": "A"
        },
        "A": {
            "U": "D",
            "C": "D",
            "A": "E",
            "G": "E"
        },
        "G": {
            "U": "G",
            "C": "G",
            "A": "G",
            "G": "G"
        }
    }
}

dna = input("Enter DNA sequence: ")
dna = dna.upper()
rna = []
for p in range(0,len(dna)):
	if dna[p] == 'A':
			rna.append('A')
	if dna[p] == 'T':
		rna.append('U')
	if dna[p] == 'G':
		rna.append('G')
	if dna[p] == 'C':
		rna.append('C')
else:
	print("The RNA sequence is: ", "".join(rna))

n = 3
out = [(rna[i:i+n]) for i in range(0, len(rna), n) if len(rna[i:i+n])==3] 

protein = []
for i in out:
	protein.append(convert[i[0]][i[1]][i[2]])
print("The protein sequence is: ", "".join(protein))
