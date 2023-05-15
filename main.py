""" Small project that making a txt file with some
protein data after giving amino, rna or dna sequence
 works only for sequences from NCBI database """


# func that translate dna to amino acids
def translation_dna_to_amino(seq):
    nucleoToAmino = {'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
                     'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
                     'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
                     'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
                     'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
                     'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
                     'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
                     'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
                     'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
                     'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
                     'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
                     'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
                     'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
                     'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
                     'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
                     'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W', }
    protein = ""
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein += nucleoToAmino[codon]
    else:
        print('your sequence is not dividing by 3.')
    return protein


# function to check type of sequence for example. Is it RNA DNA or AMINO


def checking_sequence(sequence):
    if set(sequence).issubset(set(nucleoOfDNA)):
        TypeOfSeq = 'DNA'
    elif set(sequence).issubset(set(nucleoOfRNA)):
        TypeOfSeq = 'RNA'
    elif set(sequence).issubset(set(massesOfAmino)):
        TypeOfSeq = 'amino'
    else:
        TypeOfSeq = 'Error'
    return TypeOfSeq

# new function removing special symbols from string


def removing_symbols(string):
    newString = ""
    listOfSpecialSymbols = ['!', '@', '#', '$', '%', '^', '&', '*',
                            '(', ')', "'", "[", ']', '<', '.', '>', "/", '\n']
    for letter in string:
        if letter not in listOfSpecialSymbols:
            newString += letter
    return newString


# fuction that getting name of protein from fasta file


def get_name_of_protein(string):
    firstLine = string[0]
    name = removing_symbols(firstLine)
    return name


# changing fasta file to just and aminoacids sequence


def editing_fasta_to_read_format(arg):
    arg.pop(0)
    newFastaFormat = removing_symbols(arg)
    return newFastaFormat


def sum_masses(sequence):  # function that summarize masses of each aminoacid
    sumOfMassess = 0
    for amino in sequence:
        for mass in amino:
            sumOfMassess += massesOfAmino[mass]

    return sumOfMassess

# function that counting how much amino we have in sequence


def how_many_amino(sequence):
    howMany = 0
    for aminoacid in sequence:
        for symbol in aminoacid:
            if symbol in massesOfAmino:
                howMany += 1
    return howMany


def find_stops(seq):  # finding stop codons
    i = 0
    while i < len(seq) - 2 and seq[i:i + 3] not in ('TAA', 'TAG', 'TGA'):
        i += 3
    if i < len(seq) - 2:
        return i + 3
    else:
        return -1


def find_orfs(seq):  # finding all ORF's in sequence
    orfs = []
    stopIndx = []
    PotentialStops = find_stops(seq[start:])
    if PotentialStops != -1:
        stop = start + PotentialStops
        if stop not in stopIndx:
            orfs.append((PotentialStops, start, stop))
            stopIndx.append(stop)
    orfs = sorted(orfs, reverse=True)
    for i, orf in enumerate(orfs):
        orfs[i] = (orf[1], orf[2] - 3)
    return list(orfs)


nucleoOfRNA = ['U', 'C', 'G', 'A']
nucleoOfDNA = ['A', 'G', 'T', 'C']
massesOfAmino = {
    "A": 71.0788,
    "R": 156.1875,
    "N": 114.1038,
    "D": 115.0886,
    "C": 103.1388,
    "E": 129.1155,
    "Q": 128.1307,
    "G": 57.0519,
    "H": 137.1411,
    "I": 113.1594,
    "L": 113.1594,
    "K": 128.1741,
    "M": 131.1926,
    "F": 147.1766,
    "P": 97.1167,
    "S": 87.0782,
    "T": 101.1051,
    "W": 186.2132,
    "Y": 163.1760,
    "V": 99.1326
}

userSequence = input("Please write a file fasta name with extension: ")
with open(userSequence, "r+", encoding="UTF-8")as file:
    sequenceOfAmino = file.readlines()
    nameOfProt = get_name_of_protein(sequenceOfAmino)
    sequenceAfterTrans = removing_symbols(
        editing_fasta_to_read_format(sequenceOfAmino)).upper()

# for aminoacid sequence
if (checking_sequence(sequenceAfterTrans) == 'amino'):  # it happens when you give aminoacid sequence
    proteinMassKDA = int(sum_masses(sequenceAfterTrans) / 1000)
    howManyAmino = how_many_amino(sequenceAfterTrans)
    with open(f"{nameOfProt}.txt", 'w', encoding="UTF-8") as file:
        file.write(
            f"{nameOfProt} protein mass is  {proteinMassKDA} kDa and contains {howManyAmino} aminoacids")


# for RNA sequence

elif (checking_sequence(sequenceAfterTrans) == 'RNA'):  # it happens when you give RNA sequence
    dnaSequence = sequenceAfterTrans.replace("U", "T")
    start = dnaSequence.find("ATG")
    orf = find_orfs(dnaSequence)[0]
    aminoSeq = translation_dna_to_amino(dnaSequence[orf[0]: orf[1]])
    proteinMassKDA = int(sum_masses(aminoSeq) / 1000)
    howManyAmino = how_many_amino(aminoSeq)
    with open(f"{nameOfProt}.txt", 'w', encoding="UTF-8") as file:
        file.write(
            f"{nameOfProt} protein mass is  {proteinMassKDA} kDa and contains {howManyAmino} aminoacids")


# for DNA sequence
elif (checking_sequence(sequenceAfterTrans) == 'DNA'):  # it happens when you give DNA sequence
    start = sequenceAfterTrans.find("ATG")
    orf = find_orfs(sequenceAfterTrans)[0]
    aminoSeq = translation_dna_to_amino(sequenceAfterTrans[orf[0]: orf[1]])
    proteinMassKDA = int(sum_masses(aminoSeq) / 1000)
    howManyAmino = how_many_amino(aminoSeq)
    with open(f"{nameOfProt}.txt", 'w', encoding="UTF-8") as file:
        file.write(
            f"{nameOfProt} protein mass is  {proteinMassKDA} kDa and contains {howManyAmino} aminoacids")
else:
    print('Error, try again')
