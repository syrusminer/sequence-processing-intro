import Bio
from Bio import SeqIO
from Bio.Seq import Seq
import sys
import os
import subprocess


#Opening FASTA file and returning open file for use later.
def OpenFile():
    if len(sys.argv) < 2:
        FASTAfile = input("Please enter your fasta file name: ")
    else:
        FASTAfile = sys.argv[1]

    return open(FASTAfile, 'r')

#Separates the DNA and the Headers and returns each as a new variable.
def FindDNA(theDNA):
    num = 1
    sequencesFU = {}
    dates = {}
    for record in SeqIO.parse(theDNA, "fasta"):
        DNAsequence = (str(record.seq))
        everydate = (str(record.id))
        sequencesFU[num] = DNAsequence
        dates[num] = everydate 
        num +=1
    return sequencesFU,dates

#Simply returns the number of sequences in the file.    
def Samples(theseqs):
    ttlsamples = list(theseqs)[-1]
    return {ttlsamples}

#Counting unique dates
def UniqueDates(dates):
    uqdates = set(dates.values())
    thenumber_ofuniquedates = len(uqdates)
    return thenumber_ofuniquedates

#This one is creating a list for every position.
#At each position it turns the list into a set, if the
#set is larger than one it is a SNP since no duplicates
#are placed into a set. (also used for proteins)
def findingSNP(DNAstrands):
    snp_counter = 0
    strand_length = len(next(iter(DNAstrands.values())))

    for i in range(strand_length):
        nucleotides = []

        for strand in DNAstrands.values():
            currentNuc = strand[i]
            nucleotides.append(currentNuc)

        if len(set(nucleotides)) > 1:
            snp_counter +=1

    return snp_counter


def makingProteins(theSequence):
    counter = 1
    proteinDictionary = {}
    for i in theSequence.values():
        dna_seq = Seq(i)
        protein = dna_seq.translate()
        proteinDictionary[counter] = protein
        counter += 1
    return proteinDictionary

#This function allows the user to choose the name of their new protein file
#The new name is put as the file name and the file is opened, a for loop is used
#with the zip method to allow for both variables to be placed into the new file.
def protesFASTA(seq,date):
    desired_name = input("What do you want to call your new protein file? ")
    with open(desired_name,'w') as theFile:
        for i,e in zip(seq.values(),date.values()):
            theFile.write(f">{e}\n{i}\n")

def MillionDirectories(dna,proteins):
    counter = 1
    for i,e in zip(dna.values(),proteins.values()):
        new_directory = f"Sequence_{counter}"
        os.makedirs(new_directory, exist_ok= True)

        dna_sequence = os.path.join(new_directory, "DNA_sequence.txt")
        with open(dna_sequence, "w") as dna_file:
            dna_file.write(f"DNA Sequence: {i}")

        protein_sequence = os.path.join(new_directory, "Protein_sequence.txt")
        with open(protein_sequence,"w") as protein_file:
            protein_file.write(f"Protein Sequence: {e}")

        counter += 1

def main():
    FASTAfu = OpenFile()
    
    # Use sed to extract unique sampling dates
    unique_dates_process = subprocess.run(['sed', 's/Sample\([0-9]*\)_\([0-9]*-[A-Za-z]*\) \([0-9]*\) \([a-z]*\)/\2 /', FASTAfu], capture_output=True, text=True)
    unique_dates = unique_dates_process.stdout.strip()

    sequences, dates = FindDNA(FASTAfu)
    totalsamples = str(Samples(sequences))
    num_dates = UniqueDates(unique_dates)
    SNPs = findingSNP(sequences)
    Proteins = makingProteins(sequences)
    proteinSNP = findingSNP(Proteins)

    with open("log.txt", "w") as new_file:
        new_file.write(f"Total samples: {totalsamples}\n")
        new_file.write(f"Unique sampling dates: {num_dates}\n")
        new_file.write(f"Your FASTA file contains:\n")
        new_file.write(f"     DNA variable sites:{SNPs}\n")
        new_file.write(f"     Protein variable sites:{proteinSNP}")

    protesFASTA(Proteins, dates)
    MillionDirectories(sequences, Proteins)

    FASTAfu.close()

if __name__ == "__main__":
    main()