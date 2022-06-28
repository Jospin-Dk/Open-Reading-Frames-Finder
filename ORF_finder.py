##################################################################################################
#PACKAGES
##################################################################################################
import argparse
import os.path

##################################################################################################
#SUPPORTING OBJECTS
##################################################################################################
#a dictionary containing all amino acids and their respective codons.
#an extra codon containing N on the 3rd nucleotide position was added for each of the L, V, S, P, T, A, R and G amino acids. for these specific amino acids no matter what the 3rd nucleotide was in the codon, it would still code for these specific amino acids.
DNA_to_amino={"TAA":"stop", "TAG":"stop", "TGA":"stop",
"TTT":"F", "TTC":"F",
"TTA":"L", "TTG":"L", "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L", "CTN":"L",
"ATT":"I", "ATC":"I", "ATA":"I",
"ATG":"M",
"GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V", "CTN":"V",
"TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S", "TCN":"S",
"CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P", "CCN":"P",
"ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T", "ACN":"T",
"GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A", "GCN":"A",
"TAT":"Y", "TAC":"Y",
"CAT":"H", "CAC":"H",
"CAA":"Q","CAG":"Q",
"AAT":"N", "AAC":"N",
"AAA":"K", "AAG":"K",
"GAT":"D","GAC":"D",
"GAA":"E","GAG":"E",
"TGT":"C", "TGC":"C",
"TGG":"W",
"CGT":"R","CGC":"R","CGA":"R", "CGG":"R", "AGA":"R", "AGG":"R", "AGN":"R",
"AGT":"S","AGC":"S",
"GGT":"G", "GGC":"G","GGA":"G","GGG":"G","GGN":"G"}

##################################################################################################
#METHODS
##################################################################################################
#identify ORFs
#takes as input a string containing the entire sequence or the genome from the input file, the start position from which the ORF starts from (could be 0, 1 or 2) and the minimum length of a coding sequence
#the method outputs a dictionary having the starting position of the ATG in the sequence or genome as the key and the coding sequence (ending with a stop codon) as a string as the value
def identifyORF(genome, start_pos,ORF_length=50):
    start_codon=False # will be set to true when a start codon is encountered
    stop_codon=False #will be set to true when a stop codon is encountered
    inframe=False #will be set to true when we are inside a coding sequence (encountering a new ATG within that same sequence won't be considered as a new coding sequence)
    frame={} #dictionary in which coding frames will be stored
    count=0 #keeping track of the number of nucleotides. Will also work as the key for each coding sequence in order to to know where it started in the genome
    ORF="" #storing the translated protein of a certain coding sequence. will be emptied whe a stop codon is encountered
    for ind in range (start_pos, len(genome),3):
        if not((ind+3)>len(genome)): #no index out of bound and enough nucleotides to make a codon
            codon=genome[ind]+genome[ind+1]+genome[ind+2] #creating the codon
            if codon.upper()== "ATG" and inframe==False: #checking for a start codon that starts a new coding sequence
                count=ind
                start_codon=True
                stop_codon=False
                inframe=True
            if (codon.upper() in DNA_to_amino) and (DNA_to_amino[codon.upper()]=="stop"): #checking for a stop codon
                start_codon=False
                inframe=False
                stop_codon=True
                if len(ORF)>=ORF_length: #after encountering a stop codon, check for the length of the coding sequence with respect to the set sequence length by the user
                    frame[count]=ORF
                ORF=""
            if start_codon==True and inframe==True and stop_codon==False:#we are in a coding frame
                if codon.upper() in DNA_to_amino:
                    ORF=ORF + DNA_to_amino[codon.upper()]
                else:
                    ORF=ORF + "X" #X represents an unknown codon that could be could not be translated due to the presence of N in the nucleotide sequence
    return frame

#creates the compliment of the orignal sequence and returns it as a string
def reverseStrand(genome):
    rev_genome=""
    for i in genome:
        if i =="a":
            rev_genome=rev_genome+"t"
        elif i=="t":
            rev_genome=rev_genome+"a"
        elif i=="c":
            rev_genome=rev_genome+"g"
        elif i=="g":
            rev_genome=rev_genome+"c"
        else:
            rev_genome=rev_genome+"n"
    return rev_genome

#method embedded in argparse
#checking whether an input file is a .fasta file.
def isFasta(value):
    if os.path.isfile(value)==False: #checking whether the file exists first
        raise argparse.ArgumentTypeError("file %s does not exist" % value)
    name_extension=value.split("/") #separating the file's name from the path(in case the file is not present in the same working directory as the script)
    name_split=name_extension[len(name_extension)-1].split(".") #separating the file name from the extension
    file_extension=name_split[1]
    if file_extension !="fasta":
        raise argparse.ArgumentTypeError("%s does not have a .fasta extension" % value)
    return value

#method embedded in argparse
#checking if the integer inputted for the -size argument is a positive integer (>0)
def check_positive(value):
    ivalue = int(value)
    if ivalue <= 0:
        raise argparse.ArgumentTypeError("%s is an invalid positive int value" % value)
    return ivalue

#method embedded in argparse
#checking if the output file name specified is a .fa file
def checkFasta(value):
    filename= value.split(".")
    if "fa" != filename[1]:
        raise argparse.ArgumentTypeError("%s does not have a .fa extension" % value)
    return value

#printing the ORFs to the output file
def printOutput(fileName, f_ORFsList, r_ORFsList):
    count=1
    frame=1
    head_name=genome_name.split()
    #printing the forward strand
    for dic in f_ORFsList:
        for key in dic: #coding sequences within ORFs containing an X were disregarded and not outputted to the output file
            if "X" in dic[key]:
                continue
            else:
                output=">"+head_name[0]+"_"+"F"+str(frame)+"_"+str(count)+" "+str(len(dic[key]))+" "+str(key)
                print("{0}\n{1}".format(output,dic[key]),file=of)
                count+=1
        frame+=1
    #printing the reverse strand
    for dic in r_ORFsList:
        for key in dic: #coding sequences within ORFs containing an X were disregarded and not outputted to the output file
            if "X" in dic[key]:
                continue
            else:
                output=">"+head_name[0]+"_"+"F"+str(frame)+"_"+str(count)+" "+str(len(dic[key]))+" "+str(key)
                print("{0}\n{1}".format(output,dic[key]),file=of)
                count+=1
        frame+=1
##################################################################################################
#ARGPARSE
##################################################################################################
parser= argparse.ArgumentParser()
parser.add_argument("input", help="name of a .fasta file containing a genome or a sequence", type=isFasta) #this argument is compulsory and must be mentioned by the user as the 1st argument (after the name of the script)! while the rest are all optional and could be used for extra control of the output
parser.add_argument("-size", help="the minimum length of a coding sequence in an ORF. must be a positive integer", default="50", type=check_positive)
parser.add_argument("-nbF", help="restriction on the number of frames ", choices= ["1","2","3","4","5","6"],default="6") #nbF -> number of ORFs
parser.add_argument("-output", help="name of a .fa file in which the ORFs generated from the genome will be saved", default="output.fa", type=checkFasta)
args=parser.parse_args()

##################################################################################################
#MAIN
##################################################################################################
file = args.input #catching the file name argument
f = open(file, "r")

line = f.readline().rstrip() #first line of the fasta file containing the header of the sequence or genome
nlines = 0
header = line.split(">")
genome_name=header[1] #saving the genome name

########################################################
#READING THE GENOME
########################################################
genome="" #storing the entire sequence or genome as a string
line = f.readline().rstrip() #start looping from the first line of the sequence or genome
while line: #looping though all the lines in the file
    #the file might contain some empty lines in the end ->stop
    if line!="":
        genome=genome+line #saving the entire genome in 1 line
        line=f.readline().rstrip()
        nlines+=1

if ">" in genome: #meaning that there are more than 1 sequence in the fasta file -> terminate
    print("the script can only process 1 sequence per input file")
    exit()

########################################################
#REVERSE STRAND CREATION
########################################################
rev_gen=reverseStrand(genome)
rev_genome=rev_gen[::-1] #storing the reverse of the compliment of the orginial sequence so it is possible to read from the beginning

########################################################
#IDENTIFYING ORFs
########################################################
ORFs_counter=0 #keepinig track of how many ORFs the user asked to be outputted
fow_ORFs_list=[] #list countaining a seperate dictionary for each forward strand ORF
start=0
#forward strand
while ORFs_counter<int(args.nbF):
    genome_ORF=identifyORF(genome,start,int(args.size))
    fow_ORFs_list.append(genome_ORF)
    ORFs_counter+=1
    start+=1
    if start==3:
            break

#reverse strand
rev_ORFs_list=[] #list countaining a seperate dictionary for each reverse strand ORF
start=0
while ORFs_counter<int(args.nbF):
    rev_genome_ORF=identifyORF(rev_genome,start,int(args.size))
    rev_ORFs_list.append(rev_genome_ORF)
    ORFs_counter+=1
    start+=1

########################################################
#WRITING TO FILE
########################################################
of= open(args.output, "w")
printOutput(args.output,fow_ORFs_list, rev_ORFs_list)
of.close()

##################################################################################################
#END OF PROGRAM
##################################################################################################