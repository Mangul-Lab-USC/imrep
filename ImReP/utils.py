import os
import pysam

beans = "TGACTGTGTTTCTGAACAATAAATGACTTAAACCAGGTATGGCTGCCGATGGTTATCTT"

gencode = {
      'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
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
      'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
      'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W'}


read_map = {"IGH": ["chr14", "14", 105586437, 106879844, 106032614, 107288051, "mouse", "chr12", "12", 113258768, 116009954],
            "IGK": ["chr2", "2", 88857361, 90235368, 89156874, 89630436, "mouse", "chr6", "6", 67555636, 70726754],
            "IGL": ["chr22", "22", 22026076, 22922913, 22380474, 23265085, "mouse", "chr16", "16", 19026858, 19260844],
            "TRA": ["chr14", "14", 21621904, 22552132, 22090057, 23021075, "mouse", "chr14", "14", 52427967, 54224198],
            "TRB": ["chr7", "7", 142299011, 142813287, 141998851, 141998851, "mouse", "chr6", "6", 40891296, 41558371],
            "TRG": ["chr7", "7", 38240024, 38368055, 38279625, 38407656, "mouse", "chr3", "3", 19178042, 19356476]}

basepairs = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

def getOutDir(outputfile):

    outDir = os.path.dirname(outputfile)
    if outDir == "":
        outDir = "."
    if not os.path.exists(outDir):
        os.mkdir(outDir)
    return outDir


def translate_frameshifted(sequence):
    """
    :param sequence:
    :type sequence:
    :return:
    :rtype:
    """
    translate = ''.join([gencode.get(sequence[3*i:3*i+3], 'X')
                        for i in range(len(sequence)//3)])
    return translate

def getKmers(str, k):
    """
    Creates kmers with lenght k from string str
    :param str: string to extract kmers
    :type str: string
    :return: list of kmers with length k generated from the inout string
    :rtype: list
    """
    kmrs = []
    for i in range(len(str) - k + 1):
        kmrs.append(str[i: i + k])
    return kmrs


def reverse_complement(sequence):
    """
    :param sequence:
    :type sequence:
    :return:
    :rtype:
    """
    reversed_sequence = (sequence[::-1])
    rc = ''.join([basepairs.get(reversed_sequence[i], 'X')
                 for i in range(len(sequence))])
    return rc


def nucleotide2protein2(inString):
    """
    :param inString:
    :type inString:
    :return:
    :rtype:
    """
    frames = []
    framesTemp = []

    framesTemp.append(translate_frameshifted(inString[0:]))  # 1st frame
    framesTemp.append(translate_frameshifted(inString[1:]))  # 2nd frame
    framesTemp.append(translate_frameshifted(inString[2:]))  # 3rd frame
    framesTemp.append(translate_frameshifted(
                      reverse_complement(inString)))  # negative 1st frame
    framesTemp.append(translate_frameshifted(
                      reverse_complement(inString)[1:]))  # negative 2nd frame
    framesTemp.append(translate_frameshifted(
                      reverse_complement(inString)[2:]))  # negative 3rd frame

    frs = [1, 2, 3, -1, -2, -3]

    for f, frame in zip(framesTemp, frs):
        if "_" not in f and "*" not in f:
            processedRead = f.replace(' ', '')
            frames.append((processedRead, frame))
    return frames


def dumpClones(clones, outFile):
    """
    :param clones:
    :type clones:
    :param outFile:
    :type outFile:
    :return:
    :rtype:
    """
    with open(outFile, "w") as f:
        for clone in clones:
            f.write((",".join(["%s"] * len(clone)) + "\n") % tuple(clone))


def dumpClones2(clones, outFile):
    """
    :param clones:
    :type clones:
    :param outFile:
    :type outFile:
    :return:
    :rtype:
    """
    header_line = "CDR3_AA_Seq,Chain_type,Read_count,V_chains,D_chains,\
                   J_chains\n"
    with open(outFile, "w") as f:
        f.write(header_line)
        for clone in clones:
            f.write(clone)


def getGeneType(geneName):
    """
    :param geneName:
    :type geneName:
    :return:
    :rtype:
    """
    geneType = geneName.split("|")[1].split("*")[0]
    if "-" in geneType:
        geneType = geneType.split("-")[0]
    return geneType


def getGeneType2(geneName):
    """
    :param geneName:
    :type geneName:
    :return:
    :rtype:
    """
    geneType = geneName.split("|")[1].split("*")[0]
    if "-" in geneType:
        geneType = geneType.split("-")[0]
    return geneType[:4]

def saveMappedReads(samfileName, file, chr, x, y):
    """
    Saves mapped reads from samfile into file object
    """
    k=0
    samfile = pysam.AlignmentFile(samfileName, "rb")

    for read in samfile.fetch(chr, x, y):
        rl = read.infer_query_length()
        c = read.cigartuples

        if c: # to avoid CIGAR of type None
                type = read.cigartuples[0][0]
                length = read.cigartuples[0][1]

                # if read is fully mapped, than we don't take it
                if not (type == int(0) and length == rl):
                    k+=1
                    file.write(">"+str(read.query_name))
                    file.write("\n")
                    file.write(str(read.query_sequence))
                    file.write("\n")
        else:
            k+=1
            file.write(">"+str(read.query_name))
            file.write("\n")
            file.write(str(read.query_sequence))
            file.write("\n")
    samfile.close()
    return k
