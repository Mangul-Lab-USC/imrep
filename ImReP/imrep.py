import sys
import argparse
import os
from collections import Counter
import gzip
import pysam
from Bio import SeqIO
from intervaltree import IntervalTree
import jellyfish

from ImReP.utils import *
from ImReP.settings import *
from ImReP.cast import *
from ImReP import info

if sys.version_info.major == 2:
    str = unicode

try:
    from StringIO import StringIO # Python 2
except ImportError:
    from io import StringIO # Python 3

cd = os.path.dirname(os.path.realpath(__file__))


class ImReP(object):

    def __init__(self):
        self._settings = Settings()
        self.fastq_handle = None
        self.vi_pieces = {}
        self.d_seqs = {}
        self.jay_pieces = {}
        self.pSeq_read_map = {}
        self.just_v = []
        self.just_j = []
        self.just_v_dict = {}
        self.just_j_dict = {}
        self.cdr3_dict = {}
        self.hashV = {}
        self.hashJ = {}
        self.v_chain_type = {}
        self.j_chain_type = {}
        self.debug_info = {}
        self.clonotype_CDR3_count_dict = {}
        self.read_names = {}

    def _bam2fasta(self):
        """
        Creates a fasta file with reads extracted from bam inputfile
        """
        if self._settings.format == "bam":
            print("Parse bam file with mapped and unmapped reads")
            newInput = self._settings.outputDir+"/"+self._settings.sampleName+"_input.fasta"

            with open(newInput, "w") as file:
                for tag in self._settings.chains:
                    k = self._extract_mapped(tag, file)
                    if tag != "TRD":
                        print("Number of reads extacted from ", tag, "locus : ",k)

                self._extract_unmapped(file)

                # sets fasta as new inputfile
                self._settings.inputfile = newInput

    def _extract_mapped(self, tag, file):
        """
        Extracts mapped reads tagged by tag from file and
        saves them into samfile
        """

        if tag not in read_map:
            return 0

        values = read_map.get(tag)
        chr = None
        x = None
        y = None

        if  not self._settings.species:
            if self._settings.is_chrFormat2:
                chr = values[0]
            else:
                chr = values[1]
            if self._settings.is_hg38:
                x = values[2]
                y = values[3]
            else:
                x = values[4]
                y = values[5]
        elif  self._settings.species == values[6]:
            if self._settings.is_chrFormat2:
                chr = values[7]
            else:
                chr = values[8]
            x = values[9]
            y = values[10]

        return saveMappedReads(self._settings.inputfile, file, chr, x, y)

    def _extract_unmapped(self, file):
        k=0
        samfile = pysam.AlignmentFile(self._settings.inputfile, "rb")
        for read in samfile.fetch(until_eof=True):
            if read.is_unmapped:
                k+=1
                file.write(">"+str(read.query_name))
                file.write("\n")
                file.write(str(read.query_sequence))
                file.write("\n")
        print("Number of unmapped reads extracted",k)
        samfile.close()

    def _extract_unmapped_digGold(self, file, mSet):
        """
        Writes unmapped reads extacted from fastqfile settings.allReads
        :param file: output file where reads are written
        :type file: file object
        :param mSet: set of reads
        :type mSet: set
        """
        k=0

        fastq_parser = SeqIO.parse(self.settings.allReads, "fastq")
        for record in fastq_parser:
            if record.id not in mSet:
                file.write(">"+str(record.id))
                file.write("\n")
                file.write(str(record.seq))
                file.write("\n")
                k+=1
        print("Number of unmapped reads extracted", k)

    def _populate_v(self):
        """
        """
        global cd
        chains_v = map(lambda x: cd + "/db/%s/%sV.faa" %
                       (self._settings.species, x), self._settings.chains)
        for ch_v_file in chains_v:
            for record in SeqIO.parse(ch_v_file, "fasta"):
                if "partial in 3'" not in record.description:
                    Vend = str(record.seq)[-20:]
                    kmrs = getKmers(Vend, self._settings.kmer_len)
                    for k in kmrs:
                        if k not in self.hashV:
                            self.hashV[k] = set()
                        self.hashV[k].add(record.id)
                    self.v_chain_type[record.id] = getGeneType2(record.id)
                    posC = Vend.rfind("C")
                    if posC != -1:
                        anchor = Vend[:posC]
                        rest = Vend[posC + 1:]
                        self.vi_pieces[record.id] = (anchor, rest)

    def _populate_d(self):
        """
        """
        global cd
        for chain in self._settings.chains:
            if chain in ["IGH", "TRB", "TRD"]:
                for record in SeqIO.parse(cd + "/db/%s/%sD.faa" %
                    (self._settings.species, chain), "fasta"):

                    if chain not in self.d_seqs:
                        self.d_seqs[chain] = {}
                    self.d_seqs[chain][record.id] = str(record.seq)

    def _populate_j(self):
        """
        """
        global cd
        chains_j = map(lambda x: cd + "/db/%s/%sJ.faa" % (self._settings.species, x), self._settings.chains)
        for ch_j_file in chains_j:
            for record in SeqIO.parse(ch_j_file, "fasta"):
                beginJ = str(record.seq)[:20]
                kmrs = getKmers(beginJ, self._settings.kmer_len)
                for k in kmrs:
                    if k not in self.hashJ:
                        self.hashJ[k] = set()
                    self.hashJ[k].add(record.id)
                self.j_chain_type[record.id] = getGeneType2(record.id)
                letter = "FG"
                if "IGHJ" in ch_j_file:
                    letter = "WG"
                posW = beginJ.find(letter)
                if posW != -1:
                    anchor = beginJ[:posW]
                    rest = beginJ[posW + 1:]
                    self.jay_pieces[record.id] = (anchor, rest)




    def _read_reads(self):
        """
        """

        fastqfile = self._settings.inputfile
        formatFile = "fasta"

        if self._settings.format == "fastq":
            formatFile = "fastq"
        if fastqfile.endswith(".gz"):
            with gzip.open(fastqfile, 'rb') as f:
                file_content = f.readlines()
        else:
            with open(fastqfile) as f:
                file_content = f.readlines()

        # serghei's trick
        a_read_line = file_content[1].strip()
        readlen = len(a_read_line)
        if readlen != 50:
            if self._settings.overlapStep is True:
                self._settings.overlapStep = False

        if formatFile == "fasta":
            while not file_content[0][0] == ">":
                file_content = file_content[1:]
            while not file_content[-2][0] == ">":
                file_content = file_content[:-1]
        elif formatFile == "fastq":
            while not file_content[0][0] == "@":
                file_content = file_content[1:]
            while not file_content[-4][0] == "@" or len(file_content[-1]) != len(file_content[-3]):
                file_content = file_content[:-1]
        else:
            raise Exception("Unrecognized file format: %s!!!" % fastqfile)
        self.fastq_handle = SeqIO.parse(StringIO("".join(file_content)), formatFile)

    def _full_cdr3(self):
        """
        """

        if not self.fastq_handle:
            return []
        vkeys = set(self.hashV.keys())
        jkeys = set(self.hashJ.keys())
        full_cdr3 = []
        for record in self.fastq_handle:
            # If we have paired-end reads,
            # then we have to distinguish them
            if record.id not in self.read_names:
                self.read_names[record.id] = 1
                record.id += "___1"
            else:
                count_existing = self.read_names[record.id]
                record.id += "___%s" % count_existing
                self.read_names[record.id] = count_existing + 1

            self.debug_info[record.id] = {}
            pSequences = nucleotide2protein2(str(record.seq))
            if pSequences:
                for pSeq, frame in pSequences:
                    pos1 = [pSeq.find("C"), pSeq.find("C")]
                    pos2 = [pSeq.rfind("FG"), pSeq.rfind("WG")]
                    v_overlap = "NA"
                    j_overlap = "NA"
                    vtypes = {}
                    jtypes = {}
                    if pos1 != [-1, -1]:
                        if pos1[0] != -1:
                            kmrs1 = getKmers(pSeq[:pos1[0] + 5], self._settings.kmer_len)
                            interV = set(kmrs1) & vkeys
                            vlist = []
                            for v in interV:
                                vlist.extend(list(self.hashV[v]))
                            if vlist:
                                vc = [x for x, y in Counter(vlist).items()]
                            else:
                                vc = []
                            v_cl = {}
                            for v in vc:
                                if self.v_chain_type[v] != "IGHV" and self.v_chain_type[v] not in v_cl:
                                    v_cl[self.v_chain_type[v]] = []
                                if self.v_chain_type[v] != "IGHV":
                                    v_cl[self.v_chain_type[v]].append(v)
                            f, s = pSeq[:pos1[0]], pSeq[pos1[0] + 1:]
                            v_overlap = len(f) + len(s) + 1
                            for v1, v2 in v_cl.items():
                                for v3 in v2:
                                    if v3 not in self.vi_pieces:
                                        continue
                                    v, vv = self.vi_pieces[v3]
                                    minlen1 = min(len(f), len(v))
                                    minlen2 = min(len(s), len(vv))
                                    if minlen1 > 0:
                                        mismatch1 = jellyfish.levenshtein_distance(str(f[-minlen1:]), str(v[-minlen1:]))
                                    else:
                                        mismatch1 = 0
                                    if minlen2 > 0:
                                        mismatch2 = jellyfish.levenshtein_distance(str(s[:minlen2]), str(vv[:minlen2]))
                                    else:
                                        mismatch2 = 0
                                    if (minlen1 <= 3 and mismatch2 <= 1) or (minlen1 >= self._settings.minlen1 and mismatch1 <= self._settings.mismatch1 and minlen2 >= self._settings.minlen2 and mismatch2 <= self._settings.mismatch2):
                                        vtypes[v3] = (minlen1 + minlen2 + 1, mismatch1 + mismatch2)
                        if pos1[1] != -1:
                            kmrs1 = getKmers(pSeq[:pos1[1] + 5], self._settings.kmer_len)
                            interV = set(kmrs1) & vkeys
                            vlist = []
                            for v in interV:
                                vlist.extend(list(self.hashV[v]))
                            if vlist:
                                vc = [x for x, y in Counter(vlist).items()]
                            else:
                                vc = []
                            v_cl = {}
                            for v in vc:
                                if self.v_chain_type[v] == "IGHV" and self.v_chain_type[v] not in v_cl:
                                    v_cl[self.v_chain_type[v]] = []
                                if self.v_chain_type[v] == "IGHV":
                                    v_cl[self.v_chain_type[v]].append(v)

                            f, s = pSeq[:pos1[1]], pSeq[pos1[1] + 1:]

                            v_overlap = len(f) + len(s) + 1
                            for v1, v2 in v_cl.items():
                                for v3 in v2:
                                    if v3 not in self.vi_pieces:
                                        continue
                                    v, vv = self.vi_pieces[v3]
                                    minlen1 = min(len(f), len(v))
                                    minlen2 = min(len(s), len(vv))
                                    if minlen1 > 0:
                                        mismatch1 = jellyfish.levenshtein_distance(str(f[-minlen1:]), str(v[-minlen1:]))
                                    else:
                                        mismatch1 = 0
                                    if minlen2 > 0:
                                        mismatch2 = jellyfish.levenshtein_distance(str(s[:minlen2]), str(vv[:minlen2]))
                                    else:
                                        mismatch2 = 0
                                    if (minlen1 <= 3 and mismatch2 <= 1) or (minlen1 >= self._settings.minlen1 and mismatch1 <= self._settings.mismatch1 and minlen2 >= self._settings.minlen2 and mismatch2 <= self._settings.mismatch2):

                                        vtypes[v3] = (minlen1 + minlen2 + 1, mismatch1 + mismatch2)

                    if pos2 != [-1, -1]:
                        if pos2[0] != -1:
                            if True: #pos2[0] + 3 < len(pSeq) and pSeq[pos2[0] + 3] == "G":
                                if pos2[0] > 10:
                                    offset = pos2[0] - 10
                                else:
                                    offset = 0
                                kmrs2 = getKmers(pSeq[offset:], self._settings.kmer_len)
                                interJ = set(kmrs2) & jkeys
                                jlist = []
                                for j in interJ:
                                    jlist.extend(list(self.hashJ[j]))
                                if jlist:
                                    jc = [x for x, y in Counter(jlist).items()]
                                else:
                                    jc = []
                                j_cl = {}
                                for j in jc:
                                    if self.j_chain_type[j] != "IGHJ" and self.j_chain_type[j] not in j_cl:
                                        j_cl[self.j_chain_type[j]] = []
                                    if self.j_chain_type[j] != "IGHJ":
                                        j_cl[self.j_chain_type[j]].append(j)
                                f, s = pSeq[:pos2[0]], pSeq[pos2[0] + 1:]
                                j_overlap = len(f) + len(s) + 1
                                for j1, j2 in j_cl.items():
                                    for j3 in j2:
                                        if j3 not in self.jay_pieces:
                                            continue
                                        j, jj = self.jay_pieces[j3]
                                        minlen1 = min(len(f), len(j))
                                        minlen2 = min(len(s), len(jj))
                                        if minlen2 > 0:
                                            mismatch2 = jellyfish.levenshtein_distance(str(s[:minlen2]), str(jj[:minlen2]))
                                        else:
                                            mismatch2 = 0
                                        if minlen1 > 0:
                                            mismatch1 = jellyfish.levenshtein_distance(str(f[-minlen1:]), str(j[-minlen1:]))
                                        else:
                                            mismatch1 = 0
                                        if (minlen2 <= 3 and mismatch1 <= 1) or (minlen2 >= self._settings.minlen1 and mismatch2 <= self._settings.mismatch1 and minlen1 >= self._settings.minlen2 and mismatch1 <= self._settings.mismatch2):
                                            jtypes[j3] = (minlen1 + minlen2 + 1, mismatch1 + mismatch2)
                        if pos2[1] != -1:
                            if pos2[1] > 10:
                                offset = pos2[1] - 10
                            else:
                                offset = 0
                            kmrs2 = getKmers(pSeq[offset:], self._settings.kmer_len)
                            interJ = set(kmrs2) & jkeys
                            jlist = []
                            for j in interJ:
                                jlist.extend(list(self.hashJ[j]))
                            if jlist:
                                jc = [x for x, y in Counter(jlist).items()]
                            else:
                                jc = []
                            j_cl = {}
                            for j in jc:
                                if self.j_chain_type[j] == "IGHJ" and self.j_chain_type[j] not in j_cl:
                                    j_cl[self.j_chain_type[j]] = []
                                if self.j_chain_type[j] == "IGHJ":
                                    j_cl[self.j_chain_type[j]].append(j)
                            f, s = pSeq[:pos2[1]], pSeq[pos2[1] + 1:]
                            j_overlap = len(f) + len(s) + 1
                            for j1, j2 in j_cl.items():
                                for j3 in j2:
                                    if j3 not in self.jay_pieces:
                                        continue
                                    j, jj = self.jay_pieces[j3]
                                    minlen1 = min(len(f), len(j))
                                    minlen2 = min(len(s), len(jj))
                                    if minlen2 > 0:
                                        mismatch2 = jellyfish.levenshtein_distance(str(s[:minlen2]), str(jj[:minlen2]))
                                    else:
                                        mismatch2 = 0
                                    if minlen1 > 0:
                                        mismatch1 = jellyfish.levenshtein_distance(str(f[-minlen1:]), str(j[-minlen1:]))
                                    else:
                                        mismatch1 = 0
                                    if (minlen2 <= 3 and mismatch1 <= 1) or (minlen2 >= self._settings.minlen1 and mismatch2 <= self._settings.mismatch1 and minlen1 >= self._settings.minlen2 and mismatch1 <= self._settings.mismatch2):
                                        jtypes[j3] = (minlen1 + minlen2 + 1, mismatch1 + mismatch2)
                    if vtypes or jtypes:
                        vt = {}
                        vscore = {}
                        jt = {}
                        jscore = {}
                        for x in vtypes:
                            chaint = self.v_chain_type[x]
                            if chaint[:3] not in vt:
                                vt[chaint[:3]] = []
                                vscore[chaint[:3]] = []
                            vt[chaint[:3]].append(x)
                            entry = [x] + list(vtypes[x])
                            vscore[chaint[:3]].append(entry)
                        for x in jtypes:
                            chaint = self.j_chain_type[x]
                            if chaint[:3] not in jt:
                                jt[chaint[:3]] = []
                                jscore[chaint[:3]] = []
                            jt[chaint[:3]].append(x)
                            entry = [x] + list(jtypes[x])
                            jscore[chaint[:3]].append(entry)
                        self.debug_info[record.id] = {"vscore": vscore, "jscore": jscore}
                        common = set(vt.keys()) & set(jt.keys())
                        if common:
                            if "IGH" in common:
                                full_cdr3.append(pSeq[pos1[1]: pos2[1] + 1])
                                cdr3 = pSeq[pos1[1]: pos2[1] + 1]
                            else:
                                full_cdr3.append(pSeq[pos1[0]: pos2[0] + 1])
                                cdr3 = pSeq[pos1[0]: pos2[0] + 1]
                            if cdr3 not in self.cdr3_dict:
                                self.cdr3_dict[cdr3] = []
                            self.cdr3_dict[cdr3].append(record.id)
                            if cdr3 not in self.pSeq_read_map or (cdr3 in self.pSeq_read_map and ("v" not in self.pSeq_read_map[cdr3].keys() or "j" not in self.pSeq_read_map[cdr3].keys())):
                                v_t = []
                                j_t = []
                                chtype = {}
                                for key, ch in vt.items():
                                    if key in common:
                                        v_t.extend(ch)
                                        if key not in chtype:
                                            chtype[key] = []
                                        chtype[key].extend(ch)
                                for key, ch in jt.items():
                                    if key in common:
                                        j_t.extend(ch)
                                        if key not in chtype:
                                            chtype[key] = []
                                        chtype[key].extend(ch)
                                self.pSeq_read_map[cdr3] = {"v": map(getGeneType, v_t), "j": map(getGeneType, j_t), "chain_type": chtype}
                        elif vtypes and not jtypes:
                            #if "IGH" in vtypes:
                            #    vi_partial = pSeq[pos1[1]:]
                            #else:
                            #    vi_partial = pSeq[pos1[0]:]
                            vi_partial = pSeq[pos1[1]:]
                            if vi_partial not in full_cdr3:
                                self.just_v.append(vi_partial)
                                if vi_partial not in self.just_v_dict:
                                    self.just_v_dict[vi_partial] = []
                                self.just_v_dict[vi_partial].append(record.id)
                            if vi_partial not in self.pSeq_read_map and vi_partial not in full_cdr3:
                                self.pSeq_read_map[vi_partial] = {"v": map(getGeneType, vtypes), "chain_type": vt}
                        elif jtypes and not vtypes:
                            if "IGH" in jt:
                                jay_partial = pSeq[:pos2[1] + 1]
                            else:
                                jay_partial = pSeq[:pos2[0] + 1]
                            if jay_partial not in full_cdr3:
                                self.just_j.append(jay_partial)
                                if jay_partial not in self.just_j_dict:
                                    self.just_j_dict[jay_partial] = []
                                self.just_j_dict[jay_partial].append(record.id)
                            if jay_partial not in self.pSeq_read_map and jay_partial not in full_cdr3:
                                self.pSeq_read_map[jay_partial] = {"j": map(getGeneType, jtypes), "chain_type": jt}
        return full_cdr3

    def _vj_handshakes(self):
        handshakes = []
        just_v = Counter(self.just_v)
        just_j = Counter(self.just_j)

        itree = IntervalTree()

        just_v_keys = map(lambda x: x[0], sorted(just_v.items(), key=lambda z:z[1], reverse=True))

        start = 0
        for v in just_v_keys:
            end = start + len(v) + 1
            itree.addi(start, end, v)
            start = end

        all_v_suf = "|".join(just_v_keys)
        stree = IgorSuffixTree(all_v_suf)

        for j, jj in just_j.items():
            overlap, index, terminal = stree.search_stree(j)
            if terminal and len(j[:overlap]) >= self._settings.overlapLen:
                overlapping_v = itree.search(index)

                common_chains = set(self.pSeq_read_map[list(overlapping_v)[0].data]["chain_type"].keys()) & set(self.pSeq_read_map[j]["chain_type"].keys())
                if common_chains:
                    v_t = []
                    j_t = []
                    chtype = {}
                    for key, ch in self.pSeq_read_map[list(overlapping_v)[0].data]["chain_type"].items():
                        if key in common_chains:
                            v_t.extend(map(getGeneType, ch))
                            if key not in chtype:
                                chtype[key] = []
                            chtype[key].extend(ch)
                    for key, ch in self.pSeq_read_map[j]["chain_type"].items():
                        if key in common_chains:
                            j_t.extend(map(getGeneType, ch))
                            if key not in chtype:
                                chtype[key] = []
                            chtype[key].extend(ch)
                    if len(j[overlap:]) > 0:
                        newly_born_cdr3 = list(overlapping_v)[0].data + j[overlap:]
                    else:
                        position_of_j_in_v = list(overlapping_v)[0].data.rfind(j)
                        newly_born_cdr3 = list(overlapping_v)[0].data[:position_of_j_in_v + len(j)]
                    if newly_born_cdr3 not in self.cdr3_dict:
                        self.cdr3_dict[newly_born_cdr3] = []
                    if list(overlapping_v)[0].data in self.just_v_dict:
                        self.cdr3_dict[newly_born_cdr3].extend(self.just_v_dict[list(overlapping_v)[0].data])
                    if j in self.just_j_dict:
                        self.cdr3_dict[newly_born_cdr3].extend(self.just_j_dict[j])
                    if list(overlapping_v)[0].data in self.just_v_dict:
                        del self.just_v_dict[list(overlapping_v)[0].data]
                    if j in self.just_j_dict:
                        del self.just_j_dict[j]
                    countV = just_v[list(overlapping_v)[0].data]
                    countJ = just_j[j]
                    countVJ = countV + countJ
                    for x in range(countVJ):
                        handshakes.append(newly_born_cdr3)
                    self.pSeq_read_map[newly_born_cdr3] = {"v": v_t, "j": j_t, "chain_type": chtype, "overlap": overlap}
        return handshakes

    def _map_d(self, seq, chain_type):
        d_types = set()
        for d_t, d_seq in self.d_seqs[chain_type].items():
            if seq.find(d_seq) != -1:
                d_types.add(getGeneType(d_t))
        if not d_types:
            return set(["NA"])
        return sorted(list(d_types))

    def _generateOutputRecord(self,read):

        if "___" in read:
            readId = read.split("___")[0]
        else:
            readId = read
        dinfo_v = self.debug_info[read].get("vscore", {})
        dinfo_j = self.debug_info[read].get("jscore", {})
        di_v = []
        di_j = []
        uniq_v, uniq_j = 0, 0
        for xx, yy in sorted(dinfo_v.items()):
            if yy:
                for u in yy:
                    geneName = u[0].split("|")[1]
                    di_v.append(geneName + ":" +
                                ":".join(map(str, u[1:])))
        for xx, yy in sorted(dinfo_j.items()):
            if yy:
                for u in yy:
                    geneName = u[0].split("|")[1]
                    di_j.append(geneName + ":" +
                                ":".join(map(str, u[1:])))

        if len(di_v) == 1:
            uniq_v = 1
        if len(di_j) == 1:
            uniq_j = 1
        uniq_vj = uniq_v & uniq_j
        di_v = sorted(di_v)
        di_v = ";".join(di_v)
        if not di_v:
            di_v = "NA"
        di_j = sorted(di_j)
        di_j = ";".join(di_j)
        if not di_j:
            di_j = "NA"


        return [readId, di_v, di_j, uniq_v, uniq_j, uniq_vj]

    def _saveOutput(self, clustered_clones):

        final_clones = []
        if self._settings.extendedOutput is True:
            with open(self._settings.outputDir + "/" + "full_cdr3_%s.txt" %self._settings.sampleName, "w") as f:
                header_line = "Read_name,Full_CDR3_AA_Seq,V_genes,D_genes,J_genes,V_allele_name:overlap_aminoacids:mismatches_aminoacids,J_allele_name:overlap_aminoacids:mismatches_aminoacids,Is_V_allele_uniq,Is_V_allele_uniq,Are_both_V_and_J_alleles_uniq\n"
                f.write(header_line)
                for cl in clustered_clones:
                    for clon in self.clone_dict[cl[0]]:
                        for read in self.cdr3_dict[clon]:

                            entry = self._generateOutputRecord(read)

                            f.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % (entry[0], cl[0], cl[1].replace(",",";"), cl[2].replace(",",";"), cl[3].replace(",",";"), entry[1], entry[2], entry[3], entry[4], entry[5]))

                            final_clones.append(cl[0] + ",%s" % cl[1] + ",%s," + "%s,%s,%s\n" % (cl[2].replace(",",";"), cl[3].replace(",",";"), cl[4].replace(",",";")))


            self._savePartialOutput(self._settings.outputDir + "/" + "partial_cdr3_%s.txt" %self._settings.sampleName)
        else:
            for cl in clustered_clones:
                for clon in self.clone_dict[cl[0]]:
                    for read in self.cdr3_dict[clon]:

                        final_clones.append(cl[0] + ",%s" % cl[1] + ",%s," + "%s,%s,%s\n" % (cl[2].replace(",",";"), cl[3].replace(",",";"), cl[4].replace(",",";")))

        final_clones = Counter(final_clones)
        print("%s partial-V CDR3 found" % len(self.just_v_dict))
        print("%s partial-J CDR3 found" % len(self.just_j_dict))
        if len(final_clones):
            print("%s full CDR3 found:" % len(final_clones))
            for x in ['IGH','IGK','IGL','TRA','TRB','TRD','TRG']:
                y = self.clonotype_CDR3_count_dict.get(x, 0)
                print("\t- %s of type %s" % (y, x))
        else:
            print("No full CDR3 found")

        clones = []
        for x, y in final_clones.items():
            clones.append(x % y)

        header_line = "CDR3_AA_Seq,Chain_type,Read_count,V_chains,D_chains,J_chains\n"
        with open(self._settings.outputDir + "/" + "_cdr3_%s.txt" %self._settings.sampleName, "w") as f:
            f.write(header_line)
            for clone in clones:
                f.write(clone)

    def _savePartialOutput(self):


        with open(self._settings.outputDir + "/" + "partial_cdr3_%s.txt" % self._settings.sampleName, "w") as f:
            header_line = "Read_name,Partial_CDR3_AA_Seq,V_genes,D_genes,J_genes,V_allele_name:overlap_aminoacids:mismatches_aminoacids,J_allele_name:overlap_aminoacids:mismatches_aminoacids,Is_V_allele_uniq,Is_V_allele_uniq,Are_both_V_and_J_alleles_uniq\n"
            f.write(header_line)
            for x, y in self.just_v_dict.items():
                for read in y:
                    v = ";".join(sorted(set(self.pSeq_read_map[x].get("v", ["NA"])))[:3])
                    j = ";".join(sorted(set(self.pSeq_read_map[x].get("j", ["NA"])))[:3])
                    entry = self._generateOutputRecord(read)

                    f.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % (entry[0], x, v,
                            "NA", j, entry[1], entry[2], entry[3], entry[4], entry[5]))
            for x, y in self.just_j_dict.items():
                for read in y:
                    v = ";".join(sorted(set(self.pSeq_read_map[x].get("v", ["NA"])))[:3])
                    j = ";".join(sorted(set(self.pSeq_read_map[x].get("j", ["NA"])))[:3])
                    entry = self._generateOutputRecord(read)

                    f.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % (entry[0], x, v,
                            "NA", j, entry[1], entry[2], entry[3], entry[4], entry[5]))

    def info(self):
        release = info.info.get("release", "0.1")
        print("Starting ImReP-%s (developped by %s)" % (release, ", ".join(info.info.get("contributors", ""))))
        print(info.info.get("hello"))

    def computeClones(self, inputfile=None, format="fastq", outputDir=None):

        # storing input information and sample name
        self._settings.sampleName = os.path.splitext(
                                     os.path.basename(inputfile))[0]

        self._settings.inputfile = inputfile

        self._settings.outputDir = outputDir
        if self._settings.outputDir == "":
            self._settings.outputDir = "."
        if not os.path.exists(self._settings.outputDir):
            os.mkdir(self._settings.outputDir)

        self._settings.format = format


        self._settings.checkValues()

        # check file format and converts bam to fasta
        if self._settings.format == "bam":
            self._bam2fasta( )



        if self._settings.is_digGold:
            print("Parse fastq file with orignal raw reads (all reads) \
                   and extract the unmapped reads")

            newInput = self._settings.outputDir+"/"+self._settings.sampleName+"_input.fasta"
            with open(newInput,"w") as file:
               for tag in dict["chains"]:
                   k = 0
                   extract_mapped(tag, file, k)
                   print("Number of reads extacted from ", tag, "locus : ", k)

               samfile = pysam.AlignmentFile(self._settings.inputfile, "rb")
               print  ("Parse bam file with mapped reads")
               mReads = set()
               for read in samfile.fetch():
                   if not read.is_unmapped:
                       mReads.add(read.query_name)
               extract_unmapped_digGold(file, mReads)

               self._settings.inputfile = newInput

        self._populate_v()
        self._populate_d()
        self._populate_j()
        self._read_reads()

        clones = self._full_cdr3()

        # removed for suffixtree implementation
        #if self._settings.overlapStep:
            #clones2 = self._vj_handshakes()
            #clones.extend(clones2)
        clones = Counter(clones)
        for x, y in clones.items():
            if x.endswith("G"): # cleaning of TRA
                del clones[x]
                clones[x[:-1]] = y   ## what if clones[x[:-1]] already exists????
        # here we have to cluster each chain type separately
        clones_by_type = {}

        for cdr3, count in clones.items():
            chtypes = list(map(lambda xx: (xx[0], len(xx[1])),
                           self.pSeq_read_map[cdr3]["chain_type"].items()))
            chtype = [xx for xx, yy in chtypes
                      if yy == max(chtypes, key=lambda zz: zz[1])[1]][0]
            if chtype not in clones_by_type:
                clones_by_type[chtype] = {}
            clones_by_type[chtype][cdr3] = count

        clustered_clones = []
        for chtype, clones in clones_by_type.items():   # suggestion: change clones variable name
            clustered = []

            if self._settings.noCast is True:
                clustered = []
                for clone, count in clones.items():
                    clustered.append([clone, count, [clone]])

            else:
                # execute CAST clustering
                cast_clustering = Cast(clones)
                clustered = cast_clustering.doCast(
                            self._settings.castThreshold[chtype])

            # filter out garbage
            clustered = [cclone for cclone in clustered if cclone[1] >
                         self._settings.filterThreshold]

            for cl in clustered:
                cl.append(chtype)
            self.clonotype_CDR3_count_dict[chtype] = len(clustered)
            clustered_clones.extend(clustered)

        self.clone_dict = {}
        for clone in clustered_clones:

            self.clone_dict[clone[0]] = clone[2]
            chain_type = clone[3]
            del clone[2]
            del clone[1] # remove counts for now
            j_types = None
            if chain_type in ["IGH", "TRB", "TRD"]:
                j_types = self._map_d(clone[0], chain_type)
            types = [",".join(sorted(set(self.pSeq_read_map[clone[0]]["v"]))[:3])]
            if j_types:
                types.append(",".join(j_types))
            else:
                types.append("NA")
            types.append(",".join(sorted(set(self.pSeq_read_map[clone[0]]["j"]))[:3]))
            clone.extend(types)

        self._saveOutput(clustered_clones)

    def setFormat(self, format):
        self._settings.format = format

    def setSpecies(self, species='human'):
        self._settings.species = species

    #def setOverlapStep(self, overlapStep=False):
    #    self._settings.overlapStep = overlapStep

    def setOverlapLen(self, overlapLen=5):
        self._settings.overlapLen = overlapLen

    def setFilterThreshold(self, filterThreshold=1):
        self._settings.filterThreshold = filterThreshold


    def setExtendedOutput(self, extendedOutput=False):
        self._settings.extendedOutput = extendedOutput

    def setCast(self, noCast=True):
        self._settings.noCast = noCast


    def setChains(self, chains = ['IGH','IGK','IGL','TRA','TRB','TRD','TRG']):
        self._settings.chains = chains.split(",")

    def setMinLen(self, minlen1 = 2, minlen2 = 1):
        self._settings.minlen1 = minlen1
        self._settings.minlen2 = minlen2

    def setMismatch(self, mismatch1 = 2, mismatch2 = 2):
        self._settings.mismatch1 = mismatch1
        self._settings.mismatch2 = mismatch2

    def setDigGold(self, is_digGold=False):
        self.is_digGold = is_digGold


    def setChrFormat2(self, is_chrFormat2 = False):
        self._settings.is_chrFormat2 = is_chrFormat2

    def setHg38(self, is_hg38 = False):
        self._settings.is_hg38 = is_hg38

    def setAllReads(self, allReads = None):
        self._settings.allReads = allReads


    def configure(self, args=None):

        if args is not None:
            set_dict={}

            if args.chains:
                set_dict["chains"]=args.chains.split(",")

            if args.species:
                if args.species in ["human", "mouse"]:
                    set_dict["species"] = args.species
                else:
                    raise Exception("Species must be either human or mouse")
            if args.overlapLen:
                set_dict["overlapLen"] = args.overlapLen

            #if args.noOverlapStep is not None:
            #    set_dict["noOverlapStep"] = args.noOverlapStep

            if args.filterThreshold:
                set_dict["filterThreshold"] = args.filterThreshold

            if args.extendedOutput is not None:
                set_dict["extendedOutput"] = args.extendedOutput

            if args.noCast is not None:
                set_dict["noCast"] = args.noCast

            if args.minOverlap1:
                set_dict["minlen1"] = args.minOverlap1

            if args.minOverlap2:
                set_dict["minlen2"] = args.minOverlap2

            if args.misMatch1:
                set_dict["mismatch1"] = args.misMatch1

            if args.misMatch2:
                set_dict["mismatch2"] = args.misMatch2

            if args.is_digGold:
                set_dict["is_digGold"] = args.is_digGold

            if args.allReads:
                set_dict["allReads"] = args.allReads

            if args.is_hg38:
                set_dict["is_hg38"] = args.is_hg38

            if args.is_chrFormat2:
                set_dict["is_chrFormat2"] = args.is_chrFormat2

            if args.isBAM:
                set_dict["format"] = "bam"

            if args.isFastq:
                set_dict["format"] = "fastq"

            self._settings.configure(**set_dict)



if __name__ == "__main__":

    ap = argparse.ArgumentParser("python imrep.py")

    necessary_arguments = ap.add_argument_group("Necessary Inputs")
    necessary_arguments.add_argument("reads_file", help="unmapped reads in .fasta (default) or .fastq (if flag --fastq is set) or .bam (if --bam or --digGold is set)")
    necessary_arguments.add_argument("output_dir", help="output directory to save file with CDR3 clonotypes")

    optional_arguments = ap.add_argument_group("Optional Inputs")
    optional_arguments.add_argument("--fastq", help="a binary flag used to indicate that the input file with unmapped reads is in fastq format", dest="isFastq", action="store_true")
    optional_arguments.add_argument("--bam", help="a binary flag used to indicate that the input file is a BAM file mapped and  unmapped reads", dest="isBAM", action="store_true")
    optional_arguments.add_argument("--chrFormat2", help="a binary flag used to indicate that the format of chromosome name in the bam file is in this format : chr1, chr2,..,chrX. This options is only compatible with --bam option. By default we asssume chromosmes names are indicated only by numbers :1,2,3,...,X", dest="is_chrFormat2", action="store_true")


    optional_arguments.add_argument("--hg38", help="a binary flag used to indicate that reads were mapped to hg38 rellease. The default is hg19. For mouse we support only mm10 (default). ", dest="is_hg38", action="store_true")
    optional_arguments.add_argument("-a", "--allReads", help="Original raw reads (all reads). Needs to be used with --digGold option", type=str, dest="allReads")


    optional_arguments.add_argument("--digGold", help="a binary flag used to indicate that the input file is FASTQ file with original raw reads (all reads). And unmapped reads needs to be extracted from the raw reads ( original raw reads are provided using --reads_file option). Use this option only if unmapped reads were not saved. Needs to be used with -m option", dest="is_digGold", action="store_true")


    optional_arguments.add_argument("-s", "--species", help="species (human or mouse, default human)", type=str, dest="species")
    optional_arguments.add_argument("-o", "--overlapLen", help="the minimal length to consider between reads overlapping with a V gene and reads overlapping with a J gene. Default value is 5 amino acids.", type=int)
    optional_arguments.add_argument("--noOverlapStep", help="a binary flag used in case if the user does not want to run the second stage of the ImReP assembly.", dest="noOverlapStep", action="store_true")
    optional_arguments.add_argument("--extendedOutput", help="extended output: write information read by read", dest="extendedOutput", action="store_true")
    optional_arguments.add_argument("-c", "--chains", help="chains: comma separated values from IGH,IGK,IGL,TRA,TRB,TRD,TRG", type=str)
    optional_arguments.add_argument("--noCast", help="specify this option if you want to disable CDR3 clustering", dest="noCast", action="store_true")
    optional_arguments.add_argument("-f", "--filterThreshold", help="filter out clonotypes with readcount less or equal than filterThreshold (remove outliers), default is 1", type=int)

    advanced_arguments = ap.add_argument_group("Advanced Inputs")
    advanced_arguments.add_argument("--minOverlap1", help="minimal overlap between the reads and A) the left part of V gene (before C amino acid) and B) the right part of J gene (after W for IGH and F for all other chains), default is 4", type=int)
    advanced_arguments.add_argument("--minOverlap2", help="minimal overlap between the reads and A) the right part of V gene (after C amino acid) and B) the left part of J gene (before W for IGH and F for all other chains), default is 1", type=int)
    advanced_arguments.add_argument("--misMatch1", help="maximal number of mismatches between the reads and A) the left part of V gene (before C amino acid) and B) the right part of J gene (after W for IGH and F for all other chains), default is 2", type=int)
    advanced_arguments.add_argument("--misMatch2", help="maximal number of mismatches between the reads and A) the right part of V gene (after C amino acid) and B) the left part of J gene (before W for IGH and F for all other chains), default is 2", type=int)


    args = ap.parse_args()

    imrep = ImReP()

    imrep.configure(args)

    imrep.computeClones(inputfile=args.reads_file, format=args.format, outputDir=args.output_dir)
