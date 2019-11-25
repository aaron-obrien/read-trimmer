import string
import collections
import subprocess
import itertools
import re
import edlib
import os
import pickle
import multiprocessing

import pyximport; pyximport.install(language_level=3)
from _utils import quality_trim

class PrimerDataStruct(object):
    '''
    '''
    def __init__(self,*args,**kwargs):
        '''
        '''
        self.primer_file = kwargs["primer_file"]
        self.k = kwargs["k"]
        self.ncpu = kwargs["ncpu"]
        self.seqtype = kwargs["seqtype"]
        self.primer_col = kwargs["primer_col"]

        self.save_cache = False if not "save_cache" in kwargs else self.save_cache
        self.load_cache = False if not "load_cache" in kwargs else self.load_cache

    def _pairwise_align(self,idx):
        ''' Pairwise align two primers
        :param tuple idx: (i,j) , the indices of the two primers to align
        :rtype tuple
        :returns The shorter and longer primer, (-1,-1) if no alignment
        '''
        i,j = idx
        if i == j:
            return(-1,-1)
        shorter_primer,longer_primer = (primers[i],primers[j]) \
                                       if len(primers[i]) < len(primers[j]) \
                                          else (primers[j],primers[i])
        editdist_cutoff = int(0.1*len(shorter_primer))
        temp = edlib.align(shorter_primer,longer_primer,mode="HW",k=editdist_cutoff) # align in infix mode,
        if temp["editDistance"] != -1:
            return (shorter_primer,longer_primer)
        else:
            return (-1,-1)
        
        
    def _cluster_primer_seqs(self):
        ''' Pairwise comparison of primers , store information
        for other closley matching primers for each primer
        '''
        global primers
        primers = []
        self._primer_clusters = collections.defaultdict(set)
        with open(self.primer_file,"r") as IN:
            for line in IN:
                if line.startswith("#"):
                    continue
                primer = line.strip("\n\r").split("\t")[self.primer_col]
                if primer not in primers:
                    primers.append(primer)

        pairwise_indices = ((i,j) for j in range(0,len(primers)) \
                            for i in range(0,len(primers)))

        p = multiprocessing.Pool(self.ncpu)
        chunksize = int(len(primers)/self.ncpu)
        if chunksize < 1:
            chunksize = 2

        for res in p.imap(self._pairwise_align,pairwise_indices,chunksize=chunksize):
            if res[0] == -1:
                continue
            shorter_primer,longer_primer = res
            self._primer_clusters[shorter_primer].add(longer_primer)
            self._primer_clusters[longer_primer].add(shorter_primer)
        p.close()
        p.join()           
            
    def _create_primer_search_datastruct(self):
        ''' Create k-mer index on the primers
        '''
        self._primer_kmers = collections.defaultdict(set)
        self._primer_info = collections.defaultdict(list) 

        primer_id=0
        with open(self.primer_file,"r") as IN:
            for line in IN:
                if line.startswith("#"):
                    continue
                contents = line.strip("\n\r").split("\t")
                primer   = contents[self.primer_col]
                if primer in self._primer_info:
                    temp = [str(primer_id),contents,[]]
                    self._primer_info[primer].append(temp) # store exact matching primers in the same value bucket
                else:
                    temp = [str(primer_id),contents,[]] # empty list for clustered primers from cd-hit
                    self._primer_info[primer] = [temp]

                    # create k-mer index
                    k = 8
                    kmers = set("".join(itertools.islice(primer,i,i+k)) for i in range(len(primer)+1-k))
                    for kmer in kmers:
                        self._primer_kmers[kmer].add(primer)

                primer_id+=1

        # add close primer sequences for each primer based on pairwise clustering
        for p1 in self._primer_clusters:
            for p2 in self._primer_clusters[p1]:
                self._primer_info[p1][0][2].append(p2)

    @property
    def primer_search_datastruct(self):
        '''
        '''
        if self.load_cache:
            if os.path.exists(self.primer_file+".index.cache"):
                with open(self.primer_file+".index.cache","rb") as cache:
                    return pickle.load(cache)
        self._cluster_primer_seqs()
        self._create_primer_search_datastruct()
        if self.save_cache:
            with open(self.primer_file+".index.cache.temp","wb") as cache:
                pickle.dump((self._primer_info,self._primer_kmers),cache)
            
        return (self._primer_info,self._primer_kmers)

class Trimmer(object):
    '''
    '''
    def __init__(self,*args,**kwargs):
        ''' Class constructor
        '''        
        # user must provide
        self.is_nextseq                = kwargs["is_nextseq"]
        self.max_mismatch_rate_primer  = kwargs["max_mismatch_rate_primer"]
        self.max_mismatch_rate_overlap = kwargs["max_mismatch_rate_overlap"]
        self.umi_len                   = kwargs["umi_len"]
        self.umi_len_alt               = kwargs["umi_len_alt"]
        self.common_seq_len            = kwargs["common_seq_len"]
        self.check_primer_side         = kwargs["check_primer_side"]
        self.overlap_check_len         = kwargs["overlap_check_len"]        
        self.min_primer_side_len       = kwargs["min_primer_side_len"]
        self.min_umi_side_len          = kwargs["min_umi_side_len"]
        self.umi_filter_min_bq         = kwargs["umi_filter_min_bq"]
        self.umi_filter_max_lowQ_bases = kwargs["umi_filter_max_lowQ_bases"]
        self.umi_filter_max_Ns         = kwargs["umi_filter_max_Ns"]
        self.primer3_R1                = kwargs["primer3_R1"]
        self.primer3_R2                = kwargs["primer3_R2"]
        self.tagname_umi               = kwargs["tagname_umi"]
        self.tagname_primer            = kwargs["tagname_primer"]
        self.tagname_primer_error      = kwargs["tagname_primer_error"]
        self.tag_separator             = kwargs["tag_separator"]
        self.field_separator           = kwargs["field_separator"]
        self.no_tagnames               = kwargs["no_tagnames"]
        self.trim_custom_seq_adapter   = kwargs["trim_custom_seq_adapter"]
        self.custom_seq_adapter        = kwargs["custom_seq_adapter"]
        self.poly_tail_primer_side     = kwargs["poly_tail_primer_side"]
        self.poly_tail_umi_side        = kwargs["poly_tail_umi_side"]
        self.drop_alt_seqtype          = kwargs["drop_alt_seqtype"]
        self.include_common_seq_tag    = kwargs["include_common_seq_tag"]
        
        # user can overide these defaults if needed
        self._k = 8  # kmer size
        self._r = 30 # how much of the read to use to build kmers

        # other constants
        self._revcomp_table             = bytes.maketrans(b"ACTG", b"TGAC")
        self._padding                   = 5
        self._poly_type_motif = {
            "polyA"        : re.compile(b"^([ACGTN]*?[CGTN])([A]{8,}[ACGNT]*$)"),
            "polyT"        : re.compile(b"^([ACGTN]*?[CGAN])([T]{8,}[ACGNT]*$)"),
            "polyT_5prime" : re.compile(b"^([T]{5,})([ACGTN]*$)")
        }

        # now possible to specify comma delimeted list of custom sequencing adapters
        self.custom_seq_adapter = self.custom_seq_adapter.replace(b" ", b"").split(b",")

    # kmer size
    @property
    def k(self):
        return self._k
    @k.setter
    def k(self,val):
        assert isinstance(val,int),"Please specify an integer value for k"
        self._k = val
        
    # how much of read sequence to use
    @property
    def r(self):
        return self._r        
    @r.setter
    def r(self,val):
        assert isinstance(val,int),"Please specify an integer value for k"
        self._r = val

    '''
    Assumption made is that R1 is the Primer side , R2 is the UMI side.
    If the structure is the opposite , one can do so by flipping R1 and R2 setters
    '''
    # R1 info
    @property
    def r1_info(self):
        return self._r1_info
    @r1_info.setter
    def r1_info(self,info_tuple):
        self._r1_info = info_tuple
    # R2 info
    @property
    def r2_info(self):
        return self._r2_info
    @r2_info.setter
    def r2_info(self,info_tuple):
        self._r2_info = info_tuple

    # some checks on fastq format
    def check_fastq_format(self):
        ''' Check to make sure reads adhere to FastQ format
        '''
        # unpack R1,R2 info
        r1_id,r1_seq,r1_qual = self._r1_info
        r2_id,r2_seq,r2_qual = self._r2_info

        if not r1_id.split(b" ")[0] == r2_id.split(b" ")[0]:
            raise UserWarning("R1: {} R2: {} read ids are different".format(
                r1_id.decode("ascii"), r2_id.decode("ascii")))

        if not r1_id.startswith(b"@"):
            raise UserWarning("Read ID : {} does not startwith @".format(
                r1_id.decode("ascii")))

        if not len(r1_qual) == len(r1_seq):
            raise UserWarning("R1 Qual: {} and Seq: {} lengths differ".format(
                r1_qual.decode("ascii"), r1_seq.decode("ascii")))

        if not len(r2_qual) == len(r2_seq):
            raise UserWarning("R2 Qual: {} and Seq: {} lengths differ".format(
                r2_qual.decode("ascii"), r2_seq.decode("ascii")))
        
        if len(r2_seq) == 0:
            raise UserWarning("Empty Fastq sequence for Read : {}".format(
                r1_id.decode("ascii")))
        
    def custom_sequencing_adapter_check(self,r1_seq):
        ''' Check for custom sequencing adapter on r1
        :param bytes r1_seq: R1 sequence
        :rtype int
        :returns end pos of adapter, -1 if not found
        '''
        # NOTE: Will check multiple adapters for all reads
        # Should the adapter be identified first and only that sequence used going forward ?
        # Maybe add this logic in the trim_custom_sequencing_adapter function.
        best_adapter = {"seq":None, "align":None, "score":None}
        for a in self.custom_seq_adapter:
            alignment = edlib.align(a, r1_seq[0:len(self.custom_seq_adapter)+3],
                                    mode="SHW",task="locations")
            score = float(alignment["editDistance"])/len(self.custom_seq_adapter)
            if best_adapter["seq"] is None or score < best_adapter["score"]:
                best_adapter["seq"]   = a
                best_adapter["align"] = alignment
                best_adapter["score"] = score

        if best_adapter["score"] < 0.18:
            return best_adapter["align"]["locations"][-1][1]
        else:
            return -1

    def primer_trim(self,primer_datastruct,r1_seq):
        ''' Trim primer to said bases if present on read, account for certain mismatch
        :param tuple primer_datastruct: (primer_info,primer_kmer)
        :param str   r1_seq: r1 read sequence
        '''
        # initialize variables
        best_editdist       = None
        best_score          = None
        best_primer         = None
        best_primer_len     = None
        candidates          = set()
        
        # unpack argument
        primer_info, primer_kmer = primer_datastruct
        
        temp = r1_seq[0:self._r].decode("ascii")
        for oligo in set("".join(itertools.islice(temp,i,i+self._k))
                         for i in range(self._r+1-self._k)):
            if oligo in primer_kmer:
                for c in primer_kmer[oligo]:
                    candidates.add(c)
                    for similar_primer in primer_info[c][0][2]: # iterate over all similar primers to this primer
                        candidates.add(similar_primer)
        
        if len(candidates) == 0: # no hits in the index, exhaustive search over all primers
            candidates = primer_info
        num_candidates = len(candidates)
        
        for primer in candidates:
            primer_len = len(primer)
            alignment = edlib.align(primer,r1_seq[0:len(primer)+self._padding],mode="SHW")
            editdist  = alignment["editDistance"]
            score     = float(editdist)/primer_len

            if best_score is None or score < best_score:
                best_align      = alignment
                best_primer     = primer
                best_score      = score
                best_editdist   = editdist
                best_primer_len = primer_len
            elif score == best_score:
                if best_primer_len < primer_len:
                    best_align      = alignment
                    best_primer     = primer
                    best_score      = score
                    best_editdist   = editdist
                    best_primer_len = primer_len

        assert best_score is not None

        if best_score <= self.max_mismatch_rate_primer:
            return (best_align["locations"][-1][1], best_primer, best_editdist, best_score) # return 0 based position where primer ends on the read
        else:
            return (-1,None,-1,-1)
        

    '''                   Overlap checks

         Primer  
         ======--------------------------=====>        R1
               xxxxxx
                                   xxxxxx 
           <===--------------------------========      R2
                                      UMI (synthetic side)

           ===   : the synthetic part we want to trim
           xxx   : the regions checked for overlap of R1 and R2
           ---   : endogenous read sequence

    '''    
        
    def synthetic_side_check(self,r1_seq,r2_seq):
        ''' Check for overlap between R1 and R2 on the synthetic/UMI side
        :param bytes/str r1_seq: R1 read sequence
        :param bytes/str r2_seq: R2 read sequence
        :rtype tuple
        :returns (start,end) pos on R1 (0-based coordinates) 
                 (-1,-1) if no alignment within mismatch rate
        '''  
        query = self.revcomp(r2_seq[self.synthetic_oligo_len:self.synthetic_oligo_len+self.overlap_check_len])
        if len(query) == 0:
            return (-1,-1)
        alignment = edlib.align(query,r1_seq,mode="HW",task="locations")
        if float(alignment["editDistance"])/self.overlap_check_len <= self.max_mismatch_rate_overlap:
            return alignment["locations"][-1]
        else:
            return (-1,-1)
        
    def primer_side_check(self,primer_end,r1_seq,r2_seq):
        '''
        '''
        query = self.revcomp(r1_seq[primer_end+1:primer_end+1+self.overlap_check_len])
        if len(query) == 0:
            return (-1,-1)
        alignment = edlib.align(query,r2_seq,mode="HW",task="locations")
        if float(alignment["editDistance"])/self.overlap_check_len <= self.max_mismatch_rate_overlap:
            return alignment["locations"][-1]
        else:
            return (-1,-1)        
        
    def quality_trim_(self,qual_string,seq_string,cutoff,base=33):
        return quality_trim(qual_string,seq_string,cutoff,self.is_nextseq,base)

    def revcomp(self,seq):
        ''' Reverse complement a sequence
        '''
        return seq.translate(self._revcomp_table)[::-1]

    def poly_trim(self,seq, poly_type):
        ''' Return start pos of polyA/T tail/5'end , if present
        returns -1 if not found

        :param bytes seq: The read sequence to trim
        :rtype int
        :returns the polyA start position
        '''
        match = self._poly_type_motif[poly_type].match(seq)
        if match: # found polyA
            return match.start(2)
        else:
            return -1
