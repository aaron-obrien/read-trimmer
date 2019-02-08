import multiprocessing
import gzip
import functools
import logging
import sys

import edlib

from pprint import pprint

from trimmer import PrimerDataStruct, Trimmer
from _utils import two_fastq_heads


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler(sys.stdout)
logger.addHandler(ch)

class QiaSeqTrimmer(Trimmer):
    '''
    '''
    def __init__(self,*args,**kwargs):
        super(QiaSeqTrimmer,self).__init__(*args,**kwargs)
        
        # r1,r2 info 
        self._r1_info = None
        self._r2_info = None

        self.seqtype = kwargs["seqtype"]
        
        # some duplex specific params
        self.tagname_duplex = kwargs["tagname_duplex"]
        self.is_duplex = kwargs["is_duplex"]
        if self.is_nextseq:
            self._duplex_adapters = [b"TTCTGAGCGAYYATAGGAGTCCT",b"GTTCTGAGCGAYYATAGGAGTCCT",
                                     b"CGTTCTGAGCGAYYATAGGAGTCCT",b"ACGTTCTGAGCGAYYATAGGAGTCCT"]
        else:
            self._duplex_adapters = [b"TTCTGAGCGAYYATAGGAGTCCT"]
            
        # some multimodal specific params
        self.tagname_multimodal = kwargs["tagname_multimodal"]
        self.is_multimodal = kwargs["is_multimodal"]
        self._multimodal_UMIend_adapters = [(b"ATTGGAGTCCT", "dna", b"B"), (b"ACGTTTTTTTTTTTTTTTTTTVN", "rna", b"RT"), (b"ATCTGCGGG", "rna", b"TSO")]
        if self.umi_len_alt is None:
            self.umi_len_alt = self.umi_len

        # some duplex or multimodal specific params
        self.synthetic_oligo_len = None if self.is_duplex or self.is_multimodal else self.umi_len + self.common_seq_len

    # boolean variables requested from this class
    @property
    def is_r1_primer_trimmed(self):
        return self._is_r1_primer_trimmed
    @property
    def is_r1_syn_trimmed(self):
        return self._is_r1_syn_trimmed
    @property
    def is_r2_primer_trimmed(self):
        return self._is_r2_primer_trimmed
    @property
    def is_r1_r2_overlap(self):
        return self._is_r1_r2_overlap
    @property
    def is_too_short(self):
        return self._is_too_short
    @property
    def is_odd_structure(self):
        return self._is_odd_structure
    @property
    def is_bad_umi(self):
        return self._is_bad_umi
    @property
    def is_r1_qual_trim(self):
        return self._is_r1_qual_trim
    @property
    def is_r2_qual_trim(self):
        return self._is_r2_qual_trim
    @property
    def is_r1_poly_tail_trim(self):
        return self._is_r1_poly_tail_trim
    @property
    def is_r2_poly_tail_trim(self):
        return self._is_r2_poly_tail_trim
    

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
        
        
    # duplex specific vars
    @property
    def is_duplex_adapter_present(self):
        return self._is_duplex_adapter_present
    @property
    def duplex_tag(self):
        return self._duplex_tag

    def _id_duplex_tag(self,r2_seq):
        ''' Identify duplex tag
        :param bytes r2_seq: R2 sequence

        :rtype tuple
        :returns (adapter_len -1,duplex_tag)
        '''
        best_editdist = None        
        start_pos = self.umi_len + 1 if r2_seq[0] == b"N" else self.umi_len
        for a in self._duplex_adapters:
            subseq_r2 = r2_seq[start_pos:start_pos+len(a)+3]            
            align = edlib.align(a,subseq_r2,mode="SHW",
                                   additionalEqualities=[("Y", "C"), ("Y", "T")]) # prefix mode
            editdist = align["editDistance"]
            score    = float(editdist)/len(a)
            if editdist == 0:
                best_alignment = align
                best_editdist  = editdist
                best_score     = score
                best_adapter   = a
                best_subseq    = subseq_r2
                break
            elif best_editdist == None or score < best_score:
                best_alignment = align
                best_editdist  = editdist
                best_score     = score
                best_adapter   = a
                best_subseq    = subseq_r2
        if best_score <= 0.18:
            end_pos = best_alignment["locations"][-1][1] # multiple alignments with same scores can be there , return end pos of last alignment(furthest end pos)
            duplex_tag = b"CC" if best_subseq[end_pos - 13:end_pos - 9].find(b"CC") != -1 \
                    else b"TT" if best_subseq[end_pos - 13:end_pos - 9].find(b"TT") != -1 \
                    else b"NN"
            return (end_pos,duplex_tag)
        else:
            return (-1,"-1")
        
    # multimodal specific vars
    @property
    def is_multimodal_adapter_present(self):
        return self._is_multimodal_adapter_present
    @property
    def multimodal_adapter_name(slef):
        return self._multimodal_adapter_name
        
    def _id_multimodal_adapter_name(self,r2_seq):
        ''' Identify multimodal adapter sequence on UMI end
        :param bytes r2_seq: R2 sequence

        :rtype tuple
        :returns (adapter_len -1,multimodal_adapter_name,umi_len)
        '''
        best_editdist = None
        for a, b, c in self._multimodal_UMIend_adapters:            
            umi_len = self.umi_len if b == self.seqtype else self.umi_len_alt
            start_pos = umi_len + 1 if r2_seq[0] == b"N" else umi_len
            subseq_r2 = r2_seq[start_pos:start_pos+len(a)+3]            
            align = edlib.align(a,subseq_r2,mode="SHW",additionalEqualities=[("V", "A"), ("V", "C"), ("V", "G")]) # prefix mode
            editdist = align["editDistance"]
            score    = float(editdist)/len(a)
            if editdist == 0:
                best_alignment = align
                best_editdist  = editdist
                best_score     = score
                best_adapter   = a
                best_subseq    = subseq_r2
                best_umi_len   = umi_len
                best_adpt_name = c
                break
            elif best_editdist == None or score < best_score:
                best_alignment = align
                best_editdist  = editdist
                best_score     = score
                best_adapter   = a
                best_subseq    = subseq_r2
                best_umi_len   = umi_len
                best_adpt_name = c
        if best_score <= 0.2:
            end_pos = best_alignment["locations"][-1][1] # multiple alignments with same scores can be there , return end pos of last alignment(furthest end pos)                        
            return (end_pos,best_adpt_name,best_umi_len)
        else:
            return (-1,"-1",self.umi_len)

    def _reformat_readid(self,read_id,umi,primer_id,primer_error):
        ''' update read id with umi and primer info
        :param read_id: bytes
        :param umi: bytes
        :param primer_id: bytes
        :param primer_error: bytes
        '''
        idx          = read_id.find(b" ") if read_id.find(b" ") != -1 else len(read_id)
        primer_id    = primer_id.encode("ascii")
        primer_error = primer_error.encode("ascii")

        if self.no_tagnames:
            umi_info = umi
            primer_info = primer_id
            primer_error_info = primer_error
        else:
            umi_info    = b":".join([self.tagname_umi,b"Z",umi])
            primer_info = b":".join([self.tagname_primer,b"Z",primer_id])
            primer_error_info = b":".join([self.tagname_primer_error,b"Z",primer_error])
        
        read_id_info = [read_id[0:idx],umi_info,primer_info,primer_error_info]
        
        if self.is_duplex:
            if self.is_multimodal:
                raise Exception("Cannot handle with duplex tags from multimodal reads!")
            elif self.no_tagnames:
                duplex_info = self.duplex_tag
            else:
                duplex_info = b":".join([self.tagname_duplex,b"Z",self.duplex_tag])
            read_id_info.append(duplex_info)
        
        if self.is_multimodal:
            if self.no_tagnames:
                multimodal_adaper_info = self.multimodal_adapter_name
            else:
                multimodal_adaper_info = b":".join([self.tagname_multimodal,b"Z",self.multimodal_adapter_name])
            read_id_info.append(multimodal_adaper_info)
        
        return self.tag_separator.join(read_id_info)

    def _umi_filter_rna(self,umi,umi_qual):
        ''' filter umi based on sequence and base qualities        
        '''
        num_Ns = umi.count(b"N")
        if num_Ns > self.umi_filter_max_Ns:
            return True
        counter = 0
        base = 33
        for q in umi_qual.decode("ascii"):
            qual = ord(q) - base
            if qual < self.umi_filter_min_bq:
                counter+=1
                
        if counter < self.umi_filter_max_lowQ_bases:
            return False
        
    def qiaseq_trim(self,primer_datastruct):
        """ Trim QiaSeq DNA/RNA reads
        """
        # init some bools
        self._is_r1_primer_trimmed      = False
        self._is_r1_syn_trimmed         = False
        self._is_r2_primer_trimmed      = False
        self._is_r1_r2_overlap          = False
        self._is_too_short              = False
        self._is_odd_structure          = False
        self._is_bad_umi                = False # only used for RNA

        self._is_r1_qual_trim           = False
        self._is_r2_qual_trim           = False
        
        self._is_r1_poly_tail_trim      = False # only used for RNA
        self._is_r2_poly_tail_trim      = False # only used for RNA
        
        self._is_duplex_adapter_present = False
        self._duplex_tag                = None
        self._multimodal_adapter_name   = None
        
        # init local variables
        primer_side_overlap = False
        syn_side_overlap    = False
        primer_id = "-1"
        primer_error = "-1"
        umi_len = self.umi_len # for multimodal reads, UMI length depends on sequence type (in case of leakage)
        
        # unpack R1,R2 info
        r1_id,r1_seq,r1_qual = self._r1_info
        r2_id,r2_seq,r2_qual = self._r2_info

        # trim custom sequencing adapter if needed
        if self.trim_custom_seq_adapter:
            r1_trim_start = self.custom_sequencing_adapter_check(r1_seq)
            # update r1
            if r1_trim_start != -1:
                r1_seq = r1_seq[r1_trim_start+1:]
                r1_qual = r1_qual[r1_trim_start+1:]

        # identify duplex adapter and tag if needed
        if self.is_duplex:
            adapter_end_pos,self._duplex_tag = self._id_duplex_tag(r2_seq)
            if adapter_end_pos == -1: # drop read                
                return            
            # update synthetic oligo len
            self._is_duplex_adapter_present = True
            self.synthetic_oligo_len = umi_len + (adapter_end_pos + 1)
            
        # identify multimodal UMI end adapter, if needed
        if self.is_multimodal:
            adapter_end_pos,self._multimodal_adapter_name,umi_len = self._id_multimodal_adapter_name(r2_seq)
            if adapter_end_pos == -1: # drop read                
                return            
            # update synthetic oligo len
            self._is_multimodal_adapter_present = True
            self.synthetic_oligo_len = umi_len + (adapter_end_pos + 1)
            
        # get umi
        synthetic_oligo_len = self.synthetic_oligo_len        
        if not r2_seq.startswith(b"N"):
            umi = r2_seq[0:umi_len]
            umi_qual = r2_qual[0:umi_len]
        else:
            umi = r2_seq[1:umi_len+1]
            umi_qual = r2_qual[1:umi_len+1]
            self.synthetic_oligo_len += 1
            
        if self.seqtype == "rna":
            self._is_bad_umi = self._umi_filter_rna(umi,umi_qual)            
            if self._is_bad_umi: # drop read
                self.synthetic_oligo_len = synthetic_oligo_len
                return
        
        # quality trimming        
        r1_qual_end = self.quality_trim_(r1_qual,r1_seq,14)
        r2_qual_end = self.quality_trim_(r2_qual,r2_seq,14)
        temp = len(r1_seq)  - r1_qual_end
        self.r1_qual_trim_len = temp
        if temp > 0:
            self._is_r1_qual_trim = True
        temp = len(r2_seq)  - r2_qual_end
        if temp > 0:
            self._is_r2_qual_trim = True
        self.r2_qual_trim_len = temp
        
        # update r1,r2
        r1_seq  = r1_seq[0:r1_qual_end]
        r1_qual = r1_qual[0:r1_qual_end]
        r2_seq  = r2_seq[0:r2_qual_end]
        r2_qual = r2_qual[0:r2_qual_end]        
                    
        if len(r1_seq) < self.min_primer_side_len or len(r2_seq) < self.min_umi_side_len: # skip reads too short after qual trimming
            self._is_too_short = True
            self.synthetic_oligo_len = synthetic_oligo_len
            return
            
        # track qual trimmed lengths
        r1_len = len(r1_seq)
        r2_len = len(r2_seq)

        # trim primer on R1
        r1_primer_end_pos,primer,editdist,score = self.primer_trim(primer_datastruct,r1_seq)

        if r1_primer_end_pos == -1:
            r1_id = self._reformat_readid(r1_id,umi,primer_id,primer_error)
            r2_id = self._reformat_readid(r2_id,umi,primer_id,primer_error)
            # just trim synthetic part on R2
            r2_seq = r2_seq[self.synthetic_oligo_len:]
            r2_qual = r2_qual[self.synthetic_oligo_len:]
            # update read info tuple
            self._r1_info = (r1_id,r1_seq,r1_qual)
            self._r2_info = (r2_id,r2_seq,r2_qual)
            # reset var
            self.synthetic_oligo_len = synthetic_oligo_len

            # check r2 length
            if len(r2_seq) < self.min_umi_side_len:
                self._is_too_short = True
                return
            if len(r2_seq) == 0: # min_umi_side_len is 0
                # add dummy sequence so as not to fail downstream
                r2_seq  = b"N"
                r2_qual = b"!"
                self._r1_info = (r1_id,r1_seq,r1_qual)
                self._r2_info = (r2_id,r2_seq,r2_qual)
                return
            
            return
        
        # set annotation for primer tag
        temp = primer_datastruct[0][primer][0]
        if self.seqtype == "dna":
            chrom,pos,strand,seq = temp[1]
            primer_info = chrom+"-"+strand+"-"+pos
        elif self.seqtype == "rna":
            pr_id = []
            # could have 2 or more primers with same sequence
            # use primer id in info field
            for p in primer_datastruct[0][primer]:
                pr_id.append(p[0])
            primer_info = ",".join(pr_id)
            
        primer_error = str(editdist)

        # look for overlap on synthetic side , i.e. align endo seq from R2 on R1
        syn_side_overlap_start,syn_side_overlap_end = self.synthetic_side_check(r1_seq,r2_seq)
        if syn_side_overlap_start != -1: # synthetic side overlap was successful
            syn_side_overlap = True

        if self.check_primer_side or not syn_side_overlap:
            # look for overlap on primer side , i.e. align endo seq from R1 on R2
            primer_side_overlap_start,primer_side_overlap_end = self.primer_side_check(r1_primer_end_pos,r1_seq,r2_seq)
            if primer_side_overlap_start != -1:
                primer_side_overlap = True

        num_primer_bases_R2 = r1_primer_end_pos + 1 if self.primer3_R2 == -1 else self.primer3_R2 # keep all primer bases if primer3_R2 is -1
        if syn_side_overlap:
            r1_trim_end = syn_side_overlap_end + 1
            if self.check_primer_side and primer_side_overlap : # use primer side coordinates for trimming
                if self.primer3_R2 == -1: # no primer bases to trim
                    r2_trim_end = primer_side_overlap_end + 1 + num_primer_bases_R2
                else:
                    r2_trim_end = primer_side_overlap_end + 1 + num_primer_bases_R2 # python will truncate to string end pos if we overflow the length
            else: # use syn side coordinates
                r2_trim_end = self.synthetic_oligo_len + (syn_side_overlap_end - r1_primer_end_pos) + num_primer_bases_R2
        else:
            if primer_side_overlap:
                # use primer side coordinates for trimming            
                r2_trim_end = primer_side_overlap_end + 1 + num_primer_bases_R2
                r1_trim_end = r1_primer_end_pos + 1 + primer_side_overlap_end + 1 - self.synthetic_oligo_len
            else: # no overlap on syn or primer side
                r1_trim_end = r1_len
                r2_trim_end = r2_len
                
        # update r1,r2 qual and sequence
        if self.primer3_R1 == -1: # no primer bases to trim
            r1_trim_start = 0
        else:
            r1_trim_start = r1_primer_end_pos + 1 - self.primer3_R1
            
        r2_trim_start = self.synthetic_oligo_len

        if r1_trim_start >= r1_trim_end: # weird reads , mostly only primer on R1 , no endo seq
            self._is_odd_structure = True
            self.synthetic_oligo_len = synthetic_oligo_len
            return            
        
        r1_seq  = r1_seq[r1_trim_start:r1_trim_end]
        r1_qual = r1_qual[r1_trim_start:r1_trim_end]
        r2_seq  = r2_seq[r2_trim_start:r2_trim_end]
        r2_qual = r2_qual[r2_trim_start:r2_trim_end]

        # update bools
        self._is_r1_primer_trimmed = True
        if r1_trim_end < r1_len:
            self._is_r1_syn_trimmed = True
        if r2_trim_end < r2_len:
            self._is_r2_primer_trimmed = True        
        self._is_r1_r2_overlap = primer_side_overlap or syn_side_overlap
        
        # polyA/T trimming if speRNA
        if self.seqtype == "rna":
            if self.poly_tail_primer_side != "none":
                r1_poly_tail_pos = self.poly_trim(r1_seq,self.poly_tail_primer_side)
                if r1_poly_tail_pos != -1:
                    self.r1_poly_tail_trim_len = len(r1_seq) - r1_poly_tail_pos
                    self._is_r1_poly_tail_trim = True
                    
                    r1_seq  = r1_seq[0:r1_poly_tail_pos]
                    r1_qual = r1_qual[0:r1_poly_tail_pos]
                    
            if self.poly_tail_umi_side != "none":
                r2_poly_tail_pos = self.poly_trim(r2_seq,self.poly_tail_umi_side)
                if r2_poly_tail_pos != -1:
                    self.r2_poly_tail_trim_len = len(r2_seq) - r2_poly_tail_pos
                    self._is_r2_poly_tail_trim = True
                    
                    r2_seq  = r2_seq[0:r2_poly_tail_pos]
                    r2_qual = r2_qual[0:r2_poly_tail_pos]
                    
        # final read length check
        if len(r1_seq) < self.min_primer_side_len or len(r2_seq) < self.min_umi_side_len:
            self._is_too_short = True
            self.synthetic_oligo_len = synthetic_oligo_len
            return
        
        # update read ids
        r1_id = self._reformat_readid(r1_id,umi,primer_info,primer_error)
        r2_id = self._reformat_readid(r2_id,umi,primer_info,primer_error)
        
        # reset variable
        self.synthetic_oligo_len = synthetic_oligo_len

        # final check for 0 --min_umi_side_len
        if len(r2_seq) == 0:
            r2_seq  = b"N"
            r2_qual = b"!"
            
        # update read info tuple
        self._r1_info = (r1_id,r1_seq,r1_qual)
        self._r2_info = (r2_id,r2_seq,r2_qual)
        
        assert len(r1_seq) != 0 or len(r2_seq) != 0 or len(r1_qual) != 0 or len(r2_qual) != 0,"id:{readid}\tR1:{r1}\tR2:{r2}".format(readid=r1_id.decode("ascii"),r1=r1_seq.decode("ascii"),r2=r2_seq.decode("ascii"))
        

def trim_custom_sequencing_adapter(args,buffers):
    ''' Returns whether to trim custom sequencing adapter or not
    :param buffers: 2-tuple of R1,R2 lines as byte strings
    :rtype tuple
    :returns (num_reads,True/False)
    '''
    trim_obj = QiaSeqTrimmer(
        is_nextseq                = args.is_nextseq,
        is_duplex                 = args.is_duplex,
        is_multimodal             = args.is_multimodal,
        seqtype                   = args.seqtype,
        max_mismatch_rate_primer  = args.max_mismatch_rate_primer,
        max_mismatch_rate_overlap = args.max_mismatch_rate_overlap,
        custom_seq_adapter        = args.custom_seq_adapter,
        umi_len                   = args.umi_len,
        umi_len_alt               = args.umi_len_alt,
        common_seq_len            = args.common_seq_len,
        overlap_check_len         = args.overlap_check_len,
        min_primer_side_len       = args.min_primer_side_len,
        min_umi_side_len          = args.min_umi_side_len,
        umi_filter_min_bq         = args.umi_filter_min_bq,
        umi_filter_max_lowQ_bases = args.umi_filter_max_lowQ_bases,
        umi_filter_max_Ns         = args.umi_filter_max_Ns,
        check_primer_side         = args.check_primer_side,
        trim_custom_seq_adapter   = False,
        primer3_R1                = args.primer3_bases_R1,
        primer3_R2                = args.primer3_bases_R2,
        poly_tail_primer_side     = args.poly_tail_primer_side,
        poly_tail_umi_side        = args.poly_tail_umi_side,
        tagname_duplex            = args.tagname_duplex,
        tagname_multimodal        = args.tagname_multimodal,
        tagname_umi               = args.tagname_umi,
        tagname_primer            = args.tagname_primer,
        tagname_primer_error      = args.tagname_primer_error,
        tag_separator             = args.tag_separator,
        field_separator           = args.field_separator,
        no_tagnames               = args.no_tagnames
    )
    
    buff_r1,buff_r2 = buffers
    r1_lines = buff_r1.split(b"\n")
    r2_lines = buff_r2.split(b"\n")

    num_reads = 0
    num_reads_have_adapter = 0
    i = 1
    for line in zip(r1_lines,r2_lines):
        if line[0] == "": # last element is empty because of the split("\n") above
            continue        
        if i % 4 == 2: # seq
            r1_seq,r2_seq = line
        elif i % 4 == 0: # qual
            num_reads+=1
            if trim_obj.custom_sequencing_adapter_check(r1_seq) != -1:
                num_reads_have_adapter += 1                
        i+=1
    return (num_reads,float(num_reads_have_adapter)/num_reads > 0.95)
    
def wrapper_func(args,buffer_):
    ''' This is the function called in parallel.
    Initilializes trimmer object and counters. Makes function calls
    and accumulates results for each batch/buffer
    
    :param ArgparseObject args: Command line arguments
    :param bytestring buffers_: A bytestring of length buffersize returned by the iterate_fastq function

    :rtype tuple
    :returns trimmed R1,R2 lines and metrics as a tuple
    '''
    trim_obj = QiaSeqTrimmer(
        is_nextseq                = args.is_nextseq,
        is_duplex                 = args.is_duplex,
        is_multimodal             = args.is_multimodal,
        seqtype                   = args.seqtype,
        max_mismatch_rate_primer  = args.max_mismatch_rate_primer,
        max_mismatch_rate_overlap = args.max_mismatch_rate_overlap,
        custom_seq_adapter        = args.custom_seq_adapter,
        umi_len                   = args.umi_len,
        umi_len_alt               = args.umi_len_alt,
        common_seq_len            = args.common_seq_len,
        overlap_check_len         = args.overlap_check_len,
        min_primer_side_len       = args.min_primer_side_len,
        min_umi_side_len          = args.min_umi_side_len,
        umi_filter_min_bq         = args.umi_filter_min_bq,
        umi_filter_max_lowQ_bases = args.umi_filter_max_lowQ_bases,
        umi_filter_max_Ns         = args.umi_filter_max_Ns,
        check_primer_side         = args.check_primer_side,
        trim_custom_seq_adapter   = args.to_trim_custom_adapter,
        primer3_R1                = args.primer3_bases_R1,
        primer3_R2                = args.primer3_bases_R2,
        poly_tail_primer_side     = args.poly_tail_primer_side,
        poly_tail_umi_side        = args.poly_tail_umi_side,
        tagname_duplex            = args.tagname_duplex,
        tagname_multimodal        = args.tagname_multimodal
        tagname_umi               = args.tagname_umi,
        tagname_primer            = args.tagname_primer,
        tagname_primer_error      = args.tagname_primer_error,
        tag_separator             = args.tag_separator,
        field_separator           = args.field_separator,        
        no_tagnames               = args.no_tagnames
    )
    
    # unpack input byte string                                 
    buff_r1,buff_r2 = buffer_
    r1_lines        = buff_r1.split(b"\n")
    r2_lines        = buff_r2.split(b"\n")
    
    # init some counters
    # general bookkeeping
    num_reads_after_trim  = 0    
    num_too_short         = 0
    num_odd               = 0
    num_bad_umi           = 0
    num_no_duplex         = 0
    num_no_UMIend_adapter = 0
    num_reads             = 0
    num_r1_primer_trimmed = 0
    num_r1_syn_trimmed    = 0
    num_r2_primer_trimmed = 0
    num_r1_r2_overlap     = 0    
    # duplex specific
    num_CC = 0
    num_TT = 0
    num_NN = 0
    # multimodal specific
    num_UMIend_B   = 0
    num_UMIend_RT  = 0
    num_UMIend_TSO = 0
    # quality trimming
    num_qual_trim_bases_r1 = 0
    num_qual_trim_bases_r2 = 0
    num_qual_trim_r1       = 0
    num_qual_trim_r2       = 0
    # polyA/T trimming
    num_poly_trim_bases_primer = 0
    num_poly_trim_bases_umi    = 0
    num_poly_trim_primer       = 0
    num_poly_trim_umi          = 0

    # store r1 and r2 lines
    out_lines_r1_bucket = []
    out_lines_r2_bucket = []
    
    i = 1
    for line in zip(r1_lines,r2_lines):
        if line[0] == "": # last element is empty because of the split("\n") above
            continue        
        if i % 4 == 1: # header
            r1_readid,r2_readid = line
        elif i % 4 == 2: # seq
            r1_seq,r2_seq = line
        elif i % 4 == 3: # placeholder
            pass
        elif i % 4 == 0: # qual
            num_reads += 1

            r1_qual,r2_qual = line
            # have R1 and R2 ready to process now
            if args.is_r2_primer_side:
                trim_obj.r1_info = (r2_readid,r2_seq,r2_qual)
                trim_obj.r2_info = (r1_readid,r1_seq,r1_qual)                
            else:
                trim_obj.r1_info = (r1_readid,r1_seq,r1_qual)
                trim_obj.r2_info = (r2_readid,r2_seq,r2_qual)

            trim_obj.qiaseq_trim(primer_datastruct)

            if trim_obj.is_r1_poly_tail_trim:
                num_poly_trim_bases_primer += trim_obj.r1_poly_tail_trim_len
                num_poly_trim_primer += 1
            if trim_obj.is_r2_poly_tail_trim:
                num_poly_trim_bases_umi += trim_obj.r2_poly_tail_trim_len
                num_poly_trim_umi += 1
                
            if trim_obj.is_too_short:
                num_too_short += 1
                i += 1
                continue
            elif trim_obj.is_odd_structure:
                num_odd += 1
                i += 1
                continue
            elif trim_obj.is_bad_umi:
                assert args.seqtype == "rna","UMI filter only applicable for speRNA reads!"
                num_bad_umi+=1
                i += 1
                continue
            
            if args.is_duplex and not trim_obj.is_duplex_adapter_present:
                num_no_duplex += 1
                i += 1
                continue
            
            if args.is_multimodal and not trim_obj.is_multimodal_adapter_present:
                num_no_UMIend_adapter += 1
                i += 1
                continue

            if trim_obj.is_r1_qual_trim:
                num_qual_trim_bases_r1 += trim_obj.r1_qual_trim_len
                num_qual_trim_r1 += 1
            if trim_obj.is_r2_qual_trim:
                num_qual_trim_bases_r2 += trim_obj.r2_qual_trim_len
                num_qual_trim_r2 += 1            
            
            # retrieve trimmed sequences
            if args.is_r2_primer_side:
                trimmed_r2_info = trim_obj.r1_info
                trimmed_r1_info = trim_obj.r2_info
            else:
                trimmed_r1_info = trim_obj.r1_info
                trimmed_r2_info = trim_obj.r2_info
            
            trimmed_r1_lines = trim_obj.field_separator.join([trimmed_r1_info[0],trimmed_r1_info[1],b"+",trimmed_r1_info[2]])
            trimmed_r2_lines = trim_obj.field_separator.join([trimmed_r2_info[0],trimmed_r2_info[1],b"+",trimmed_r2_info[2]])

            if trim_obj.is_r1_primer_trimmed:
                num_r1_primer_trimmed += 1
            if trim_obj.is_r1_syn_trimmed:
                num_r1_syn_trimmed += 1                
            if trim_obj.is_r2_primer_trimmed:
                num_r2_primer_trimmed += 1
            if trim_obj.is_r1_r2_overlap:
                num_r1_r2_overlap += 1
                
            if args.is_duplex:
                duplex_tag = trim_obj.duplex_tag
                if duplex_tag == b"CC":
                    num_CC += 1
                elif duplex_tag == b"TT":
                    num_TT += 1
                elif duplex_tag == b"NN":
                    num_NN += 1
                else:
                    raise Exception("Invalid duplex tag : {}".format(duplex_tag.decode("ascii")))

            if args.is_multimodal:
                multimodal_adapter_name = trim_obj.multimodal_adapter_name
                if multimodal_adapter_name == b"B":
                    num_UMIend_B += 1
                elif multimodal_adapter_name == b"RT":
                    num_UMIend_RT += 1
                elif multimodal_adapter_name == b"TSO":
                    num_UMIend_TSO += 1
                else:
                    raise Exception("Invalid UMI side common sequence (multimodal) : {}".format(multimodal_adapter_name.decode("ascii")))

            out_lines_r1_bucket.append(trimmed_r1_lines)
            out_lines_r2_bucket.append(trimmed_r2_lines)
            num_reads_after_trim += 1
            
        i += 1
        
    metrics = (num_r1_primer_trimmed,num_r1_syn_trimmed,num_r2_primer_trimmed,num_r1_r2_overlap,num_too_short,num_odd,num_bad_umi,num_reads,num_qual_trim_bases_r1,num_qual_trim_bases_r2,num_qual_trim_r1,num_qual_trim_r2,num_poly_trim_bases_primer,num_poly_trim_primer,num_poly_trim_bases_umi,num_poly_trim_umi,num_CC,num_TT,num_NN,num_no_duplex,num_UMIend_B,num_UMIend_RT,num_UMIend_TSO,num_no_UMIend_adapter,num_reads_after_trim)
    out_lines_R1 = b"\n".join(out_lines_r1_bucket)
    out_lines_R2 = b"\n".join(out_lines_r2_bucket)

    return (out_lines_R1,out_lines_R2,metrics)

            
def iterate_fastq(f,f2,ncpu,buffer_size=4*4*1024**2):
    ''' Copied from cutadapt, 
    added logic to yield a list of buffers equal to the number of CPUs
    :param file_handle f:  R1 fastq
    :param file_handle f2: R2 fastq
    :param int ncpu: length of buffer list 
    :param int buffer_size: size of each element of buffer list

    :yields list of length ncpu
    '''
    buf1 = bytearray(buffer_size)
    buf2 = bytearray(buffer_size)
    
    # Read one byte to make sure we are processing FASTQ
    start1 = f.readinto(memoryview(buf1)[0:1])
    start2 = f2.readinto(memoryview(buf2)[0:1])

    if (start1 == 1 and buf1[0:1] != b'@') or (start2 == 1 and buf2[0:1] != b'@'):
        raise Exception('Paired-end data must be in FASTQ format when using multiple cores')

    to_yield = []
    nchunks = 0
    while True:
        bufend1 = f.readinto(memoryview(buf1)[start1:]) + start1
        bufend2 = f2.readinto(memoryview(buf2)[start2:]) + start2
        if start1 == bufend1 and start2 == bufend2:
            break

        end1, end2 = two_fastq_heads(buf1, buf2, bufend1, bufend2)
        assert end1 <= bufend1
        assert end2 <= bufend2

        if end1 > 0 or end2 > 0:
            nchunks+=1
            to_yield.append((memoryview(buf1)[0:end1].tobytes(), memoryview(buf2)[0:end2].tobytes()))
            
        start1 = bufend1 - end1
        assert start1 >= 0
        buf1[0:start1] = buf1[end1:bufend1]
        start2 = bufend2 - end2
        assert start2 >= 0
        buf2[0:start2] = buf2[end2:bufend2]
        
        if nchunks == ncpu:
            yield to_yield
            to_yield = []
            nchunks = 0

    if start1 > 0 or start2 > 0:
        to_yield.append((memoryview(buf1)[0:start1].tobytes(), memoryview(buf2)[0:start2].tobytes()))
    if len(to_yield) > 0:
        yield to_yield

def init(l):
    global lock
    lock = l

def open_fh(fname1,fname2,read=True):
    ''' Return appropriate file handles
    :param str fname1: R1 fastq file name
    :param str fname2: R2 fastq file name
    :param bool read: rb/wb mode
    :rtype tuple
    :returns tuple of file handles
    '''
    mode = "rb" if read else "wb"
    if fname1.endswith(".gz"):
        return (gzip.open(fname1,mode),gzip.open(fname2,mode))
    else:
        return (open(fname1,mode),open(fname2,mode))
    
def close_fh(fh1,fh2):
    ''' Closes file handles
    :param str fh1: R1 file handle
    :param str fh2: R2 file handle
    '''
    fh1.close()
    fh2.close()
    
def main(args):
    ''' Calls wrapper func for trimming in parallel using multiprocessing,
    aggregates metrics and writes trimmed fastq and metric files
    '''
    # unpack some params
    primer_file = args.primer_file
    r1          = args.r1
    r2          = args.r2
    out_r1      = args.out_r1
    out_r2      = args.out_r2
    out_metrics = args.out_metrics
    # convert to bytestring
    args.tagname_primer       = args.tagname_primer.encode("ascii")
    args.tagname_primer_error = args.tagname_primer_error.encode("ascii")    
    args.tagname_umi          = args.tagname_umi.encode("ascii")
    args.tagname_duplex       = args.tagname_duplex.encode("ascii")
    args.tagname_multimodal   = args.tagname_multimodal.encode("ascii")
    args.tag_separator        = args.tag_separator.encode("ascii")
    args.field_separator      = args.field_separator.encode("ascii")
    args.custom_seq_adapter   = args.custom_seq_adapter.encode("ascii")

    global primer_datastruct
    
    primer_datastruct = PrimerDataStruct(k=8,primer_file=args.primer_file,
                                         ncpu=args.ncpu,seqtype=args.seqtype,
                                         primer_col=args.primer_col).primer_search_datastruct
    l = multiprocessing.Lock()

    logger.info("\n{}\n".format("--"*10))
    logger.info("Created Primer Datastruct\n")

    # check custom sequencing adapter in the first chunk of the fastq ~ 11000 reads
    if not args.trim_custom_seq_adapter:
        f,f2 = open_fh(r1,r2)
        for buffers in iterate_fastq(f,f2,1):
            assert len(buffers) == 1
            num_reads, to_trim_custom_adapter = trim_custom_sequencing_adapter(args,buffers[0])
            break
    
        logger.info("Checked first {n} reads for custom sequencing adapter.".format(n=num_reads))
        if to_trim_custom_adapter:
            logger.info("Custom sequencing adapter present in > 95% reads")
        args.to_trim_custom_adapter = to_trim_custom_adapter
        close_fh(f,f2)
    else: # force trimming of custom sequencing adapter , useful for ion reads
        args.to_trim_custom_adapter = args.trim_custom_seq_adapter
        
    logger.info("\nRunning program with args : {}\n".format(args))

    num_r1_primer_trimmed = 0
    num_r1_syn_trimmed    = 0
    num_r2_primer_trimmed = 0
    num_r1_r2_overlap     = 0
    num_too_short         = 0
    num_odd               = 0
    num_bad_umi           = 0
    num_no_duplex         = 0
    num_no_UMIend_adapter = 0
    total_reads           = 0
    num_after_trim        = 0
    
    num_qual_trim_r1_bases = 0
    num_qual_trim_r2_bases = 0
    num_qual_trim_r1 = 0
    num_qual_trim_r2 = 0

    num_poly_trim_bases_primer  = 0
    num_poly_trim_primer        = 0
    num_poly_trim_bases_umi     = 0
    num_poly_trim_umi           = 0 
    
    num_CC = 0
    num_TT = 0
    num_NN = 0    
    
    num_UMIend_TSO = 0
    num_UMIend_RT  = 0
    num_UMIend_B   = 0
    
    nchunk = 1

    f,f2 = open_fh(r1,r2)
    f_out_r1,f2_out_r2 = open_fh(out_r1,out_r2,read=False)
    
    p = multiprocessing.Pool(args.ncpu,initializer=init, initargs=(l,))
    func = functools.partial(wrapper_func,args)
    
    for chunks in iterate_fastq(f,f2,args.ncpu):
        res = p.map(func,chunks)

        for trimmed_r1_lines,trimmed_r2_lines,counter_tup in res:

            # unpack counters
            num_r1_primer_trimmed   +=  counter_tup[0]
            num_r1_syn_trimmed      +=  counter_tup[1]
            num_r2_primer_trimmed   +=  counter_tup[2]
            num_r1_r2_overlap       +=  counter_tup[3]
            num_too_short           +=  counter_tup[4]
            num_odd                 +=  counter_tup[5]
            num_bad_umi             +=  counter_tup[6]
            total_reads             +=  counter_tup[7]
            
            num_qual_trim_r1_bases  +=  counter_tup[8]
            num_qual_trim_r2_bases  +=  counter_tup[9]        
            num_qual_trim_r1        +=  counter_tup[10]
            num_qual_trim_r2        +=  counter_tup[11]

            num_poly_trim_bases_primer += counter_tup[12]
            num_poly_trim_primer       += counter_tup[13]
            num_poly_trim_bases_umi    += counter_tup[14]
            num_poly_trim_umi          += counter_tup[15]            
            
            num_CC         +=  counter_tup[16]
            num_TT         +=  counter_tup[17]
            num_NN         +=  counter_tup[18]
            num_no_duplex  +=  counter_tup[19]
            
            num_UMIend_B          += counter_tup[20]
            num_UMIend_RT         += counter_tup[21]
            num_UMIend_TSO        += counter_tup[22]
            num_no_UMIend_adapter += counter_tup[23]

            num_after_trim +=  counter_tup[24]

            f_out_r1.write(trimmed_r1_lines)
            f_out_r1.write(b"\n")
            f2_out_r2.write(trimmed_r2_lines)
            f2_out_r2.write(b"\n")
            
            nchunk += 1
            logger.info("Processed : {n} reads".format(n=total_reads))
            
    p.close()
    p.join()

    thousand_comma = lambda x: "{:,}".format(x) if args.thousand_comma else str(x) 
    
    out_metrics = [
        thousand_comma(total_reads) + "\tTotal read fragments",
        thousand_comma(num_r1_primer_trimmed) + "\tNum Primer side reads with primer identified",
        thousand_comma(num_r1_syn_trimmed) + "\tNum Primer side reads with 3' synthetic oligo trimmed",
        thousand_comma(num_r2_primer_trimmed) + "\tNum UMI side reads with 3' synthetic oligo trimmed (including primer)",
        thousand_comma(num_r1_r2_overlap) + "\tNum read fragments overlapping",
        thousand_comma(num_qual_trim_r1) + "\tNum Primer side reads qual trimmed",
        thousand_comma(num_qual_trim_r2) + "\tNum UMI side reads qual trimmed",        
        "{qual_trim_r1}\tAvg num Primer side bases qual trimmed",
        "{qual_trim_r2}\tAvg num UMI side bases qual trimmed",
    ]
    out_metrics_dropped = [
        thousand_comma(num_too_short) + "\tNum read fragments dropped too short",
        thousand_comma(num_odd) + "\tNum read fragments dropped odd structure (usually Primer side has only primer sequence)"
    ]
    total_dropped = num_too_short + num_odd
    
    if args.is_duplex:
        out_metrics.extend(
            [thousand_comma(num_CC) + "\tNum CC reads",
             thousand_comma(num_TT) + "\tNum TT reads",
             thousand_comma(num_NN) + "\tNum NN reads"])
        out_metrics_dropped.append(thousand_comma(num_no_duplex) + "\tNum read fragments dropped (no duplex adapter)")
        total_dropped += num_no_duplex

    if args.is_multimodal:
        out_metrics.extend(
            [thousand_comma(num_UMIend_B)   + "\tNum UMI end reads with common sequence B",
             thousand_comma(num_UMIend_RT)  + "\tNum UMI end reads with common sequence RT",
             thousand_comma(num_UMIend_TSO) + "\tNum UMI end reads with common sequence TSO"])
        num_dropped_leakage = num_UMIend_B if args.seqtype == "rna" else num_UMIend_RT + num_UMIend_TSO
        out_metrics_dropped.append(thousand_comma(num_dropped_leakage) + "\tNum read fragments dropped (no/wrong common adapter)")
        total_dropped += num_dropped_leakage

    if args.seqtype == "rna":
        if args.poly_tail_primer_side != "none":
            out_metrics.extend(
                ["{p1_trim}\tAvg num bases {p1} trimmed Primer side".format(p1 = args.poly_tail_primer_side,
                                                                            p1_trim = 0 if num_poly_trim_bases_primer == 0 else \
                                                                            round(float(num_poly_trim_bases_primer)/num_poly_trim_primer,2)),
                 thousand_comma(num_poly_trim_primer) + "\tNum reads {p1} trimmed Primer side".format(p1 = args.poly_tail_primer_side)])
        if args.poly_tail_umi_side != "none":
            out_metrics.extend(
                ["{p2_trim}\tAvg num bases {p2} trimmed UMI side".format(p2 = args.poly_tail_umi_side,
                                                                         p2_trim = 0 if num_poly_trim_bases_umi == 0 else \
                                                                         round(float(num_poly_trim_bases_umi)/num_poly_trim_umi,2)),
                 thousand_comma(num_poly_trim_umi) + "\tNum reads {p2} trimmed UMI side".format(p2 = args.poly_tail_umi_side)])


        out_metrics_dropped.append(thousand_comma(num_bad_umi) + "\tNum read fragments dropped bad UMI")

        total_dropped += num_bad_umi       

    out_metrics.extend(out_metrics_dropped)                                                                          
    out_metrics.append(thousand_comma(num_after_trim) + "\tNum read fragments after trimming")
    out_metrics.append("{pct_after_trim} \tPct read fragments after trimming".format(pct_after_trim = round(100*num_after_trim/total_reads, 2) if total_reads else 0.00))
    
    out_metrics_lines = "\n".join(out_metrics).format(qual_trim_r1 = 0 if num_qual_trim_r1_bases == 0 else \
                                                      round(float(num_qual_trim_r1_bases)/(num_qual_trim_r1),2),
                                                      qual_trim_r2 = 0 if num_qual_trim_r2_bases == 0 else \
                                                      round(float(num_qual_trim_r2_bases)/(num_qual_trim_r2),2))
    
    f_out_metrics = open(args.out_metrics,"w")                                                      
    f_out_metrics.write(out_metrics_lines)
    f_out_metrics.write("\n")   
    
    close_fh(f,f2)
    close_fh(f_out_r1,f2_out_r2)
    f_out_metrics.close()
        
    assert num_after_trim + total_dropped == total_reads, "Read accouting failed !"
    
if __name__ == '__main__':
    main()
