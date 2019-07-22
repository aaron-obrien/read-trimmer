import unittest
import sys
import os
from argparse import Namespace

sys.path.append(os.path.join(
    os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))) , "trimmer/"))

from trimmer import PrimerDataStruct
import helper

# to do :
# 1. Need to add more test cases to all routines here
# 2. Need example fastq files to simulate command line running
# 3. Add some rna and duplex specific test cases


class TestQiaSeqTrimDna(unittest.TestCase):
    '''
    '''
    def setUp(self):
        '''
        '''
        # primer datastruct
        self.primer_datastruct = PrimerDataStruct(k=8,primer_file=os.path.join(os.path.dirname(__file__),"test_data/test_primers_dna.txt"),
                                                  ncpu=1,seqtype="dna",primer_col=3).primer_search_datastruct

        # argparse obj
        self.args = helper.helper_return_args()
        # trimmer obj
        self.trimmer_obj = helper.helper_return_qiaseq_obj(self.args)

    def test_overlap_synthetic_side_check(self):
        ''' Check the default R1/R2 overlap using the Synthetic/UMI side sequence
        '''
        # edit dist = 0
        r1_seq = b"CCCAGGTCACCATCAAATACATCGGAGCCAGCCCCTTCGGAGGGTGCCAGTGGAGACCTGGGGGCCTCCTCTTCAGAGGGCTCCAGCCCTAGTGTCAGGTCCCCACCGCAGGACTCCAATTATCCCGCTTAT"
        r2_seq = b"ATAAGCGGGATAATTGGAGTCCTGCGGTGGGGACCTGACACTAGGGCTGGAGCCCTCTGAAGAGGAGGCCCCCAGGTCTCCACTGGCACCCTCCGAAGGGGCTGGCTCCGATGTATTTGATGGTGACCTGGGCAAAACGCAATACTGTACA"        
        self.assertEqual(self.trimmer_obj.synthetic_side_check(r1_seq,r2_seq),(84,108)) # checking (start,end) pos of overlap
        r1_seq = b"CTGCAGCCCCAACCGACAGAAGCCCGTGGTAGACCATTCTGTGCGGATCAATTCTGTCGGCAGCACCGCTTCCTCCTCCCAGCCTCTGCTTGTGCACGACGATGTCTGAGCAGAATCAGTGTTTGGGTCACCCCTCCAGGAATGATCAGGA"
        r2_seq = b"GTTGTGTATATAATTGGAGTCCTGATCATTCCTGGAGGGGTGACCCAAACACTGATTCTGCTCAGACATCGTCGTGCACAAGCAGAGGCTGGGAGGAGGAAGCGGTGCTGCCGACAGAATTGATCCGCACAGAATGGTCTACCACGGGCTT"
        self.assertEqual(self.trimmer_obj.synthetic_side_check(r1_seq,r2_seq),(122,146))
        
        # edit dist = 1
        r1_seq = b"TGGAGATGTGATAATTTCAGGAAACAAAAATTTGTGCTATGCAAATACAATAAACTGGGAAAAACTGTTTGGGACCTCCGGAGGACTCCAATTACATCCCTTTT"
        r2_seq = b"AAAAGGGATGTAATTGGAGTCCTCCGGAGGTCCCAAACAGTTTTTTCCAGTTTATTGTATTTGCATAGCACAAATTTTTGTTTCCTGAAATTATCACATCTCCACAAAACGCAATACTGTACATT"
        self.assertEqual(self.trimmer_obj.synthetic_side_check(r1_seq,r2_seq),(56,80))
        r1_seq = b"CCAGGCCAAAGTCACAGATCTTCACAATTTTTCCTTGTGCCAGGAGAACGTTGCGAGCAGCCAGATCACGGTGGACACACTGCATGGAAAAGGAAGAAATGACTCAGGATCAAGCCATCTGTAGCTGAGTAGCCCAACCTATGCGATC"
        r2_seq = b"CGGAACGACTAGATTGGAGTCCTCAGCTACAGATGGCTTGAGCCTGAGTCATTTCTTCCTTTTCCATGCAGTGCGTCCACCGATATCTGGCTGCTCGCAACGTTCTCCTGGCACAAGGAAAAATTGTGAAGATCTGTGACTTTGGCCTGG"
        self.assertEqual(self.trimmer_obj.synthetic_side_check(r1_seq,r2_seq),(102,126))

        # edit dist = 3
        r1_seq = b"CTTGTTTCCCACTAGCACCATAGGTACATCATCCGAGTCTTTTACTCGCTTAATCTGCTCCCTAAAAACGGGAATATATTATCAGAACATAAGAAAAACAAGATTAGGCTGGGTACAGTGGCTCATGCCTGCAATCTCAGCACTTTGGGAG"
        r2_seq = b"GTTACCCGAAGAATTGGAGTCCTTCCCAAAGTGCTGGGAATACAGGCATGAGCCACTGCACCCAGCCTAATCTTTGTATTTTTTTGTAGAGACCGGGTTTTGCCATGTTGCCCAGGCTAATCTCAAACTCCTGGGTTCAAGCAGTCTGCC"
        self.assertEqual(self.trimmer_obj.synthetic_side_check(r1_seq,r2_seq),(125,149))

        # no overlap on synthetic side
        r1_seq = b"ATCAAGTGGATGGCGCTGGAGTCCCAGCTGAGCTTGGAAACTATAAATCTTAGCTTCCTCTTAGTTTCAAGTCACTTATAAATGAAGAGGAGGAAGAACCTGGGAACCTGAAGGGAAGATGGGGGCTGTGCATTAGGG"
        r2_seq = b"GCAGTTTATTGGATTGGAGTCCTAGACAGTCACTTGGGTTGCCACATTGGCACTTGTACTGCAGAAGTCCCTTGAGCACCTATGACAACAGGCTGGAGCCGCCTATTCTTGCCCCCGCCCTAATGCACAGCCCCCATCTTCCCTTCAGGTT"
        self.assertEqual(self.trimmer_obj.synthetic_side_check(r1_seq,r2_seq),(-1,-1))

        # edit dist = 4, do not count this as an overlap as cut off is 0.12 mismatch rate
        r1_seq = b"ATCACTTTGCGCGGTGTAGATATGATCAAAAAGGGATTCAATCACCATCCATTTAACTGGCATCCGAACGACTCCACTTTCCTAATTTAC"
        r2_seq = b"GTAAATTAGGAAATTGGAGTCCTTCGGATTCCAGCTAAATGGATGGCAATTGAATCCCTTTTTGTTCATATCTACACCACGCAAAGTGATCAAAACGCAATACTGTACATT"
        self.assertEqual(self.trimmer_obj.synthetic_side_check(r1_seq,r2_seq),(-1,-1))

    def test_primer_trim(self):
        ''' Some checks done as part of the test_overlap_primer_side_check below. Need seperate tests here.
        '''
        pass

    def test_multimodal_polyT_5prime(self):
        ''' Check 5' polyT trim of UMI side reads
        '''
        temp_args                            = self.args
        temp_args.is_multimodal              = True
        temp_args.trim_polyT_5prime_umi_side = True
        temp_trimmer_obj = helper.helper_return_qiaseq_obj(temp_args)

        temp_trimmer_obj._is_r2_polyT_5prime_trim = False
        r2_seq = b"GTTAATCGCCACGTTTTTTTTTTTTTTTTTTGTTTTTTTTTTTCTTTTTTTCCTAACACGTGCGCTCCTGCTCCACCCCCCCGGGCCGGCGGGCGGGACGGGCCGGTGGGGCGCCCTTGGGGGGCTTGGGAGGCCCCGGGGTCCCACCTC"
        temp_trimmer_obj.synthetic_oligo_len = 10 + len(b"ACGTTTTTTTTTTTTTTTTTTVN")
        temp_trimmer_obj._multimodal_adapter_name = b"RT"
        temp_trimmer_obj._r2_polyT_5prime_trim_wrap(r2_seq)
        self.assertEqual(temp_trimmer_obj._is_r2_polyT_5prime_trim, True)
        self.assertEqual(temp_trimmer_obj.synthetic_oligo_len, len(b"GTTAATCGCCACGTTTTTTTTTTTTTTTTTTGTTTTTTTTTTT"))
        self.assertEqual(temp_trimmer_obj.r2_polyT_5prime_trim_len, len(b"TTTTTTTTTT"))

        temp_trimmer_obj._is_r2_polyT_5prime_trim = False
        temp_trimmer_obj.synthetic_oligo_len = 10 + len(b"ACGTTTTTTTTTTTTTTTTTTVN")
        temp_trimmer_obj._multimodal_adapter_name = b"TSO"
        temp_trimmer_obj._r2_polyT_5prime_trim_wrap(r2_seq)
        self.assertEqual(temp_trimmer_obj._is_r2_polyT_5prime_trim, False)

        # only 4 polyTs after adapter - should not trim this
        temp_trimmer_obj._is_r2_polyT_5prime_trim = False
        temp_trimmer_obj.synthetic_oligo_len = 10 + len(b"ACGTTTTTTTTTTTTTTTTTTVN")
        temp_trimmer_obj._multimodal_adapter_name = b"RT"
        r2_seq = b"GTTAATCGCCACGTTTTTTTTTTTTTTTTTTGTTTTTCTTTTTTTCCTAACACGTGCGCTCCTGCTCCACCCCCCCGGGCCGGCGGGCGGGACGGGCCGGTGGGGCGCCCTTGGGGGGCTTGGGAGGCCCCGGGGTCCCACCTC"
        temp_trimmer_obj._r2_polyT_5prime_trim_wrap(r2_seq)
        self.assertEqual(temp_trimmer_obj._is_r2_polyT_5prime_trim, False)

        # 5 Ts - should trim this
        temp_trimmer_obj._is_r2_polyT_5prime_trim = False
        temp_trimmer_obj.synthetic_oligo_len = 10 + len(b"ACGTTTTTTTTTTTTTTTTTTVN")
        temp_trimmer_obj._multimodal_adapter_name = b"RT"
        r2_seq = b"GTTAATCGCCACGTTTTTTTTTTTTTTTTTTGTTTTTTCTTTTTTTCCTAACACGTGCGCTCCTGCTCCACCCCCCCGGGCCGGCGGGCGGGACGGGCCGGTGGGGCGCCCTTGGGGGGCTTGGGAGGCCCCGGGGTCCCACCTC"
        temp_trimmer_obj._r2_polyT_5prime_trim_wrap(r2_seq)
        self.assertEqual(temp_trimmer_obj._is_r2_polyT_5prime_trim, True)
        self.assertEqual(temp_trimmer_obj.synthetic_oligo_len, len(b"GTTAATCGCCACGTTTTTTTTTTTTTTTTTTGTTTTTT"))
        self.assertEqual(temp_trimmer_obj.r2_polyT_5prime_trim_len, len(b"TTTTT"))

        # all Ts after adapter - should trim all the Read
        temp_trimmer_obj._is_r2_polyT_5prime_trim = False
        temp_trimmer_obj.synthetic_oligo_len = 10 + len(b"ACGTTTTTTTTTTTTTTTTTTVN")
        temp_trimmer_obj._multimodal_adapter_name = b"RT"
        r2_seq = b"GTTAATCGCCACGTTTTTTTTTTTTTTTTTTGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
        temp_trimmer_obj._r2_polyT_5prime_trim_wrap(r2_seq)
        self.assertEqual(temp_trimmer_obj._is_r2_polyT_5prime_trim, True)
        self.assertEqual(temp_trimmer_obj.synthetic_oligo_len, len(r2_seq))
        self.assertEqual(temp_trimmer_obj.r2_polyT_5prime_trim_len, len(b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"))


    def test_overlap_primer_side_check(self):
        ''' Test whether the --primer-side-check flag enables the R1/R2 overlap rescue
        '''
        # edit dist = 0
        r1_seq = b"GCCTCTTGCTTCTCTTTTCCTATCCTGAGTAGTGGTAATCTACTGGGACGGAACAGCTTTGAGGTGCGTGTTTGTGCCTGTCCTGGGAGAGACCGGCGCACAGAGGAAGAGAATCTCCGCAAGAAAGGGGAGCCTCACCACGAGCTGCCC"
        r2_seq = b"AGGCAACTAAACATTGGAGTCCTCTGGGGGCAGCTCGTGGTGAGGCTCCCCTTTCTTGCGGAGATTCTCTTCCTCTGTGCGCCGGTCTCTCCCAGGACAGGCACAAACACGCACCTCAAAGCTGTTCCGTCCCAGTAGATTACCACTACT"
        primer_end,primer,primer_editdist,primer_score = self.trimmer_obj.primer_trim(self.primer_datastruct,r1_seq)
        self.assertEqual(self.trimmer_obj.primer_side_check(primer_end,r1_seq,r2_seq),(119,143))
        
        # edit dist = 1
        r1_seq = b"TCCTGTGTATAAAAAGATAGCTAAATTCATGCATCATAAGCTCATTAATACTCTTCCTTACCATCCCCATTTTTAAAGATGATCTCATTGTTCTGAAACAGTAACTCTGACATGATGTCTGGGTTCTCCCAATTCAACCACAGTGGCCTTT"
        r2_seq = b"TTGCAGTAATACATTGGAGTCCTCGAATTATGTCCTCTGCAAAAAGGCCACTGTGGTTGAATTGGGAGAACCCAGACATCATGTCAGAGTTACTGTTTCAGAACAATGAGATCATCTTTAAAAATGGGGATGGTAAGGAAGAGTATTAATG"
        primer_end,primer,primer_editdist,primer_score = self.trimmer_obj.primer_trim(self.primer_datastruct,r1_seq)
        self.assertEqual(self.trimmer_obj.primer_side_check(primer_end,r1_seq,r2_seq),(127,150))
        
        # edit dist = 3
        r1_seq = b"ACAGAAAGCCCTGTAGAGCATCCATGAAATCTGGTCGCCTCATTTGCTCAACTAAAAACTTCAACTGTACCTGATTAAAAAAAAAGGATAGTCACAGTAAGAAACTGACTTAATAAAAGATAAAACTAACTGCTAATCTCTTTCTCACATG"
        r2_seq = b"AGGGCAAAGATGATTGGAGTCCTTGAGGTTCTCATGTGAGAAAGAGATTAGCAGTTAGTTTTATCTTTTATTAAGTCAGTTTCTTACTGTGACTATCCTTTTTTTTTAATCAGGTACAGTTGAAGTTTTTAGTTGAGCAAATGAGGCGAC"
        primer_end,primer,primer_editdist,primer_score = self.trimmer_obj.primer_trim(self.primer_datastruct,r1_seq)        
        self.assertEqual(self.trimmer_obj.primer_side_check(primer_end,r1_seq,r2_seq),(128,149))

        # edit dist = 4 , no overlap
        r1_seq = b"CACAATGGCACGGTTGAATGTAAGGCTTACAACGATGTGGGCAAGACTTCTGCCTATTTTAACTTTGCATTTAAAGGTAACAACAAAGGTATATTTCTTTTTAATCCAATTTAAGGGGATGTTTAGGCTCTGTCTACCATATCAGTCATGA"
        r2_seq = b"GGCAGGTATTAAATTGGAATCCTCTTAAAATCATGACTGATATGGTAGACAGAGCCTAAACATCCCCTTAAATTGGATTAAAAAGAAATATACCTTTGTTGTTACCTTTAAATGCAAAGTTAAAATAGGCAGAAGTCTTGCCCACATCGTT"
        primer_end,primer,primer_editdist,primer_score = self.trimmer_obj.primer_trim(self.primer_datastruct,r1_seq)
        self.assertEqual(self.trimmer_obj.primer_side_check(primer_end,r1_seq,r2_seq),(-1,-1))

        # edit dist > 4 , no overlap
        r1_seq = b"ACTCTTGCCTACGCCACCAGCTCCCGGGTTCAAGTGATTCTCCTGCCTCAGACTCCTGAGTAGCTGGTATTACAGGCGCACACCACCACACCCAGCTAATTTTTTTTTTTGTATTTTTAGTAGAGATGGGGTTTCACCATGTTGGCCAGG"
        r2_seq = b"GAGAGCGAAACAATTGGAGTCCTGTAATCCCAGCACTTTGGGAGGCTGAGGTGGGTGGATCATCTGAGGTCAGGAGTTTTGAGACCAGCCTGGCCAACATGGTGAAACCCCATCTCTACTAAAAATACAAAAAAAAAAATTAGCTGGGTG"
        primer_end,primer,primer_editdist,primer_score = self.trimmer_obj.primer_trim(self.primer_datastruct,r1_seq)
        self.assertEqual(self.trimmer_obj.primer_side_check(primer_end,r1_seq,r2_seq),(-1,-1))
        
    def test_qiaseq_trim(self):
        ''' Full trimming routine
        '''
        # sucessful overlap on both sides
        self.trimmer_obj.r1_info = (
            b"@M01750:502:000000000-AUK2G:1:1101:12973:1963 1:N:0:2",
            b"CCCAGGTCACCATCAAATACATCGGAGCCAGCCCCTTCGGAGGGTGCCAGTGGAGACCTGGGGGCCTCCTCTTCAGAGGGCTCCAGCCCTAGTGTCAGGTCCCCACCGCAGGACTCCAATTATCCCGCTTAT",
            b"CDCCCDFFFFFFGGGGGGGGGGHGGGGGGHHGGGGHHHHGGGGGEFFHHHHHHHHHGGHHGGGGGGGHHHHHHHHHGHGGGHHGHHHGGHHHHHHHHHHHHHHHHGGGGGGGGHHHHHHHHHHHHHGGGGGH")
        self.trimmer_obj.r2_info = (
            b"@M01750:502:000000000-AUK2G:1:1101:12973:1963 2:N:0:2",
            b"ATAAGCGGGATAATTGGAGTCCTGCGGTGGGGACCTGACACTAGGGCTGGAGCCCTCTGAAGAGGAGGCCCCCAGGTCTCCACTGGCACCCTCCGAAGGGGCTGGCTCCGATGTATTTGATGGTGACCTGGGCAAAACGCAATACTGTACA",
            b"ABBBBFBBBDAFGGGGGGGGGGHHGGGEEGGGGGGHHGHHHHHHHHHHGGGGFGGHGHHHHHHHHHGGGGGGGGGGGHHHHHHHHHHHHGGHGGGEFGGGGGGGGGGHGGGGGCHHHHHHHHHDFHHHGHHGHHGHHGGCGFGHGHHFHH0")
        self.trimmer_obj.qiaseq_trim(self.primer_datastruct)
        self.assertEqual(self.trimmer_obj.r1_info,(
            b"@M01750:502:000000000-AUK2G:1:1101:12973:1963\tmi:Z:ATAAGCGGGATA\tpr:Z:chr17-1-37883632-28\tpe:Z:0",
            b"ATCGGAGCCAGCCCCTTCGGAGGGTGCCAGTGGAGACCTGGGGGCCTCCTCTTCAGAGGGCTCCAGCCCTAGTGTCAGGTCCCCACCGC",
            b"GGHGGGGGGHHGGGGHHHHGGGGGEFFHHHHHHHHHGGHHGGGGGGGHHHHHHHHHGHGGGHHGHHHGGHHHHHHHHHHHHHHHHGGGG"))
        self.assertEqual(self.trimmer_obj.r2_info,(
            b"@M01750:502:000000000-AUK2G:1:1101:12973:1963\tmi:Z:ATAAGCGGGATA\tpr:Z:chr17-1-37883632-28\tpe:Z:0",
            b"GCGGTGGGGACCTGACACTAGGGCTGGAGCCCTCTGAAGAGGAGGCCCCCAGGTCTCCACTGGCACCCTCCGAAGGGGCTGGCTCCGAT",
            b"HGGGEEGGGGGGHHGHHHHHHHHHHGGGGFGGHGHHHHHHHHHGGGGGGGGGGGHHHHHHHHHHHHGGHGGGEFGGGGGGGGGGHGGGG"))

        # umi side overlap fail , primer side sucessful
        self.trimmer_obj.r1_info = (
            b"@M01750:502:000000000-AUK2G:1:1101:15325:1360 2:N:0:2",
            b"CTTGCTCTGATAGGAAAATGAGATCTACTGTTTTCCTTTACTTACTACACCTCAGATATATTTCTTCATGAAGACCTCATAGTAAAAATAGGTGATTTTGGTCTAGCTACAGTGAAATCTCGATGGAGTGGGTCCCATCAGTTTGAACAGT",
            b"BABBBFFFFFFFGGFFG4FFFGHGGHHHFFGF5FEGGHGFHFCHEFF5AAAGBFH3AADHGHFHFBGBFGFHHFEHBGFF5HDFGHFEGHFFDFFFBG5G??FD5FGH5GHFGHFHHBFGHHGFF?CGGFAGEHGHBGHH3?FBDBEGHG?")
        self.trimmer_obj.r2_info = (
            b"@M01750:502:000000000-AUK2G:1:1101:15325:1360 2:N:0:2",
            b"AAGTATAGGTGGATTGGAGTCCTCAGACAACTGTTCAAACTGATGGGACCCACTCCATCGAGATTTCACTGTAGCTAGACCAAAATCACCTATTTTTACTATGAGGTCTTCATGAAGAAATATATCTGAGGAGTAGTAAGTAAAGGAAAA",
            b"1>A1>D@FFB1>GGGGGGGGGGHHH1110BG1AA33DGGFEHH211EFGG?GG00B1F1AF//AEBGHGFF22D2AGHGGHGHC?F1110F1B122F/DFFFGEEGHHH1FHGHEBG1>BF1@BFEDF>F/1?BEFD2BFE2G2F1BBGE")
        self.trimmer_obj.qiaseq_trim(self.primer_datastruct)
        self.assertEqual(self.trimmer_obj.r1_info,(
            b"@M01750:502:000000000-AUK2G:1:1101:15325:1360\tmi:Z:AAGTATAGGTGG\tpr:Z:chr7-1-140453212-36\tpe:Z:0",
            b"TGTTTTCCTTTACTTACTACACCTCAGATATATTTCTTCATGAAGACCTCATAGTAAAAATAGGTGATTTTGGTCTAGCTACAGTGAAATCTCGATGGAGTGGGTCCCATCAGTTTGAACAGT",
            b"FFGF5FEGGHGFHFCHEFF5AAAGBFH3AADHGHFHFBGBFGFHHFEHBGFF5HDFGHFEGHFFDFFFBG5G??FD5FGH5GHFGHFHHBFGHHGFF?CGGFAGEHGHBGHH3?FBDBEGHG?"))
        self.assertEqual(self.trimmer_obj.r2_info,(
            b"@M01750:502:000000000-AUK2G:1:1101:15325:1360\tmi:Z:AAGTATAGGTGG\tpr:Z:chr7-1-140453212-36\tpe:Z:0",
            b"CAGACAACTGTTCAAACTGATGGGACCCACTCCATCGAGATTTCACTGTAGCTAGACCAAAATCACCTATTTTTACTATGAGGTCTTCATGAAGAAATATATCTGAGGAGTAGTAAGTAAAGGAAAA",
            b"HH1110BG1AA33DGGFEHH211EFGG?GG00B1F1AF//AEBGHGFF22D2AGHGGHGHC?F1110F1B122F/DFFFGEEGHHH1FHGHEBG1>BF1@BFEDF>F/1?BEFD2BFE2G2F1BBGE"))        

        self.trimmer_obj.r1_info = (
            b"@M01750:502:000000000-AUK2G:1:1101:15325:1360 2:N:0:2",
            b"CTTGCTCTGATAGGAAAATGAGATCTACTGTTTTCCTTTACTTACTACACCTCAGATATATTTCTTCATGAAGACCTCATAGTAAAAATAGGTGATTTTGGTCTAGCTACAGTGAAATCTCGATGGAGTGGGTCCCATCAGTTTGAACAGT",
            b"BABBBFFFFFFFGGFFG4FFFGHGGHHHFFGF5FEGGHGFHFCHEFF5AAAGBFH3AADHGHFHFBGBFGFHHFEHBGFF5HDFGHFEGHFFDFFFBG5G??FD5FGH5GHFGHFHHBFGHHGFF?CGGFAGEHGHBGHH3?FBDBEGHG?")
        self.trimmer_obj.r2_info = (
            b"@M01750:502:000000000-AUK2G:1:1101:15325:1360 2:N:0:2",
            b"AAGTATAGGTGGATTGGAGTCCTCAGACAACTGTTCAAACTGATGGGACCCACTCCATCGAGATTTCACTGTAGCTAGACCAAAATCACCTATTTTTACTATGAGGTCTTCATGAAGAAATATATCTGAGGAGTAGTAAGTAAAGGAAAA",
            b"1>A1>D@FFB1>GGGGGGGGGGHHH1110BG1AA33DGGFEHH211EFGG?GG00B1F1AF//AEBGHGFF22D2AGHGGHGHC?F1110F1B122F/DFFFGEEGHHH1FHGHEBG1>BF1@BFEDF>F/1?BEFD2BFE2G2F1BBGE")


class TestQiaSeqTrimDuplex(unittest.TestCase):
    '''
    '''
    def setUp(self):
        '''
        '''
        # primer datastruct
        self.primer_datastruct = PrimerDataStruct(k=8,primer_file=os.path.join(os.path.dirname(__file__),"test_data/test_primers_dna.txt"),
                                                  ncpu=1,seqtype="dna",primer_col=3).primer_search_datastruct

        # argparse obj
        self.args = helper.helper_return_args()
        self.args.is_duplex          = True
        self.args.is_phased_adapters = True
        # trimmer obj
        self.trimmer_obj = helper.helper_return_qiaseq_obj(self.args)

    def testDuplexAdapterOffset(self):
        ''' Test whether the duplex adapters have the correct offsets stored
        '''
        for d in self.trimmer_obj._duplex_adapters:
            start, end = d.offset_tag
            self.assertEqual(d.seq[start:end], b"AYYA")

    def testDuplexAdapterTrim(self):
        ''' Test whether the duplex adapters are being correctly trimmed
        '''
        r2_seq = b"TACGAAATAGTACTGAGCGATTATAGGATGGCTGTCAATAATCCCCCGCCTCTGCTGGGCCCTGCGAATCACTCCCTGCCATTGATTACTGAGGAGTGTCAATTTCAGGTTGAATTCATCCCTAGTGAACAAAACGCAATACTGTA"
        self.assertEqual(self.trimmer_obj._id_duplex_tag(r2_seq), (16, b"TT", b"69-8-6")) # perfect match to an adapter

        r2_seq = b"TACGAAATAGTACTGAGCGCCTTTAGGATGGCTGTCAATAATCCCCCGCCTCTGCTGGGCCCTGCGAATCACTCCCTGCCATTGATTACTGAGGAGTGTCAATTTCAGGTTGAATTCATCCCTAGTGAACAAAACGCAATACTGTA"
        self.assertEqual(self.trimmer_obj._id_duplex_tag(r2_seq), (16, b"NN", b"69-8-6")) # imperfect match , ambigous duplex tag

        r2_seq = b"TACGAAATAGTACTGACGCCCATAGGATGGCTGTCAATAATCCCCCGCCTCTGCTGGGCCCTGCGAATCACTCCCTGCCATTGATTACTGAGGAGTGTCAATTTCAGGTTGAATTCATCCCTAGTGAACAAAACGCAATACTGTA"
        self.assertEqual(self.trimmer_obj._id_duplex_tag(r2_seq), (15, b"CC", b"69-8-6")) # imperfect match - 1 deleted base in adapter

        r2_seq = b"TACGAATAGTACTGAGCGCCCATAGGATGGCTGTCAATAATCCCCCGCCTCTGCTGGGCCCTGCGAATCACTCCCTGCCATTGATTACTGAGGAGTGTCAATTTCAGGTTGAATTCATCCCTAGTGAACAAAACGCAATACTGTA"
        self.assertEqual(self.trimmer_obj._id_duplex_tag(r2_seq), (15, b"CC", b"69-7-6")) # imperfect match - 1 deleted base in UMI , changes duplex adapter identfied

        r2_seq = b"TACGAAATAGTACTGAGCGACTATAGGATGGCTGTCAATAATCCCCCGCCTCTGCTGGGCCCTGCGAATCACTCCCTGCCATTGATTACTGAGGAGTGTCAATTTCAGGTTGAATTCATCCCTAGTGAACAAAACGCAATACTGTA"
        self.assertEqual(self.trimmer_obj._id_duplex_tag(r2_seq), (16, b"NN", b"69-8-6")) # perfect match - ambigous duplex tag

                              #GCGAYYATAGGA
        r2_seq = b"TACGAAATAGTAGCGACCATAGGATGGCTGTCAATAATCCCCCGCCTCTGCTGGGCCCTGCGAATCACTCCCTGCCATTGATTACTGAGGAGTGTCAATTTCAGGTTGAATTCATCCCTAGTGAACAAAACGCAATACTGTA"
        self.assertEqual(self.trimmer_obj._id_duplex_tag(r2_seq), (12, b"CC", b"69-4-6")) # perfect match - 69-4-6 adapter

        r2_seq = b"TACGAAATAGTACGACCATAGGATGGCTGTCAATAATCCCCCGCCTCTGCTGGGCCCTGCGAATCACTCCCTGCCATTGATTACTGAGGAGTGTCAATTTCAGGTTGAATTCATCCCTAGTGAACAAAACGCAATACTGTA"
        self.assertEqual(self.trimmer_obj._id_duplex_tag(r2_seq), (11, b"CC", b"69-3-6")) # perfect match - 69-3-6 adapter

        r2_seq = b"TACGAAATAGTACGACCATAGGATGGCTGTCAATAATCCCCCGCCTCTGCTGGGCCCTGCGAATCACTCCCTGCCATTGATTACTGAGGAGTGTCAATTTCAGGTTGAATTCATCCCTAGTGAACAAAACGCAATACTGTA"
        self.assertEqual(self.trimmer_obj._id_duplex_tag(r2_seq), (11, b"CC", b"69-3-6")) # 1-base deletion turns into 69-4-6 into 69-3-6 smaller adapter - which is fine.

if __name__ == '__main__':
    unittest.main()
        
