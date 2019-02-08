import argparse
import os
import multiprocessing

import qiaseq_trimmer

def init_parser():
    '''
    '''
    global parser
    parser = argparse.ArgumentParser(description = "A customizable adapter and primer trimming tool for QIASeq reads which relies on paired end overlap.")
    
    parser.add_argument("--r1", required = True,help = "Input R1 fastq file")
    
    parser.add_argument("--r2", required = True,help = "Input R2 fastq file")
    
    parser.add_argument("--out_r1", required = True, help = "Output file for R1 trimmed fastq")
    
    parser.add_argument("--out_r2", required = True,
                        help = "Output file for R2 trimmed fastq")
    
    parser.add_argument("--out_metrics", required = True,
                        help = "Output file path for trimming metrics")
    
    parser.add_argument("--primer_file", required = True,
                        help = "Primer file with 3' coordinates" )
    
    parser.add_argument("--primer_col", required = True, type = int,
                        help = "0-based Column number with primer sequence in the Primer file")
    
    parser.add_argument("--seqtype", default = "dna", const = "dna",
                        nargs = "?", choices = ["dna","rna"],
                        help = "Sequencing type : dna/rna. Default : %(default)s")
    
    parser.add_argument("--is_nextseq", action = "store_true",
                        help = "Whether this is a NextSeq sequencing run")
    
    parser.add_argument("--is_duplex", action = "store_true",
                         help = "Whether this is a duplex sequencing experiment")

    parser.add_argument("--is_r2_primer_side",action = "store_true", help = "Is R2 read the primer side?")    
    
    parser.add_argument("--check_primer_side", action = "store_true",
                        help = "User primer side overlap coordinates for trimming R2 3' end")
    
    parser.add_argument("--custom_seq_adapter", default = "AATGTACAGTATTGCGTTTTG",
                        help = "The custom sequencing adapter used in library preperation (SPE side). Default : %(default)s")
    
    parser.add_argument("--trim_custom_seq_adapter", action = "store_true",
                        help = "Choose this flag to trim custom sequencing adapter."\
                        "If not selected the first ~ 40000 reads will be checked heuristically to determine whether to trim or not.")
    
    parser.add_argument("--tagname_primer", default = "pr",
                        help = "Tag name for Primer ID. Default : %(default)s")
    
    parser.add_argument("--tagname_primer_error", default = "pe",
                        help = "Tag name for Primer Edit Dist. Default : %(default)s")
    
    parser.add_argument("--tagname_umi", default = "mi",
                        help = "Tag name for UMI sequence. Default : %(default)s")
    
    parser.add_argument("--tagname_duplex", default = "DU",
                        help = "Tag name for duplex tag. Default : %(default)s")
    
    parser.add_argument("--tag_separator",default = "\t",
                        help = "separator for readID, umi, and primer tags. Default : %(default)s")
    
    parser.add_argument("--field_separator",default = "\n",
                        help = "separator for readID, fastq sequence, and fastq base quality. Default : %(default)s")

    parser.add_argument("--no_tagnames", action = "store_true",
                        help = "Choose this option to have no tagnames.")
    
    parser.add_argument("--primer3_bases_R1", required = False, default = 8,
                        type = int, help = "Number of bases on 3' end of primer to keep on R1. Default : %(default)s")
    
    parser.add_argument("--primer3_bases_R2", required = False, default = 8,
                        type = int, help = "Number of bases on 3' end of primer to keep on R2. Default : %(default)s")
    
    parser.add_argument("--umi_len", required = True, type = int, help = "Length of UMI sequence")

    parser.add_argument("--common_seq_len", required = True, type = int, help = "Length of the common sequence (UMI side).")

    parser.add_argument("--min_primer_side_len", default = 50, type = int,
                        help = "Minimum length of the Primer side read. Default : %(default)s")
    
    parser.add_argument("--min_umi_side_len", default = 50, type = int,
                        help = "Minimum length of the UMI side read. Default : %(default)s")
    
    parser.add_argument("--overlap_check_len", default = 25,
                        type = int, help = "Sequence length for overlap check. Default : %(default)s")
    
    parser.add_argument("--ncpu", default = multiprocessing.cpu_count(),
                        type = int, help = "Number of CPUs to use. Default : %(default)s")
    
    parser.add_argument("--max_mismatch_rate_overlap", default = 0.12,
                        type = float, help = "Mismatch rate to tolerate for overlap check. Default : %(default)s")
    
    parser.add_argument("--max_mismatch_rate_primer", default = 0.12,
                        type = float, help = "Mismatch rate to tolerate for primer identification. Default : %(default)s")

    parser.add_argument("--poly_tail_primer_side", default = "none", const = "none",
                        nargs = "?", choices = ["polyA","polyT","none"],
                        help = "Choose whether there is a polyA/T tail on the primer side. Default : %(default)s")

    parser.add_argument("--poly_tail_umi_side", default = "none", const = "none",
                        nargs = "?", choices = ["polyA","polyT","none"],
                        help = "Choose whether there is a polyA/T tail on the UMI side. Default : %(default)s")    

    parser.add_argument("--umi_filter_min_bq", default = 20, type = int,
                        help = "Minimum Base quality below which a base in the UMI sequence is considered to be of low quality."\
                        "Only applicable to speRNA. Default : %(default)s")

    parser.add_argument("--umi_filter_max_lowQ_bases", default = 1, type = int,
                        help = "Maximum number of lowQ bases to tolerate in the UMI region. Reads having more than this number are dropped."\
                        "Only applicable to speRNA. Default : %(default)s")

    parser.add_argument("--umi_filter_max_Ns", default = 1, type = int, help = "Tolerate these many Ns in the UMI sequence."\
                        "Reads having more than this number are dropped. Only applicable to speRNA. Default : %(default)s")

    parser.add_argument("--thousand_comma", action = "store_true", help = "Show integers with thousand's comma. Default : %(default)s")

    parser.add_argument("--is_multimodal", action = "store_true", help = "Reads generated using multimodal kit.")

    parser.add_argument("--tagname_multimodal", default = "MM",
                        help = "Tag name for multimodal common sequence tag. Default : %(default)s")

    parser.add_argument("--umi_len_alt", type = int, 
                        help = "Length of alternative UMI sequence, e.g. DNA UMI length for RNA reads and vice versa, only applicable to multimodal reads.")

def main(args):
    '''
    '''
    qiaseq_trimmer.main(args)
    
if __name__ == "__main__":
    init_parser()
    args = parser.parse_args()
    main(args)

