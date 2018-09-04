# read-trimmer
Paired end trimmer for QIASeq DNA/RNA reads 

## Installation
```
pip3 install edlib
pip3 install Cython
git clone https://github.com/qiaseq/read-trimmer.git
cd read-trimmer
python3 setup.py build_ext --inplace
```

## Usage Instructions
```
python3 read-trimmer/trimmer/run.py --help

Example command :
python3 trimmer/run.py --r1 <R1 Fastq Path>  --r2 <R2 Fastq Path> --out_r1 <Trimmed R1 Fastq Path> --out_r2 <Trimmed R2 Fastq Path> --out_metrics <Output Metric File Path> --primer_file tests/test_data/test_primers_dna.txt --seqtype dna --primer_col 3 --umi_len 12 --common_seq_len 11

Run tests :
python3 trimmer/tests/test_qiaseq_trim.py -v
```

