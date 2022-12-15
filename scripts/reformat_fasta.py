import shutil
from io import TextIOWrapper

COUNT = 0


def parse_fasta(f: TextIOWrapper) -> list:
    desc = ""
    seq = ""
    for l in f:
        l = l.strip()
        if l[0] == ">":
            yield desc, seq
            desc = l
            seq = ""
        else:
            seq += l
    yield desc, seq


def filter_seqs(seqs: list, min: int) -> list:
    for desc, seq in seqs:
        if len(seq) >= min:
            global COUNT
            COUNT += 1
            yield "> c_%012d" % COUNT, seq


def write_fasta(f: TextIOWrapper, seqs: list):
    # seqs.sort(key=lambda t: -len(t[1])) # Forces everything into mem at once but doesn't seem like there's a good way around that
    for desc, seq in seqs:
        f.write(f"{desc}\n{seq}\n")


with open(snakemake.input[0]) as f:
    with open(f"{snakemake.output[0]}", "w") as g:
        write_fasta(g, list(filter_seqs(parse_fasta(f), snakemake.params.len)))
