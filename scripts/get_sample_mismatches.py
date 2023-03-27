from Bio import pairwise2, SeqIO


def create_dict(fasta):
    gene_names_dict = {}
    seq_dict = {}
    for keys in fasta:
        gene = keys.split("___")[0]
        length = keys.split("length:")[1].rstrip()
        gene_samp = keys.split("bin_id:")[1]
        gene_samp = gene_samp.split("|")[0]
        gene_samp = gene_samp.split("_contigs")[0]
        gene_samp = gene + ", " + gene_samp + "__" + length
        gene_names_dict[gene_samp] = gene
        seq_dict[gene_samp] = str(fasta[keys])
    return gene_names_dict, seq_dict


def format_mismatch(align1, align2, score, begin, end):
    mismatch = 0
    for a, b in zip(align1[begin:end], align2[begin:end]):
        if a != b and not (a == "-" or b == "-"):
            mismatch += 1
    return mismatch


all_fasta_sequences = {
    fa.description: fa.seq for fa in SeqIO.parse(snakemake.input.all, "fasta")
}
query_fasta_sequences = {
    fa.description: fa.seq for fa in SeqIO.parse(snakemake.input.sample, "fasta")
}
all_genes, all_seqs = create_dict(all_fasta_sequences)
query_genes, query_seqs = create_dict(query_fasta_sequences)
with open(snakemake.output[0], "w") as out:
    out.write("query\tsubject\tmismatches\n")
    for keys in query_genes:
        query_seq = query_seqs[keys]
        subj_names = [
            names
            for names, gene_name in all_genes.items()
            if gene_name == query_genes[keys]
        ]
        for subj in subj_names:
            subj_seq = all_seqs[subj]
            mismatches = format_mismatch(
                *pairwise2.align.globalms(
                    query_seq,
                    subj_seq,
                    5,
                    -4,
                    -10,
                    -4,
                    penalize_extend_when_opening=True,
                    penalize_end_gaps=False,
                    one_alignment_only=True,
                )[0]
            )
            out.write(keys + "\t" + subj + "\t" + str(mismatches) + "\n")
