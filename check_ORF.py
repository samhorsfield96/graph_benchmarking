def check_ORF_in_ref(ref_fasta_for, ref_fasta_rev, query_fasta, outfasta):
    from Bio import SeqIO
    ref_seq_list = []
    for ref_rec in SeqIO.parse(ref_fasta_for, "fasta"):
        ref_seq_list.append(str(ref_rec.seq))
    for ref_rec in SeqIO.parse(ref_fasta_rev, "fasta"):
        ref_seq_list.append(str(ref_rec.seq))
    with open(outfasta, "w") as f:
        for query_rec in SeqIO.parse(query_fasta, "fasta"):
            if not any(str(query_rec.seq) in ref_seq for ref_seq in ref_seq_list):
                f.write(">" + str(query_rec.description) + "\n" + str(query_rec.seq) + "\n")

def check_ref_in_query(ref_fasta, query_fasta, outfile):
    from Bio import SeqIO
    query_list = []
    for query_rec in SeqIO.parse(query_fasta, "fasta"):
        query_list.append(str(query_rec.seq))

    with open(outfile, "w") as f:
        for ref_rec in SeqIO.parse(ref_fasta, "fasta"):
            if not any(str(ref_rec.seq) in query_rec for query_rec in query_list):
                f.write(">" + str(ref_rec.description) + "\n" + str(ref_rec.seq) + "\n")
