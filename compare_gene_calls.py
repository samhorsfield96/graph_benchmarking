from Bio import SeqIO


def compare_fasta(ref_fasta, query_fasta, type, group):
    total_ref_records = 0
    total_query_records = 0
    total_correct_query_records = 0

    ref_seq_list = []

    #set group for detection of exact gene call matches
    if group == "1":
        group_list = ["CR931645.1", "CR931648.1", "CR931646.1", "CR931647.1"]
    elif group == "2":
        group_list = ["CR931719.1", "CR931658.1", "CR931659.1", "CR931660.1", "CR931717.1"]
    elif group == "3":
        group_list = ["CR931662.1", "CR931663.1", "CR931664.1", "CR931665.1", "CR931666.1"]

    #create list of all reference genes
    for ref_rec in SeqIO.parse(ref_fasta, "fasta"):
        source_id = ((ref_rec.description.strip()).split("_"))[1]
        total_ref_records += 1
        ref_seq_list.append((source_id, str(ref_rec.seq)))

    unmatched_seq_list = ref_seq_list[:]

    #iterate through query, checking if entry is in ref_seq_list, depending on whether running prodigal or ggCaller due to gene merging in ggCaller
    if type == 'ggc':
        for ggc_rec in SeqIO.parse(query_fasta, "fasta"):
            ORF_colours = (((((ggc_rec.description.strip()).split("["))[4]).replace("]", "")).replace("'", "")).replace(", ", "")
            ORF_colours = list(ORF_colours)
            for i in range(0, len(ORF_colours)):
                if ORF_colours[i] == "1":
                    ORF_pair = (group_list[i], str(ggc_rec.seq))
                    total_query_records += 1
                    if ORF_pair in unmatched_seq_list:
                        total_correct_query_records += 1
                        unmatched_seq_list.remove(ORF_pair)
    elif type == 'prod':
        for prod_rec in SeqIO.parse(query_fasta, "fasta"):
            ORF_id = ((prod_rec.description.strip()).split("_"))[0] + ".1"
            ORF_pair = (ORF_id, str(prod_rec.seq))
            total_query_records += 1
            if ORF_pair in unmatched_seq_list:
                total_correct_query_records += 1
                unmatched_seq_list.remove(ORF_pair)



    #calculate recall and precision
    recall = total_correct_query_records / total_ref_records
    precision = total_correct_query_records / total_query_records

    print(query_fasta)
    print("Total ORFs: {}".format(total_query_records))
    print("Recall: {}".format(recall))
    print("Precision: {}".format(precision))
    return(unmatched_seq_list)


if __name__ == '__main__':
    from Bio import SeqIO
    import sys

    reference_genes = sys.argv[2]
    gene_calls = sys.argv[2]
    caller_type = sys.argv[3]
    group = sys.argv[4]

    unmatched = compare_fasta(reference_genes, gene_calls, caller_type, group)
