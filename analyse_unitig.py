def analyse_unitigs(gfa_stats, rtab, output):
    unitig_dict = {}

    rtab_set = set()
    stats_set = set()
    with open(gfa_stats, "r") as g, open(rtab, "r") as r:
        g.readline() #header
        for line in g:
            split_line = (line.strip()).split("\t")
            node_id = split_line[0]
            node_len = split_line[1]
            node_deg = split_line[2]

            unitig_dict[node_id] = {}
            unitig_dict[node_id]['node_len'] = node_len
            unitig_dict[node_id]['node_deg'] = node_deg
            unitig_dict[node_id]['allele_freq'] = 0

            stats_set.add(node_id)

        r.readline() #header
        for line in r:
            split_line = (line.strip()).split("\t")
            node_id = split_line[0]
            freq_list = split_line[1:]
            freq_list = [int(i) for i in freq_list]
            allele_freq = (sum(freq_list)/len(freq_list))
            unitig_dict[node_id]['allele_freq'] = allele_freq

            rtab_set.add(node_id)

    set_diff = stats_set.difference(rtab_set)
    print(len(stats_set))
    print(len(rtab_set))
    print(set_diff)
    print(len(set_diff))

    # output new distributions of node length, degree and allele frequency
    with open(output, "w") as d:
        d.write("Node_id\tNode_len\tNode_deg\tAllele_freq\n")
        for key, item in unitig_dict.items():
            d.write(str(key) + "\t" + str(item["node_len"]) + "\t" + str(item["node_deg"]) + "\t" + str(item["allele_freq"]) + "\n")

if __name__ == '__main__':
    import sys
    gfa_stats = sys.argv[1]
    rtab = sys.argv[2]
    output = sys.argv[3]
    analyse_unitigs(gfa_stats, rtab, output)