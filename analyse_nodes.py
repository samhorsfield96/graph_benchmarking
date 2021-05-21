def flip_sign(string):
    if string == '-':
        return '+'
    elif string == '+':
        return '-'
    else:
        pass

def analyse_nodes(infile, outpref):
    import numpy as np
    from statistics import mean

    stop_codons = ["TAA", "TGA", "TAG", "TTA", "TCA", "CTA"]
    start_codons = ["ATG", "GTG", "TTG", "CAT", "CAC", "CAA"]

    node_dict = {}
    link_dict = {}
    number_nodes = 0
    number_edges = 0
    number_stops = 0
    number_starts = 0

    with open(infile, "r") as f:
        for line in f:
            split_line = (line.strip()).split("\t")
            entry_type = split_line[0]

            #count nodes
            if entry_type == 'S':
                number_nodes += 1
                node_id = split_line[1]
                node_len = len(split_line[2])

                # determine number of stop codons
                number_stops = number_stops + sum([split_line[2].count(x) for x in stop_codons])
                number_starts = number_starts + sum([split_line[2].count(x) for x in start_codons])

                if node_id not in node_dict:
                    node_dict[node_id] = {}
                    node_dict[node_id]["node_len"] = node_len
                    node_dict[node_id]["node_deg"] = 0
                else:
                    node_dict[node_id]["node_len"] = node_len

            #count non-redundant links
            elif entry_type == 'L':
                node_source = split_line[1]
                node_source_dir = split_line[2]
                node_sink = split_line[3]
                node_sink_dir = split_line[4]

                #flip signs
                flip_source = flip_sign(node_source_dir)
                flip_sink = flip_sign(node_sink_dir)

                #generate link ids
                link_id = str(node_source) + str(node_source_dir) + str(node_sink) + str(node_sink_dir)
                flip_link_id = str(node_sink) + str(flip_sink) + str(node_source) + str(flip_source)

                if link_id not in link_dict and flip_link_id not in link_dict:
                    number_edges += 1
                    link_dict[link_id] = 1
                    link_dict[flip_link_id] = 1

                    if node_source not in node_dict:
                        node_dict[node_source] = {}
                        node_dict[node_source]["node_len"] = 0
                        node_dict[node_source]["node_deg"] = 1
                    else:
                        node_dict[node_source]["node_deg"] += 1

                    if node_sink not in node_dict:
                        node_dict[node_sink] = {}
                        node_dict[node_sink]["node_len"] = 0
                        node_dict[node_sink]["node_deg"] = 1
                    else:
                        node_dict[node_sink]["node_deg"] += 1

                else:
                    pass

            else:
                pass


    #output distributions of node length and degree
    dist_file = outpref + "_dists.txt"
    with open(dist_file, "w") as d:
        d.write("Node_id\tNode_len\tNode_deg\n")
        for key, item in node_dict.items():
            d.write(str(key) + "\t" + str(item["node_len"]) + "\t" + str(item["node_deg"]) + "\n")


    #create distributions
    node_len_list = []
    node_deg_list = []

    for key, item in node_dict.items():
        node_len_list.append(item["node_len"])
        node_deg_list.append(item["node_deg"])

    #create summary results dictionary
    results_header = ["node_len", "node_deg"]
    results_lists = [node_len_list, node_deg_list]
    results_dict = {}

    total_bases = sum(node_len_list)

    for index, item in enumerate(results_header):
        results_dict[item] = {}
        results_dict[item]["mean"] = mean(results_lists[index])
        results_dict[item]["std"] = np.std(results_lists[index])
        results_dict[item]["min"] = min(results_lists[index])
        results_dict[item]["max"] = max(results_lists[index])
        q75, q50, q25 = np.percentile(results_lists[index], [75, 50, 25])
        results_dict[item]["uq"] = q75
        results_dict[item]["med"] = q50
        results_dict[item]["lq"] = q25
        results_dict[item]["iqr"] = (q75 - q25)

    #output summary file
    summmary_file = outpref + "_summary.txt"
    with open(summmary_file, "w") as s:
        s.write("total_nodes" + "\t" + str(number_nodes) + "\n" + "total_directed_edges" + "\t" + str(number_edges) + "\n" +
        "total_bases" + "\t" + str(total_bases) + "\n" + "total_stop_codons" + "\t" + str(number_stops) +
        "\n" + "total_start_codons" + "\t" + str(number_starts) + "\n")
        for key, item in results_dict.items():
            for result, stat in item.items():
                s.write(str(key) + "_" + str(result) + "\t" + str(stat) + "\n")

if __name__ == '__main__':
    import sys
    gfa = sys.argv[1]
    outpref = sys.argv[2]
    analyse_nodes(gfa, outpref)
