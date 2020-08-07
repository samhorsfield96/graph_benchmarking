def analyse_read_coverage(gaf, blast, reads, outfile):
    import numpy as np
    from statistics import mean

    print("Generating gene dictionary from blast file")

    #caculate reads generated for each gene
    gene_dict = {}
    with open(blast, "r") as b:
        for line in b:
            line_list = line.split("\t")
            gene_start = int(line_list[8])
            gene_end = int(line_list[9])
            line_list[0] = line_list[0].replace(">", "")
            line_list[0] = line_list[0].replace(",", "")
            gene_dict[line_list[0]] = {}
            gene_dict[line_list[0]]["reads"] = {}
            gene_dict[line_list[0]]["range"] = (gene_start, gene_end)
            #gene_dict[line_list[0]]["read_coverage"] = set()

    print("Generating total read dictionary from read list")

    total_read_dict = {}
    with open(reads, "r") as r:
        for line in r:
            read_id = line.strip()
            split_line = read_id.split("_")
            read_start = int(split_line[1])
            read_end = int(split_line[1]) + int(split_line[6])
            read_id = read_id.replace(">", "")
            total_read_dict[read_id] = (read_start, read_end)

    print("Adding read mappings to gene dictionary")

    #total_genic_content_reads = 0

    #work out which reads map to which genes
    for gene, gene_item in gene_dict.items():
        start_gene, end_gene = gene_item["range"]
        for read, read_item in total_read_dict.items():
            start_read, end_read = read_item
            if (start_gene <= start_read <= end_gene) or (start_gene <= end_read <= end_gene) or (start_read <= start_gene and end_read >= end_gene):
                gene_dict[gene]["reads"][read] = {}
                gene_dict[gene]["reads"][read]["num_mappings"] = 0
                gene_dict[gene]["reads"][read]["read_length"] = 0
                gene_dict[gene]["reads"][read]["read_covered"] = 0
                gene_dict[gene]["reads"][read]["range"] = read_item
                #gene_dict[gene]["read_coverage"].update(range(start_read, end_read))
        #total_genic_content_reads += len(gene_dict[gene]["read_coverage"])

    #print(total_genic_content_reads)

    print("Analysing read mappings from gaf file")

    read_dict = {}
    with open(gaf, "r") as f:
        for line in f:
            map_list = line.split("\t")
            read_name = map_list[0]
            read_length = int(map_list[1])
            read_start = int(map_list[2])
            read_end = int(map_list[3])
            read_covered = range(read_start, read_end)

            if map_list[0] not in read_dict:
                read_dict[read_name] = {}
                read_dict[read_name]['num_mappings'] = 1
                read_dict[read_name]['read_length'] = read_length
                read_dict[read_name]['read_covered'] = set(read_covered)
            else:
                read_dict[read_name]['num_mappings'] += 1
                read_dict[read_name]['read_covered'].update(read_covered)

    print(read_dict)

    print("Adding read mappings to gene dictionary")

    for read, read_item in read_dict.items():
        for gene, gene_item in gene_dict.items():
            try:
                gene_dict[gene]["reads"][read]["num_mappings"] = read_dict[read]['num_mappings']
                gene_dict[gene]["reads"][read]["read_length"] = read_dict[read]['read_length']
                gene_dict[gene]["reads"][read]["read_covered"] = read_dict[read]['read_covered']
            except KeyError:
                pass

    print(gene_dict)

    print("Calculating average per gene read mappings")

    perc_unique_gene_mapping = []
    perc_multi_gene_mapping = []
    perc_total_gene_mapping = []
    num_unique_gene_mapping = []
    num_multi_gene_mapping = []
    num_total_gene_mapping = []

    mapped_genic_read_bases_dict = {}
    unique_coverage = {}
    multi_coverage = {}
    total_genic_read_bases_dict = {}
    for gene, gene_item in gene_dict.items():
        read_count = 0
        unique_reads = 0
        multi_reads = 0
        for read, read_item in gene_item["reads"].items():
            if read_item['read_covered'] != 0:
                mapped_genic_read_bases_dict[read] = len(read_item['read_covered'])
            total_genic_read_bases_dict[read] = read_item['read_length']
            read_count += 1
            if read_item["num_mappings"] == 1:
                unique_reads += 1
                unique_coverage[read] = (len(read_item['read_covered']) / read_item['read_length']) * 100
            elif read_item["num_mappings"] > 1:
                multi_reads += 1
                multi_coverage[read] = (len(read_item['read_covered']) / read_item['read_length']) * 100
            else:
                pass
        if read_count != 0:
            perc_unique = (unique_reads / read_count) * 100
            perc_multi = (multi_reads / read_count) * 100
            perc_total_mapped = ((unique_reads + multi_reads) / read_count) * 100
        else:
            perc_unique = 0
            perc_multi = 0
            perc_total_mapped = 0

        perc_unique_gene_mapping.append(perc_unique)
        perc_multi_gene_mapping.append(perc_multi)
        perc_total_gene_mapping.append(perc_total_mapped)
        num_unique_gene_mapping.append(unique_reads)
        num_multi_gene_mapping.append(multi_reads)
        num_total_gene_mapping.append((unique_reads + multi_reads))

    print("Calculating genic-read mapping statistics")
    print(mapped_genic_read_bases_dict)

    genic_mapping_reads_bases = 0
    for key, item in mapped_genic_read_bases_dict.items():
        genic_mapping_reads_bases += item

    total_genic_read_bases = 0
    for key, item in total_genic_read_bases_dict.items():
        total_genic_read_bases += item

    unique_coverage_list = []
    for key, item in unique_coverage.items():
        unique_coverage_list.append(item)

    multi_coverage_list = []
    for key, item in multi_coverage.items():
        multi_coverage_list.append(item)

    total_coverage_list = unique_coverage_list + multi_coverage_list


    print("Generating results dictionary")

    results_header = ["unique_mapping_per_gene", "multi_mapping_per_gene", "total_mapping_per_gene", "unique_mapping_per_gene_perc", "multi_mapping_per_gene_perc", "total_mapping_per_gene_perc", "unique_mapped_genic_read_coverage", "multi_mapped_genic_read_coverage", "total_mapped_genic_read_coverage"]
    results_lists = [num_unique_gene_mapping, num_multi_gene_mapping, num_total_gene_mapping, perc_unique_gene_mapping, perc_multi_gene_mapping, perc_total_gene_mapping, unique_coverage_list, multi_coverage_list, total_coverage_list]
    results_dict = {}

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

    print("Generating output file")

    with open(outfile, "w") as o:
        o.write("genic_mapping_reads_bases" + "\t" + str(genic_mapping_reads_bases) + "\n" + \
                "total_genic_read_bases" + "\t" + str(total_genic_read_bases) + "\n")
        for key, item in results_dict.items():
            for result, stat in item.items():
                o.write(str(key) + "_" + str(result) + "\t" + str(stat) + "\n")

if __name__ == '__main__':
    import sys
    gaf = sys.argv[1]
    blast = sys.argv[2]
    reads = sys.argv[3]
    outfile = sys.argv[4]
    analyse_read_coverage(gaf, blast, reads, outfile)