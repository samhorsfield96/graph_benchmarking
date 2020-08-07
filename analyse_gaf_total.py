def analyse_gaf_total(infile, outfile):
    import re
    from statistics import mean
    import numpy as np

    total_map_count = 0
    read_dict = {}
    total_node_set = set()
    with open(infile, "r") as f:
        for line in f:
            map_list = line.split("\t")
            read_name = map_list[0]
            map_identity = float(re.split(':', map_list[14])[2])
            align_block_len = int(map_list[10])
            read_length = int(map_list[1])
            read_start = int(map_list[2])
            read_end = int(map_list[3])
            read_covered = range(read_start, read_end)

            total_map_count += 1

            nodes = re.split('[>|<]', map_list[5])
            nodes = [x for x in nodes if x]

            total_node_set.update(tuple(nodes))

            if map_list[0] not in read_dict:
                read_dict[read_name] = {}
                read_dict[read_name]['num_mappings'] = 1
                read_dict[read_name]['nodes'] = set(nodes)
                read_dict[read_name]['read_length'] = read_length
                read_dict[read_name]['read_covered'] = set(read_covered)
                read_dict[read_name]['identity'] = []
                read_dict[read_name]['identity'].append(map_identity)
                read_dict[read_name]['align_seq'] = []
                read_dict[read_name]['align_seq'].append(int(align_block_len))
            else:
                read_dict[read_name]['num_mappings'] += 1
                read_dict[read_name]['nodes'].update(tuple(nodes))
                read_dict[read_name]['read_covered'].update(read_covered)
                read_dict[read_name]['identity'].append(map_identity)
                read_dict[read_name]['align_seq'].append(int(align_block_len))

        #mapping totals
        total_unique_reads = 0
        total_multi_reads = 0

        #identity averages
        unique_identity = []
        multi_identity = []

        #read % coverage
        unique_coverage = []
        multi_coverage = []
        total_unique_read_bases = 0

        #read mapping average
        read_mapping = []

        #alignment averages
        unique_seq_alignment = []
        multi_seq_alignment = []

        #node mappings
        unique_node_mappings = 0
        multi_node_mappings = 0

        for key, items in read_dict.items():
            read_mapping.append(items['num_mappings'])
            total_unique_read_bases += len(items['read_covered'])
            if items['num_mappings'] == 1:
                total_unique_reads += 1
                unique_identity.extend(items['identity'])
                unique_seq_alignment.extend(items['align_seq'])
                unique_node_mappings += len(items['nodes'])
                unique_coverage.append((len(items['read_covered'])/items['read_length']) * 100)
            else:
                total_multi_reads += 1
                multi_identity.extend(items['identity'])
                multi_seq_alignment.extend(items['align_seq'])
                multi_node_mappings += len(items['nodes'])
                multi_coverage.append((len(items['read_covered']) / items['read_length']) * 100)

        # generate results dict for stats results

        total_identity = unique_identity + multi_identity
        total_seq_alignment = unique_seq_alignment + multi_seq_alignment
        total_coverage = unique_coverage + multi_coverage

        results_header = ["unique_identity", "multi_identity", "total_identity", "unique_alignment", "multi_alignment", "total_alignment", "unique_read_coverage", "multi_read_coverage", "total_read_coverage"]
        results_lists = [unique_identity, multi_identity, total_identity, unique_seq_alignment, multi_seq_alignment, total_seq_alignment, unique_coverage, multi_coverage, total_coverage]
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

        #non-stats results
        unique_sum_alignment = sum(unique_seq_alignment)
        multi_sum_alignment = sum(multi_seq_alignment)
        total_sum_alignment = unique_sum_alignment + multi_sum_alignment
        total_node_mappings = len(total_node_set)
        average_read_mapping = mean(read_mapping)
        total_reads = total_unique_reads + total_multi_reads


        with open(outfile, "w+") as w:
            w.write("Total_map_count" + "\t" + str(total_map_count) + "\n" + \
                    "Total_reads" + "\t" + str(total_reads) + "\n" + \
                    "Unique_read_map_count" + "\t" + str(total_unique_reads) + "\n" + \
                    "Multi_read_map_count" + "\t" + str(total_multi_reads) + "\n" + \
                    "Average_maps_per_read" + "\t" + str(average_read_mapping) + "\n" + \
                    "Total_align_len" + "\t" + str(total_sum_alignment) + "\n" + \
                    "Unique_align_len" + "\t" + str(unique_sum_alignment) + "\n" + \
                    "Multi_align_len" + "\t" + str(multi_sum_alignment) + "\n")
            for key, item in results_dict.items():
                for result, stat in item.items():
                    w.write(str(key) + "_" + str(result) + "\t" + str(stat) + "\n")
            w.write("Total_nodes_mapped" + "\t" + str(total_node_mappings) + "\n" + \
                    "Unique_node_mappings" + "\t" + str(unique_node_mappings) + "\n" + \
                    "Multi_node_mappings" + "\t" + str(multi_node_mappings) + "\n" + \
                    "Total_unique_read_bases" + "\t" + str(total_unique_read_bases))

if __name__ == '__main__':
    import sys
    gaf = sys.argv[1]
    outfile = sys.argv[2]
    analyse_gaf_total(gaf, outfile)
