from collections import Counter
from defaults import ProgDefaults as dv


# from https://stackoverflow.com/questions/312443/how-do-i-split-a-list-into-equally-sized-chunks
# authored by Ned Batchelder
def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def parse_file_head_todict(foi, key_idx):
    parse_dict = {}
    cycle = 0
    with open(foi, "r") as f:
        while True:
            # reads line-by-line to reduce memory load
            line = f.readline()
            if not line:
                break
            if line.startswith("#"):
                continue

            line = line.rstrip().split('\t')
            if cycle == 0:
                parse_dict.setdefault('header', line)
            else:
                parse_dict.setdefault(line[key_idx], line)
            cycle += 1
    return parse_dict


def parse_file_nohead_todict(foi, key_idx):
    parse_dict = {}
    cycle = 0
    with open(foi, "r") as f:
        while True:
            # reads line-by-line to reduce memory load
            line = f.readline()
            if not line:
                break
            if line.startswith("#"):
                continue

            line = line.rstrip().split('\t')
            parse_dict.setdefault(line[key_idx], line)
            cycle += 1
    return parse_dict


def parse_file_nohead_tolist(foi):
    parse_list = []
    cycle = 0
    with open(foi, "r") as f:
        while True:
            # reads line-by-line to reduce memory load
            line = f.readline()
            if not line:
                break
            if line.startswith("#"):
                continue

            line = line.rstrip().split('\t')
            parse_list.append(line)
            cycle += 1
    return parse_list


def check_gaps(gene_list):
    accessions = list(set(k[dv.ORDER_DICT["accession"]] for k in gene_list))
    accessions.sort()

    ungapped_entries = []
    orig_gapped = []
    grouped_accessions = {i: [] for i in accessions}
    for gene in gene_list:
        grouped_accessions[gene[dv.ORDER_DICT["accession"]]].append(gene)

    for i in grouped_accessions.keys():
        acc_count = Counter([j[dv.ORDER_DICT["gene"]] for j in grouped_accessions[i]])
        max_num_genes = acc_count.most_common()[0][1]
        if max_num_genes > 1:
            gapped_genes = [j[0] for j in acc_count.most_common() if j[1] == max_num_genes]
            for gene in gapped_genes:
                more_than_two = [j for j in grouped_accessions[i] if j[dv.ORDER_DICT["gene"]] == gene]

                range_nums = []
                for j in more_than_two:
                    for k in j[dv.ORDER_DICT["range"]].split(":"):
                        range_nums.append(int(k))
                range_nums.sort()

                new_entry = more_than_two[0].copy()
                new_entry[dv.ORDER_DICT["range"]] = "%s:%s" % (min(range_nums), max(range_nums))
                new_entry[dv.ORDER_DICT["length"]] = "%s" % (max(range_nums) - min(range_nums))
                print(new_entry)
                ungapped_entries.append(new_entry)

                for j in more_than_two:
                    orig_gapped.append(j)

    for i in orig_gapped:
        gene_list.remove(i)
    for i in ungapped_entries:
        gene_list.append(i)

    gene_list.sort(key=lambda i: i[dv.ORDER_DICT["species"]])


def update_gene_dict_from_file(update_file, ref_dict):
    temp_add_dict = parse_file_nohead_tolist(update_file)
    temp_add_dict = [[b.replace("\"", "") for b in m] for m in temp_add_dict if m != ['']]

    if len(set(m[1] for m in temp_add_dict)) != len(temp_add_dict):
        print("FATAL: Ambiguity in add_dict file!")
        add_dict_dups = Counter([m[1] for m in temp_add_dict])
        for y in [b for b in add_dict_dups.keys() if add_dict_dups[b] > 1]:
            print("\'%s\' appears %s times" % (y, add_dict_dups[y]))
        exit()
    final_add_dict = {}
    for m in temp_add_dict:
        try:
            final_add_dict.setdefault(m[1], ref_dict[m[0]])
        except KeyError:
            final_add_dict.setdefault(m[1], m[0])

    ref_dict.update(final_add_dict)
