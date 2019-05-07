import os
import sqlite3
import csv
from os import listdir
from os.path import isfile, join
from collections import Counter
from collections import defaultdict
from ete3 import NCBITaxa


d = defaultdict(list)
taxids = []
mypath = input("Enter path to csv files (i.e. /home/user/csv/files/): ")
allfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
# print(str(allfiles))
onlyfiles = [s for s in allfiles if '.csv' in s]
filecount = len(onlyfiles)
level = input("1: phylum ----> 2: class ----> 3: order ----> 4: family ----> 5: genus ----> 6: species\nEnter rank number: ")
ncbi = NCBITaxa()
DEFAULT_TAXADB = os.path.join(os.environ.get('HOME', '/'), '.etetoolkit', 'taxa.sqlite')
DB_VERSION = 2
rank_dict = {"1": "phylum", "2": "class", "3": "order", "4": "family", "5": "genus", "6": "species"}


def is_taxadb_up_to_date(dbfile=DEFAULT_TAXADB):
    """
    Check if a valid and up-to-date taxa.sqlite database exists
    If dbfile is not specified, DEFAULT_TAXADB is assumed
    :param dbfile:
    :return:
    """
    db = sqlite3.connect(dbfile)
    try:
        r = db.execute('SELECT version FROM stats;')
        version = r.fetchone()[0]
    except (sqlite3.OperationalError, ValueError, IndexError, TypeError):
        version = None
    db.close()
    if version != DB_VERSION:
        return False
    return True


def get_desired_ranks(taxid, desired_ranks):
    """
    Get rank entered by user for all taxids.
    :param taxid:
    :param desired_ranks:
    :return:
    """
    lineage = ncbi.get_lineage(taxid)
    names = ncbi.get_taxid_translator(lineage)
    lineage2ranks = ncbi.get_rank(names)
    ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())
    return{'{}_id'.format(rank): ranks2lineage.get(rank, '<not present>') for rank in desired_ranks}


def read_csv():
    """
    Read the csv files and build lists with available taxids and occurance.
    :return:
    """
    count = 0
    for csv in onlyfiles:
        klist = []
        taxnames = []
        acc_column = 6
        lca_column = 7
        lca_data = ['0', '1']
        taxid_column = 4
        with open(mypath + csv) as nf:
            linenr = 0
            for line in nf:
                if ',,' not in line:
                    line = line.split(',')
                    name = ""
                    if linenr > 0:
                        if float(line[acc_column]) >= 80 and line[lca_column].strip('\r\n') in lca_data:
                            if rank_dict[level] == "species":
                                taxid2name = ncbi.get_taxid_translator([int(line[taxid_column])])
                                for dummytid, tname in taxid2name.items():
                                    namelist = str(tname).split(" ")[:2]
                                for nl in namelist:
                                    name += str(nl) + " "
                                name += line[lca_column].strip('\r\n')
                                taxnames.append(str(name))
                            else:
                                taxids = [line[taxid_column]]
                                desired_ranks = [rank_dict[level]]
                                for taxid in taxids:
                                    ranks = get_desired_ranks(taxid, desired_ranks)
                                    for dummykey, rank in ranks.items():
                                        if rank != '<not present>':
                                            taxid2name = ncbi.get_taxid_translator([int(rank)])
                                            for dummytid, tname in taxid2name.items():
                                                namelist = str(tname).split(" ")[:2]
                                            for nl in namelist:
                                                name += str(nl) + " "
                                            name += line[lca_column].strip('\r\n')
                                            taxnames.append(str(name))
                linenr += 1
            # nf.close()
            namecounter = Counter(taxnames)
            for k, v in namecounter.items():
                if k not in d and count > 0:
                    for dummyc in range(0, count):
                        d[k].append(0)
                d[k].append(v)
                klist.append(k)
            for dk, dummydv in d.items():
                if dk not in klist:
                    d[dk].append(0)
            count += 1
            namecounter.clear()
    for k, v in d.items():
        if len(d[k]) < len(onlyfiles):
            d[k].append(0)
    make_rank_csv()
    print ("Finished\n")


def make_rank_csv():
    """
    Create a new csv file based on the rank provided by the user.
    The WIMP csv files will be used as input and the rank will be searched using the NCBI taxid translator.
    :return:
    """
    print ("Creating csv file based on rank...\n")
    with open(rank_dict[level] + '.csv', 'w') as tf:
        tf.write(rank_dict[level] + ',')
        i = 0
        for s in onlyfiles:
            if i < len(onlyfiles) - 1:
                tf.write(s.strip(".csv") + ',')
            else:
                tf.write(s.strip(".csv"))
            i += 1
        tf.write('\n')
        z = 0
        for k in d.items():
            x = 0
            line = ""
            line += str(k[0])
            line += ','
            for v in k[1]:
                line += str(v)
                if x < len(onlyfiles) - 1:
                    line += ','
                x += 1
            if z < len(d.items()) - 1:
                line += '\n'
            z += 1
            tf.write(line)
        total_line = "SUM"
    for i in range(1, filecount + 1):
        total = 0
        cr = csv.reader(open(rank_dict[level] + '.csv'))
        cr.__next__()
        for row in cr:
            try:
                total += int(row[i])
            except ValueError:
                pass
        total_line += ","
        total_line += str(total)
    with open(tf.name, 'a') as sfile:
        sfile.write("\n")
        sfile.write(total_line)
    print (tf.name + " created\n")
    percentage_nanopore(tf.name)


def percentage_nanopore(nanopore):
    """
    Creates a file with the occurrence in percentages of all species or genera.
    :param nanopore: Nanopore input file with rank (species or genera) and count
    :return:
    """
    print ("Reading ", nanopore, "...")
    count = 0
    output = "percentage_" + nanopore
    with open(nanopore) as linecount:
        for line in linecount:
            count += 1
    with open(nanopore) as last:
        for i, l in enumerate(last):
            if i == count-1:
                sum_row = l.split(',')
            else:
                pass
    print ("Writing ", output, "...")
    with open(output, 'w') as perc_nano:
        with open(nanopore) as nano:
            ncount = 0
            header = nano.readline()
            samples = header.split(',')
            for i, dummys in enumerate(samples):
                if i > 0:
                    perc_nano.write(samples[i].strip("\n").strip(".csv") + '_#,' + samples[i].strip("\n").strip(".csv") + '_%,')
                elif i == 0:
                    perc_nano.write(rank_dict[level] + ',')
            perc_nano.write("\n")
            sum_row.pop(0)
            for line in nano:
                if ncount < count:
                    wline = ""
                    sample = line.split(',')
                    wline += sample[0].strip("\n").strip(".csv") + ','
                    for x in range(0, len(samples)-1):
                        percentage = float(100*float(sample[x + 1])/float(sum_row[x]))
                        wline += sample[x + 1].strip("\n") + ',' + str(percentage) + ','
                    perc_nano.write(wline)
                    perc_nano.write("\n")
                    ncount += 1
    print ("Writing ", output, "... finished")


def sort_file(ifilename, ofilename):
    """
    Sorting the nanopore and/or illumina files based on the headers in the fieldnames list
    :param ifilename: Name of the input file with the Nanopore or Illumina info
    :param ofilename: The sorted output file
    :return:
    """
    ifilename = ifilename
    ofilename = ofilename
    with open(ifilename, 'r') as infile, open(ofilename, 'a') as outfile:
        # output dict needs a list for new column ordering. Enter the column names in order.
        fieldnames = [
            '', 193770, 193961, 194125, 185551,  193964, 193963
        ]
        # fieldnames = [
        #     '', 4371, 5013, 4278, 4322, 6043, 4282, 6015, 
        #     5028, 5044, 5006, 4206, 4294, 4364, 5136, 4198, 
        #     5408, 4343, 4195, 5066, 4135, 4243, 4059, 6030, 
        #     5418, 6026, 5164, 4395, 4233, 6119, 4356, 5001, 
        #     5145, 4186, 6003, 4124, 4323, 5082, 4043, 5411, 
        #     6036, 6045, 4187, 4185, 5168, 5027, 4280, 6020, 
        #     5421, 5186, 6107, 6009, 6019, 5067, 5417, 4179, 
        #     5204, 6014, 5103, 6334, 5032,5163, 5304, 5425, 
        #     6035, 4306, 5306, 5015, 6049, 6010, 6074, 
        #     'EC01', 'C1', 'C15'
        # ]
        fields = [rank_dict[level]]
        for fn in fieldnames:
            if fn != "":
                fields.append(str(fn) + '_#')
                fields.append(str(fn) + '_%')
            else:
                fields.append(str(fn))
        writer = csv.DictWriter(outfile, fieldnames=fields, delimiter=',')
        # reorder the header first
        writer.writeheader()
        for row in csv.DictReader(infile):
            # writes the reordered rows to the new file
            writer.writerow(row)


if __name__ == '__main__':
    try:
        if is_taxadb_up_to_date(DEFAULT_TAXADB):
            if rank_dict[level] not in ["species"]:
                print ("Getting " + rank_dict[level] + " names from NCBI...\n")
            else:
                print ("Getting " + rank_dict[level] + " from csv files...\n")
            read_csv()
            # Add the input and output file to sort the csv.
            sort_file("percentage_" + rank_dict[level] + ".csv", rank_dict[level] + "_sorted.csv")
        else:
            print ("Taxonomy database is updating...\n")
            ncbi.update_taxonomy_database()
    except KeyError:
        print ("Please enter one of the options.")
