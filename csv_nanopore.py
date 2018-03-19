import os
import sqlite3
import csv
from os import listdir
from os.path import isfile, join
from collections import Counter
from collections import defaultdict
from ete2 import NCBITaxa


d = defaultdict(list)
taxids = []
mypath = raw_input("Enter path to csv files (i.e. /home/user/csv/files/): ")
allfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
onlyfiles = [s for s in allfiles if '.csv' in s]
filecount = len(onlyfiles)
level = raw_input("1: phylum ----> 2: class ----> 3: order ----> 4: family ----> 5: genus ----> 6: species\nEnter rank number: ")
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
        tax = []
        tids = []
        taxnames = []
        with open(mypath + csv) as nf:
            linenr = 0
            for line in nf:
                if ',,' not in line:
                    line = line.split(',')
                    name = ""
                    if linenr > 0:
                        if float(line[3]) >= 80 and line[8] in ['1', '2']:
                            tids.append(line[4])
                            if rank_dict[level] == "species":
                                # tax.append(line[2])
                                taxid2name = ncbi.get_taxid_translator([int(line[2])])
                                for tid, tname in taxid2name.iteritems():
                                    namelist = str(tname).split(" ")[:2]
                                for nl in namelist:
                                    name += str(nl) + " "
                                taxnames.append(str(name))
                            elif rank_dict[level] == "genus":
                                # tax.append(line[7])
                                taxid2name = ncbi.get_taxid_translator([int(line[7])])
                                for tid, tname in taxid2name.iteritems():
                                    namelist = str(tname).split(" ")[:2]
                                for nl in namelist:
                                    name += str(nl) + " "
                                taxnames.append(str(name))
                            else:
                                taxids = [line[2]]
                                desired_ranks = [rank_dict[level]]
                                for taxid in taxids:
                                    ranks = get_desired_ranks(taxid, desired_ranks)
                                    for key, rank in ranks.items():
                                        if rank != '<not present>':
                                            # tax.append(rank)
                                            taxid2name = ncbi.get_taxid_translator([int(rank)])
                                            for tid, tname in taxid2name.iteritems():
                                                namelist = str(tname).split(" ")[:2]
                                            for nl in namelist:
                                                name += str(nl) + " "
                                            taxnames.append(str(name))
                linenr += 1
            # counter = Counter(tax)
            namecounter = Counter(taxnames)
            for k, v in namecounter.iteritems():
                if k not in d and count > 0:
                    for c in range(0, count):
                        d[k].append(0)
                d[k].append(v)
                klist.append(k)
            for dk, dv in d.iteritems():
                if dk not in klist:
                    d[dk].append(0)
            count += 1
            namecounter.clear()
    for k, v in d.iteritems():
        if len(d[k]) < len(onlyfiles):
            d[k].append(0)
    make_rank_csv()
    print "Finished\n"


def make_rank_csv():
    """
    Create a new csv file based on the rank provided by the user.
    The WIMP csv files will be used as input and the rank will be searched using the NCBI taxid translator.
    :return:
    """
    print "Creating csv file based on rank...\n"
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
        for k in d.iteritems():
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
        cr.next()
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
    print tf.name + " created\n"
    percentage_nanopore(tf.name)


def percentage_nanopore(nanopore):
    """
    Creates a file with the occurrence in percentages of all species or genera.
    :param nanopore: Nanopore input file with rank (species or genera) and count
    :return:
    """
    print "Reading ", nanopore, "..."
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
    print "Writing ", output, "..."
    with open(output, 'w') as perc_nano:
        with open(nanopore) as nano:
            ncount = 0
            header = nano.readline()
            samples = header.split(',')
            for i, s in enumerate(samples):
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
    print "Writing ", output, "... finished"


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
        # output dict needs a list for new column ordering
        fieldnames = [4186, 4323, 4135, 4322, 5013, 4371, 4343, 4059, 4243, 6015, 4195, 6030,
                      4278, 5066, 4395, 6119, 4198, 5136, 5028, 4364, 5006, 6043, 5044, 5408,
                      4282, 4356, 4206, 5145, 5164, 6026, 5418, 4233, 4179, 5186, 4043, 5168,
                      4294, 5015, 5411, 6009, 6019, 5001, 5067, 5417, 6049, 4280, 6036, 4185,
                      5027, 4124, 4187, 6045, 6020, 5082, 6003, 5204, 6334, 6014, 5103, 5306,
                      6107, 5163, 5425, 5304, 5032, 4306, 6035, 5421, 6010, 6074, "C1", "C15",
                      "EC01"]
        fields = ['species']
        for fn in fieldnames:
            fields.append(str(fn) + '_#')
            fields.append(str(fn) + '_%')
        writer = csv.DictWriter(outfile, fieldnames=fields, extrasaction='ignore')
        # reorder the header first
        writer.writeheader()
        for row in csv.DictReader(infile):
            # writes the reordered rows to the new file
            writer.writerow(row)


if __name__ == '__main__':
    try:
        if is_taxadb_up_to_date(DEFAULT_TAXADB):
            if rank_dict[level] not in ["species"]:
                print "Getting " + rank_dict[level] + " names from NCBI...\n"
            else:
                print "Getting " + rank_dict[level] + " from csv files...\n"
            read_csv()
            sort_file("percentage_species.csv", "species_sorted.csv")
        else:
            print "Taxonomy database is updating...\n"
            ncbi.update_taxonomy_database()
    except KeyError:
        print "Please enter one of the options."
