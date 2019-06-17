import os
import sqlite3
import csv
import re
import random
from os import listdir
from os.path import isfile, join
from collections import Counter
from collections import defaultdict
from collections import OrderedDict
from ete3 import NCBITaxa
from operator import itemgetter
from datetime import datetime, date, time


d = defaultdict(list)
ncbi = NCBITaxa()
DEFAULT_TAXADB = os.path.join(os.environ.get(
    'HOME', '/'), '.etetoolkit', 'taxa.sqlite')
DB_VERSION = 2
barcodes = ["BC01", "BC02", "BC03", "BC04", "BC05", "BC06",
            "BC07", "BC08", "BC09", "BC10", "BC11", "BC12"]


def get_input():
    """Get all user input and return all files and settings.
    
    Returns:
        Filepaths and all QC and classification files.
        Searchranks that will be added to the OTU table.
        Minimum qscore used for filtering the reads.
    """
    while True:
        mypath = input("Enter classification files path: ")
        mypathqc = input("Enter QC files path: ")
        minqscore = input("Please enter the minimun qc score per read: ")
        if mypath == "" or mypathqc == "" or minqscore == "":
            print("Invalid input.")
            get_input()
        else:
            allfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
            qcfiles = [f for f in listdir(mypathqc) if isfile(join(mypathqc, f))]
            onlyfiles = [s for s in qcfiles if '22.csv' in s]
            searchrank = ["family", "genus"]
            return mypath, mypathqc, minqscore, allfiles, qcfiles, onlyfiles, searchrank


def is_taxadb_up_to_date(dbfile=DEFAULT_TAXADB):
    """Check if a valid and up-to-date taxa.sqlite database exists.
    If dbfile is not specified, DEFAULT_TAXADB is assumed
    
    Keyword Arguments:
        dbfile: The NCBI taxonomy database file to check. (default: {DEFAULT_TAXADB})
    
    Returns:
         True if the database is up-to-date or false if it is not.

    Raises:
        sqlite3.OperationalError, ValueError, IndexError, TypeError:
        Set the NCBI taxonomy database version to None if there is no version found.

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
    """Finds the rank entered by the user for specific taxonomy IDs.
    
    Arguments:
        taxid: Taxonomy ID from EPI2ME csv file.
        desired_ranks: The desired rank to find based on the taxonomy ID.

    Returns:
        Rank based on taxid
    """
    lineage = ncbi.get_lineage(taxid)
    names = ncbi.get_taxid_translator(lineage)
    lineage2ranks = ncbi.get_rank(names)
    ranks2lineage = dict((rank, taxid)
                         for (taxid, rank) in lineage2ranks.items())
    return{'{}_id'.format(rank): ranks2lineage.get(rank, '<not present>') for rank in desired_ranks}


def read_basecalling_qc():
    """Read the QC file from EPI2ME to filter specific reads 
    from the classification file.

    Returns:
        A list of reads that passed the filtering steps.
        A dictionary with read counts per barcode.
    """
    print()
    print("Filtering reads...")
    print()
    ok_read_ids = []
    barcode_list = []
    for qcfile in onlyfiles:
        read_id_column = 1
        barcode_column = 2
        seqlen_column = 6
        mean_qscore_column = 7
        with open(mypathqc + qcfile) as qf:
            linenr = 0
            for line in qf:
                if ',,' not in line:
                    line = line.split(',')
                    if line[barcode_column] != "NA":
                        if linenr > 0:
                            if int(line[seqlen_column]) >= 1400 and int(line[seqlen_column]) <= 1700:
                                if float(line[mean_qscore_column].strip('\n')) >= int(minqscore):
                                    ok_read_ids.append(line[read_id_column])
                                    barcode_list.append(line[barcode_column])
                linenr += 1
        barcode_dict = Counter(barcode_list)
        print(len(ok_read_ids), "reads will be used")
        print()
    return ok_read_ids, barcode_dict


def read_csv():
    """Read the csv files and build lists with available taxids and occurance.

    Raises:
        IndexError: Set rank to (no rank) is none is found.
    """
    taxnames = []
    allnames, headers, randombc = get_all_species(ok_read_ids)
    for csv in onlyfiles:
        read_id_column = 1
        barcode_column = 5
        acc_column = 6
        taxid_column = 4
        for bc in barcodes:
            print("Checking", bc, "in", csv)
            with open(mypath + csv) as nf:
                linenr = 0
                for line in nf:
                    if ',,' not in line:
                        line = line.split(',')
                        name = ""
                        if line[barcode_column] == bc:
                            if linenr > 0:
                                if (
                                    float(line[acc_column]) >= 80 and
                                    line[read_id_column] in ok_read_ids
                                ):
                                    taxid2name = ncbi.get_taxid_translator([int(line[taxid_column])])
                                    bestrankdict = ncbi.get_rank([int(line[taxid_column])])
                                    bestrank = list(bestrankdict.values())
                                    for dummytid, tname in taxid2name.items():
                                        namesplit = tname.split(' ')
                                        if len(namesplit) > 2:
                                            splitnr = 0
                                            for split in namesplit[:2]:
                                                if splitnr < 1:
                                                    name += split + ' '
                                                else:
                                                    name += split
                                                splitnr += 1
                                        else:
                                            name = str(tname)
                                    if line[barcode_column] == bc:
                                        try:
                                            fullname = str(name + " (" + bestrank[0] + ")")
                                        except IndexError:
                                            fullname = str(name + " (no rank)")
                                        taxnames.append(str(fullname))
                    linenr += 1
                namecounter = Counter(taxnames)
                for k, v in namecounter.items():
                    values = []
                    values.append(v)
                    values.append(csv.strip(".csv") + bc)
                    d[k].append(values)
                allnames = list(dict.fromkeys(allnames))
                for key in allnames:
                    if key not in namecounter.keys():
                        values = []
                        values.append(0)
                        values.append(csv.strip(".csv") + bc)
                        d[key].append(values)
                namecounter.clear()
                taxnames = []
                print(bc, "done!")
    make_rank_csv(headers)
    print("Finished\n")


def get_all_species(ok_read_ids):
    """Get all species that are found in the csv files.
    
    Arguments:
        ok_read_ids: A list with all read ids that passed the filtering stage.
    
    Returns:
        All species names and header needed to create a new species file.
        Random line number per barcode for future subset functionality.
    
    Raises:
        IndexError: Set rank to (no rank) if no rank is found in the NCBI taxonomy database.
    """
    randombc = {}
    for barcode, bcount in barcode_dict.items():
        values = []
        if bcount > 400:
            for dummyx in range(420):
                values.append(random.randint(1,bcount))
        else:
            for x in range(1, bcount):
                values.append(x)
        randombc[barcode] = values
    allnameslist = []
    headers = []
    for csv in onlyfiles:
        print()
        print("Getting all available species from " + csv + "...")
        print()
        read_id_column = 1
        barcode_column = 5
        acc_column = 6
        taxid_column = 4
        for bc in barcodes:
            headers.append(csv.strip(".csv") + bc)
        with open(mypath + csv) as nf:
            print("Reading csv file...")
            print()
            linenr = 0
            for line in nf:
                if ',,' not in line:
                    line = line.split(',')
                    name = ""
                    if linenr > 0:
                        if (
                            float(line[acc_column]) >= 80 and
                            line[read_id_column] in ok_read_ids and
                            line[barcode_column] != "NA"
                        ):
                            bestrankdict = ncbi.get_rank([int(line[taxid_column])])
                            bestrank = list(bestrankdict.values())
                            taxid2name = ncbi.get_taxid_translator([int(line[taxid_column])])
                            for dummytid, tname in taxid2name.items():
                                namesplit = tname.split(' ')
                                if len(namesplit) > 2:
                                    splitnr = 0
                                    for split in namesplit[:2]:
                                        if splitnr < 1:
                                            name += split + ' '
                                        else:
                                            name += split
                                        splitnr += 1
                                else:
                                    name = str(tname)
                            if str(name) not in allnameslist:
                                try:
                                    fullname = str(name + " (" + bestrank[0] + ")")
                                except IndexError:
                                    fullname = str(name + " (no rank)")
                                allnameslist.append(str(fullname))
                linenr += 1
    allnames = list(set(allnameslist))
    print(len(allnames), "species found!")
    print()
    return allnames, headers, randombc


def make_rank_csv(headers):
    """Generate an OTU table with species occurance.

    Arguments:
        headers: List of header names with the run and barcode names.

    Raises:
        ValueError, IndexError and StopIteration: Catch errors 
        when calculating the sum of the barcodes.
    """
    print()
    print("Creating csv file based on rank...\n")
    print()
    with open('otu_count.csv', 'w') as tf:
        family = ""
        genus = ""
        tf.write('family,genus,best (rank),')
        i = 0
        for h in headers:
            tf.write(h)
            if i < len(headers) - 1:
                tf.write(',')
            i += 1
        tf.write('\n')
        for key, value in sorted(d.items(), key=itemgetter(1), reverse=True):
            x = 0
            line = ""
            keyname = re.sub(r"[\(\[].*?[\)\]]", "", key)
            name2taxid = ncbi.get_name_translator([keyname[:-1]])
            for dummytaxonomyname, taxonomyid in name2taxid.items():
                taxids = [taxonomyid[0]]
                # desired_ranks = searchrank
                for taxid in taxids:
                    for srank in searchrank:
                        ranks = get_desired_ranks(taxid, [srank])
                        for dummykey, rank in ranks.items():
                            if rank != '<not present>':
                                taxid2name = ncbi.get_taxid_translator([int(rank)])
                                for dummytid, tname in taxid2name.items():
                                    if srank == "family":
                                        family = str(tname)
                                    elif srank == "genus":
                                        genus = str(tname)
                            else:
                                if srank == "family":
                                    family = "NA"
                                elif srank == "genus":
                                    genus = "NA"
            line += str(family + ',' + genus + ',')
            line += str(key + ',')
            for v in value:
                line += str(v[0])
                if x < len(headers) - 1:
                    line += ','
                else:
                    line += '\n'
                x += 1
            tf.write(line)
        total_line = "SUM,,"
        for i in range(3, len(headers) + 3):
            try:
                total = 0
                cr = csv.reader(open('otu_count.csv'))
                cr.__next__()
                for row in cr:
                    try:
                        total += int(row[i])
                    except (ValueError, IndexError):
                        pass
            except StopIteration:
                pass
            total_line += ","
            total_line += str(total)
    with open(tf.name, 'a') as sfile:
        sfile.write(total_line)
    print(tf.name + " created\n")
    print()
    percentage_nanopore(tf.name)


def percentage_nanopore(filename):
    """Calculate percentage of species occurance and generate a 
    new OTU table with percentages.
    
    Arguments:
        filename: Filename of the OTU table.

    Raises:
        ZeroDivisionError: Set percantage to 0 if count is 0.
    """
    print("Create percentage OTU table...")
    count = 0
    output = "otu_percentage.csv"
    with open(filename) as linecount:
        for line in linecount:
            count += 1
    with open(filename) as last:
        for i, l in enumerate(last):
            if i == count-1:
                sum_row = l.split(',')
            else:
                pass
    print("Writing ", output, "...")
    with open(output, 'w') as perc_nano:
        with open(filename) as nano:
            ncount = 0
            header = nano.readline()
            samples = header.split(',')
            for i, dummys in enumerate(samples):
                headerline = ""
                if i == 0:
                    headerline += 'family,'
                elif i == 1:
                    headerline += 'genus,'
                elif i == 2:
                    headerline += 'best (rank),'
                elif i > 2:
                    headerline += samples[i].strip("\n").strip(".csv")
                if i > 2 and i < len(samples)-1:
                    headerline += ","
                perc_nano.write(headerline)
            perc_nano.write("\n")
            sum_row.pop(0)
            for line in nano:
                if ncount < count:
                    wline = ""
                    sample = line.split(',')
                    wline += sample[0].strip("\n").strip(".csv") + ','
                    wline += sample[1].strip("\n").strip(".csv") + ','
                    wline += sample[2].strip("\n").strip(".csv") + ','
                    for x in range(2, len(sample) - 1):
                        try:
                            percentage = float(
                                100*float(sample[x + 1].strip('\n'))/float(sum_row[x]))
                        except ZeroDivisionError:
                            percentage = 0
                        wline += str(percentage)
                        if x < len(sample)-2:
                            wline += ','
                        elif x > len(sample)-2:
                            wline += '\n'
                    perc_nano.write(wline)
                    perc_nano.write("\n")
                    ncount += 1
    print("Percentage OTU table finished!")
    print()


if __name__ == '__main__':
    """Check if database is up-to-date.
    If database is up-to-date start the OTU table creation.
    If database is not up-to-date the database will first be updated.
    """
    if is_taxadb_up_to_date(DEFAULT_TAXADB):
        start = datetime.now()
        print("--------------------------------------------------------------------")
        print("Starting OTU script at", datetime.now().strftime("%d-%m-%Y %H:%M:%S"))
        print("--------------------------------------------------------------------")
        mypath, mypathqc, minqscore, allfiles, qcfiles, onlyfiles, searchrank = get_input()
        ok_read_ids, barcode_dict = read_basecalling_qc()
        read_csv()
        end = datetime.now()
        runtime = end - start
        totalruntime = str(runtime).split(".")[0]
        hours = int(totalruntime.split(":")[0])
        minutes = int(totalruntime.split(":")[1])
        seconds = int(totalruntime.split(":")[2])
        days = int(hours / 24)
        print("--------------------------------------------------------------------")
        print("Finished OTU script at", datetime.now().strftime("%d-%m-%Y %H:%M:%S"))
        if days > 0:
            print("Total runtime:", str(days), "d", str(hours), "h", str(minutes), "m", str(seconds), "s")
        elif hours > 0:
            print("Total runtime:", str(hours), "h", str(minutes), "m", str(seconds), "s")
        elif minutes > 0:
            print("Total runtime:", str(minutes), "m", str(seconds), "s")
        else:
            print("Total runtime:", str(seconds), "seconds")
        print("--------------------------------------------------------------------")
    else:
        print("--------------------------------------------------------------------")
        print("Taxonomy database is updating...")
        print("--------------------------------------------------------------------")
        ncbi.update_taxonomy_database()
        print("--------------------------------------------------------------------")
        print("Update finished!!!")
        print("--------------------------------------------------------------------")