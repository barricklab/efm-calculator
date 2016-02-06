####################################################
## EVOLUTIOINARY FAILURE MODE CALCULATOR
## Version 1.0.1
####################################################

EFM_VERSION = "1.0.1"

# Import subprocess for command line stuff
import subprocess
import os
import csv
from tempfile import NamedTemporaryFile

# Import itertools
from itertools import combinations

# Import RegEx libraries
import re
import math

# Import ElementTree for parsing BioBrick XML
from xml.etree import ElementTree

# CLI functionality
import argparse
from Bio import SeqIO, Seq

SUB_RATE = float(2.2 * 10 ** (-10))


def run_mummer(fasta_file, org):
    """
    Executes the repeat-match command in MUMmer and parses the output.
    Executes the nucmer and show-coords commands in MUMmer and parses the output.
    :param fasta_file: Name of FASTA file with sequence
    :param org: Name of organism for calculating rates
    :return: Dictionary of results
    """
    # Initialize dictionary
    long_repeats = dict()
    # Execute MUMmer command
    try:
        output = subprocess.check_output(['repeat-match', '-n', '16', '-f', fasta_file])
    except:
        raise
        print "MUMmer command 'repeat-match' failed."
    output_lines = output.splitlines()
    # Loop through mummer output and skip the first 2 lines which only contain header information
    for line in output_lines[2:]:
        split_line = line.split()
        start_pos1 = int(split_line[0])
        start_pos2 = int(split_line[1])
        repeat_length = int(split_line[2])
        if repeat_length in long_repeats:
            for locations in long_repeats[repeat_length]:
                if start_pos2 in locations and start_pos1 in locations:
                    continue
                elif start_pos1 in locations:
                    locations.append(start_pos2)
                    locations.sort()
                    found = True
                elif start_pos2 in locations:
                    locations.append(start_pos1)
                    locations.sort()
                    found = True
                else:
                    found = False
            if found is False:
                # This is definitely a new repeat of repeat_length that we haven't seen before
                long_repeats[repeat_length].append(sorted([start_pos1, start_pos2]))
        else:
            # We have not seen a repeat of this length before
            long_repeats[repeat_length] = [sorted([start_pos1, start_pos2])]

    # Initialize list of dictionaries and reformat for tabular output
    output_list = []
    seen = []
    for repeat_length, location_list in long_repeats.iteritems():
        for locations in location_list:
            rate_per_repeat = float(0)
            if len(locations) > 2:
                # Must account for a homologous sequence occurring more than twice
                combos = combinations(locations, 2)
                for combo in combos:
                    rate_per_repeat += get_recombo_rate(int([repeat_length][0]), combo[0], combo[1], org)
            else:
                rate_per_repeat = get_recombo_rate(int([repeat_length][0]), locations[0], locations[1], org)
            entry = {'location': locations,
                     'sequence': '', # We don't record the sequence for long repeats
                     'length': [repeat_length],
                     'count': len(locations),
                     'raw_rate': rate_per_repeat,
                     'type': 'rmd'}
            if [locations, [repeat_length]] not in seen:
                output_list.append(entry)
                seen.append([locations, [repeat_length]])

    # Generate a temporary file for storing the '.delta' file that will be passed to show-coords
    coords_file = NamedTemporaryFile(suffix='.delta', delete=False)
    try:
        # Execute nucmer command
        subprocess.call(['nucmer',
                         '-l',
                         '16',
                         '-f',  # Forward strand only
                         '--maxmatch',
                         '--nosimplify',
                         '--prefix='+coords_file.name[:-6],
                         fasta_file,
                         fasta_file],
        )
        # Execute show-coords command
        output = subprocess.check_output(['show-coords', '-T', '-I 100', '-H', '-d', coords_file.name])
    except:
        raise
        print "MUMmer command 'nucmer' and 'show-coords' failed."
    finally:
        # Remove the '.delta' file because we're finished with it
        os.remove(coords_file.name)

    content_lines = output.splitlines()
    # Initialize list of dictionaries for tabular output

    for line in content_lines:
        clean_line = line.split()
        start_pos1 = int(clean_line[0])
        start_pos2 = int(clean_line[2])
        length1 = int(clean_line[4])
        length2 = int(clean_line[5])
        if start_pos1 != start_pos2:
            # Remove repeats by converting to set
            length_list = sorted(list(set([length1, length2])))
            location = sorted([start_pos1, start_pos2])
            rate = get_recombo_rate(int(length_list[0]), location[0], location[1], org)
            entry = {'location': location,
                     'sequence': '',
                     'length': length_list,
                     'count': len(location),
                     'raw_rate': rate,
                     'type': 'rmd'}
            if [location, length_list] not in seen:
                output_list.append(entry)
                seen.append([location, length_list])

    return output_list


def get_repeats_in_window(n, sequence, min_count, org):
    """
    Generates a list of k-mer repeats in a given window size
    :param n: Size of k-mer (2 for dinucleotides, 3 for trinucleotides, etc.)
    :param sequence: String of DNA sequence
    :param min_count: Minimum number of repeating units before recording
    :param org: Host organism for rate calculations
    :return: List of dictionaries of results for tabular output
    """
    # Initialize repeats dictionary, in the format repeats -> { 'index' : ['CG', 2] }
    repeats = dict()
    if n == 0:
        # Return an empty list, n = 0 doesn't make sense
        return []

    # Start counting repeats at 1
    repeat_count = 1
    i = 0
    if n == 1:
        # n = 1 is a special, simpler case to handle
        while i < (len(sequence) - n):
            if sequence[i] == sequence[i + 1]:
                # Record as a repeat
                if repeat_count + 1 >= min_count:
                    repeats[i + n - n * repeat_count] = [sequence[i], repeat_count + 1]
                # Jump one base ahead
                i += 1
                repeat_count += 1
            else:
                # Reset repeat_count
                i += 1
                repeat_count = 1
    else:
        # n > 1
        while i < (len(sequence) - n):
            # If the first base of the current window matches the first base of the next window
            # AND the entire window equals the next consecutive window
            if sequence[i] == sequence[i + n] and sequence[i:i + n] == sequence[(i + n):(i + 2 * n)]:
                # If the window is larger than 4, make sure it doesn't contain any smaller repeating subunits
                if check_subunits(sequence[i:i + n]) is True:
                    i += 1
                else:
                    # Record as a repeat
                    if repeat_count + 1 >= min_count:
                        repeats[i + n - n * repeat_count] = [sequence[i:i + n], repeat_count + 1]
                    # Jump forward a full window of bases
                    i = i + n
                    repeat_count += 1
            else:
                repeat_count = 1
                # Otherwise shift the window forward by one base
                i += 1

    # Reformat for tabular output
    output_list = []
    for index, contents in repeats.iteritems():
        mut_rate = get_mut_rate(contents[1], n, org)
        entry = {'location': [index],
                 'sequence': str(contents[0]),
                 'length': [len(contents[0])],
                 'count': contents[1],
                 'raw_rate': mut_rate,
                 'type': 'ssr'
        }
        output_list.append(entry)

    return output_list


def check_subunits(sequence):
    """
    Checks to see if a sequence can be broken up into smaller repeating subunits. For example, ATAT is actually just
    two instances of AT.
    :param sequence: Short DNA sequence
    :return: Boolean. True if sequence can be broken up into smaller repeating subunits. Otherwise, False.
    """
    if re.match(r"^(.+?)\1+$", sequence):
        return True
    else:
        return False


def get_mut_rate(repeat_count, unit_length, org):
    '''
    Calculates mutation rate for simple sequence repeats
    :param repeat_count: Number of times the repeating unit occurs
    :param unit_length: Length of repeating unit
    :param org: Host organism
    :return: Mutation rate
    '''
    mut_rate = float(0)
    if org == 'ecoli' or org == 'reca':
        if unit_length == 1:
            # Formula based on analysis of Lee et. al. data
            mut_rate = float(10 ** (0.72896 * repeat_count - 12.91471))
        elif unit_length > 1:
            mut_rate = float(10 ** (0.06282 * repeat_count - 4.74882))
    elif org == 'yeast':
        if unit_length == 1:
            mut_rate = float(10 ** (0.3092 * repeat_count - 7.3220))
        elif unit_length > 1:
            mut_rate = float(10 ** (0.11141 * repeat_count - 7.65810))
    return mut_rate


def get_recombo_rate(length, location1, location2, org):
    '''
    Calculate the recombination rate based on the Oliviera, et. al. formula
    :param length: Length of homologous region
    :param location1: Location of first homologous region
    :param location2: Location of second homologous region
    :param org: Host organism
    :return: Recombination rate
    '''
    spacer = abs(int(location2) - int(location1)) - int(length)
    # If the homologous sequences overlap we can't calculate a rate
    if spacer < 0:
        return 0
    if org == 'ecoli' or org == 'yeast':
        recombo_rate = float(((8.8 + spacer) ** (-29.0 / length)) * (length / (1 + 1465.6 * length)))
    elif org == 'reca':
        recombo_rate = float(
            ((200.4 + spacer) ** (-8.8 / length)) * (length / (1 + 2163.0 * length + 14438.6 * spacer)))

    return recombo_rate


def get_biobrick_features(tree):
    """
    Extracts sequence annotations from BioBrick XML files
    :param tree: parsed XML tree
    :return: List of dictionaries of features to display
    """
    features = []
    colors = {
        'promoter': 'green',
        'stop': 'red',
        'cds': 'blue',
        'rbs': 'orange',
        'binding': 'purple',
        'BioBrick': 'black'
    }
    for node in tree.iter(tag='feature'):
        if node.find('type').text in colors:
            entry = dict()
            entry['type'] = node.find('type').text
            entry['title'] = node.find('title').text
            entry['startpos'] = int(node.find('startpos').text)
            entry['length'] = int(node.find('endpos').text) - int(node.find('startpos').text)
            entry['color'] = colors[entry['type']]
            features.append(entry)

    return features


def get_genbank_features(genome):
    """
    Extracts sequence annotations a BioPython SeqRecord object
    :param genome: BioPython SeqRecord of genbank file
    :return: List of dictionaries of features to display
    """
    features = []
    for feature in genome.features:
        if feature.qualifiers:
            name = feature.qualifiers.itervalues().next()[0]
        else:
            name = 'Untitled'
        entry = dict()
        entry['type'] = str(feature.type)
        entry['title'] = name
        entry['startpos'] = int(feature.location.nofuzzy_start)
        entry['length'] = int(feature.location.nofuzzy_end) - int(feature.location.nofuzzy_start)
        features.append(entry)

    return features


def truncate_table(repeats, min_rate):
    """
    Returns a subset of repeats with a mutation rate higher than min_rate
    :param repeats: List of dictionaries of repeats
    :param min_rate: Minimum mutation rate to display in table
    :return: A sorted and truncated list of dictionaries of repeats
    """
    trunc = (repeat for repeat in repeats if
             repeat['raw_rate'] > min_rate and repeat['raw_rate'] != '' and repeat['overlap'] == True)
    trunc_sort = sorted(trunc, key=lambda k: k['raw_rate'] if k['raw_rate'] != '' else 0, reverse=True)
    return trunc_sort


def check_overlap(repeats, features, check_features):
    """
    Checks each repeat to see if it overlaps with a feature, and sets 'overlap' to True if it does.
    :param repeats: List of dictionaries of repeats
    :param features: List of dictionaries of features
    :return: List of dictionaries of repeats, with each repeat containing a boolean 'overlap' value
    """
    # If there are no features to report set all overlap to true to display everything
    if check_features is False:
        for repeat in repeats:
            repeat['overlap'] = True
        return repeats
    else:
        for repeat in repeats:
            for feature in features:
                if repeat['type'] == 'rmd':
                    repeat_range = set(xrange(repeat['location'][0], repeat['location'][-1] + repeat['length'][0]))
                else:
                    repeat_range = set(
                        xrange(repeat['location'][0], repeat['location'][0] + repeat['length'][0] * repeat['count']))
                feature_range = set(xrange(feature['startpos'], feature['startpos'] + feature['length']))
                # Magical ampersand checks if two sets overlap
                if repeat_range & feature_range:
                    repeat['overlap'] = True
                    break
                else:
                    repeat['overlap'] = False

    return repeats


def rate_sum(repeats, seq_len):
    """
    Calculates an RIP score for given sequence
    :param repeats: List of dictionaries of repeats
    :param seq_len: Length of input sequence
    :return: Total predicted RIP score for whole sequence
    """
    ssr_sum = float(0)
    rmd_sum = float(0)
    for entry in repeats:
        if entry['raw_rate'] != '' and entry['overlap'] is True:
            if entry['type'] == 'ssr':
                ssr_sum += entry['raw_rate']
            elif entry['type'] == 'rmd':
                rmd_sum += entry['raw_rate']

    base_rate = float(seq_len) * float(SUB_RATE)
    # Add in the mutation rate of an individual nucleotide
    r_sum = ssr_sum + rmd_sum + base_rate
    # Set the maximum rate sum to 1 for now.
    if r_sum > 1:
        r_sum = float(1)
    rel_rate = (float(r_sum) / float(base_rate))

    return {'rip': rel_rate, 'ssr_sum': ssr_sum, 'rmd_sum': rmd_sum, 'bps_sum': base_rate}


def process_efm(fasta_filepath, features, my_seq, org, check_features):
    """
    Takes command line arguments and finds potentially hypermutable sites in a submitted sequence(s).
    :param fasta_filepath: the path of the FASTA file
    :param features: dictionary of features
    :param org: the organism to which the sequence belongs
    :param check_features: boolean determining to check for features
    :return: A CSV format to be written to disk
    """
    # Define the paths for your input file
    input_file = fasta_filepath
    

    # Set maximum window size for get_repeats_in_window()
    unit_length = 15

    # Integrate repeats generated from nucmer and from repeat-match
    mummer_repeats = run_mummer(input_file, org)

    # Get short repeats (SSRs) in each window size up to max unit_length
    all_ssr = []
    for i in range(unit_length):
        min_count = math.ceil(8 / (i + 1))
        if min_count < 3:
            min_count = 3
        if (i + 1) == 1:
            min_count = 4
        repeat_in_window = get_repeats_in_window(i + 1, my_seq, min_count, org)
        if repeat_in_window:
            all_ssr += repeat_in_window

    # Merge repeat lists together
    merged_repeats = mummer_repeats + all_ssr

    # Check if any areas overlap annotated regions
    merged_repeats = check_overlap(merged_repeats, features, check_features)

    # Truncate repeat list based on absolute mutation rate
    merged_repeats_trunc = truncate_table(merged_repeats, 10 ** (-9))

    # Find the sum of all mutation rates for sequences.
    overall_rate = rate_sum(merged_repeats, len(my_seq))

    return {'repeats': merged_repeats_trunc if merged_repeats_trunc else '',
            'features': features,
            'seq_length': len(my_seq),
            'rate': overall_rate,
            'title': 'title',
            'check_features': check_features,
            'organism': org,
            'version': EFM_VERSION}
    


def process_file(filepath, organism):
    """Process a single file given by the user on the command line."""

    fasta_filepath = filepath

    # Determine file type
    spstring = re.split('/', filepath)
    fname = spstring[-1].lower()
    fnamesplit = re.split('\.', fname)
    ftype = fnamesplit[-1]
    if ftype == 'gb':
        ftype = 'genbank'
    check_features = (ftype == 'genbank')

    # Open the file and get metadata
    obj_file = SeqIO.read(filepath, ftype)
    features = get_genbank_features(obj_file)
    my_seq = str(obj_file.seq)
    
    # Create FASTA file if necessary
    if ftype != 'fasta':
        fasta_filepath = "/tmp/" + fnamesplit[0] + ".fasta"
        with open(fasta_filepath, 'w') as handle:
            SeqIO.convert(filepath, "genbank", handle, "fasta")

    # Process the file
    output_dict = process_efm(fasta_filepath, features, my_seq, organism, check_features)
    return output_dict

def process_dir(dirpath, organism):
    """Process several files all contained within a directory specified
    on the command line."""

    flist = list() # list of files
    metadata_lst = list() # list of dictionaries containing sequence metadata

    for dirName, subDirList, fileList in os.walk(dirpath):
        flist = fileList

    # Gather metadata on all files and store in a list
    for f in flist:
        metadata_lst.append(process_file(f, organism))
    
    # todo: finish this, lol
    for item in metadata_lst:
        print "Item in lst:\n" + str(item)

    return

def main():
    """Driver for command line version of EFM calculator."""

    """Get the input from the user via the command line. The user may
    specify (currently) either GenBank or FASTA files. GenBank files
    must be converted into FASTA format for use with repeat-match in 
    MUMmer."""

    # Define temporary directory
    tmp_path = "/tmp/"

    # Define error strings
    ERR_NO_FILE = "No file(s) specified."
    ERR_SPEC_ORG = "Must specify organism if using FASTA file."
    ERR_NO_ACCESS = "Cannot access path:"
    parser = argparse.ArgumentParser()
    parser.add_argument('file', type=str, help="file or path containing sequence(s) to analyze")
    #parser.add_argument('--format', type=str, help="specify format of sequence file(s)")
    parser.add_argument('--organism', choices=['ecoli', 'reca', 'yeast'], help="specify organism")
    args = parser.parse_args()

    # Enforce constraints on input
    organism = ""

    if not args.file:
        print "ERROR:", ERR_NO_FILE
        return

    if not(os.access(args.file, os.R_OK)):
        print "ERROR:", ERR_NO_ACCESS, args.file
        return

    if not args.organism:
        organism = 'ecoli'
    else:
        organism = args.organism

    # Process a single file
    if not (os.path.isdir(args.file)):
        d = process_file(args.file, organism)
        print "OUTPUT:\n" + str(d)

    else:
        process_dir(args.file, organism)

if __name__ == "__main__":
    main()
