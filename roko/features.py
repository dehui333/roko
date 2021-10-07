import pysam
from Bio import SeqIO
from labels import *
import argparse
from multiprocessing import Pool
import gen
from data import DataWriter

ENCODED_UNKNOWN = encoding[UNKNOWN]

GROUP_SIZE = 10000 # not used?
NUM_WORKERS = 6 # not used?
MAX_INS = 3 # not used?

'''
Returns regions of the reference. A sliding window moves across the reference sequence to demarcate the regions, there
can be overlaps.
regions are named tuples with the format (ref name, start index, end index). The indices are end-exclusive and start at 0.
'''
def generate_regions(ref, ref_name, window=100_000, overlap=300):
    length = len(ref)
    i = 0

    while i < length:
        end = i + window
        yield Region(ref_name, i, min(end, length))

        if end >= length:
            break
        else:
            i = end - overlap

'''
Check if a certain position index is between the first and last bases of any of a collection of alignments. The aligns argument is a list 
of TargetAlign objects.
'''
def is_in_region(pos, aligns):
    for a in aligns:
        if a.start <= pos < a.end:
            return True
    return False

'''
Inputs: 
    bam_X: bam file for alignment of reads to ref(the draft assembly)
    bam_Y: bam file for alignment of truth to ref
    ref: the draft sequence
    region: named tuples with the format (ref name, start index, end index) demarcating a region on the ref
Output:
    A tuple (ref name, list of lists of positions , list of feature matrices, list of  )

'''
def generate_train(args):
    bam_X, bam_Y, ref, region = args # regions are named tuples with the format (ref name, start index, end index)

    # Get mapped and primary alignments overlapping with the region
    alignments = get_aligns(bam_Y, ref_name=region.name, start=region.start, end=region.end)  
    # Filter based on some criteria
    filtered = filter_aligns(alignments)

    print(f'Finished generating labels for {region.name}:{region.start}-{region.end}.')

    if not filtered:
        print('No alignments.')
        return None

    positions, examples, labels = [], [], []

    for a in filtered: # For every alignment of truth to draft
        pos_labels = dict() # Mapping of positions to label
        n_pos = set() # set of positions with unknown

        t_pos, t_labels = get_pos_and_labels(a, ref, region) # What this alignment says about each of the positions in the region
        for p, l in zip(t_pos, t_labels):
            if l == ENCODED_UNKNOWN:
                n_pos.add(p) # No information on the position
            else:
                pos_labels[p] = l

        pos_sorted = sorted(list(pos_labels.keys()))
        
        # Region where the current alignment has information on 
        region_string = f'{region.name}:{pos_sorted[0][0]+1}-{pos_sorted[-1][0]}' # region string, 1 index and both sides inclusive

        # Generate feature matrix for this region
        # result is a tuple (list of lists of positions, list of matrices)
        result = gen.generate_features(bam_X, str(ref), region_string)

        # P is a list of positions corresponding to columns of the matrix, X is the matrix 
        for P, X in zip(*result):
            Y = [] # labels in this window
            to_yield = True

            for p in P:
                assert is_in_region(p[0], filtered) # At least one alignment should have info on this position

                if p in n_pos:
                    to_yield = False # The current alignment covers the position, but has no info on this position itself
                    break

                try:
                    y_label = pos_labels[p]
                except KeyError:
                    if p[1] != 0:
                        y_label = encoding[GAP] # If the truth has no info on this pos, probably false positive for insertion
                    else:
                        raise KeyError(f'No label mapping for position {p}.') # For some reason, though the matrix window is contained
                                                                              # in the region where the alignment has info, there is no label 
                                                                              # for some position?              
                Y.append(y_label)

            if to_yield:
                positions.append(P)
                examples.append(X)
                labels.append(Y)

    print(f'Finished generating examples for {region.name}:{region.start}-{region.end}.')
    return region.name, positions, examples, labels


def generate_infer(args):
    bam_X, ref, region = args

    region_string = f'{region.name}:{region.start+1}-{region.end}'
    result = gen.generate_features(bam_X, ref, region_string)

    positions, examples = [], []

    for P, X in zip(*result):
        positions.append(P)
        examples.append(X)

    print(f'Finished generating examples for {region.name}:{region.start}-{region.end}.')
    return region.name, positions, examples, None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('ref', type=str) # reference path
    parser.add_argument('X', type=str) # path to X bam file
    parser.add_argument('--Y', type=str, default=None) # path to Y bam file
    parser.add_argument('o', type=str) # output path
    parser.add_argument('--t', type=int, default=1) # Number of processes
    args = parser.parse_args()

    # Training if --Y present, else inference.
    inference = False if args.Y else True
    size = 0

    # Parse reference sequences into list of tuples [(id, sequence), (id2, sequence2), ...]
    with open(args.ref, 'r') as handle:
        refs = [(str(r.id), str(r.seq)) for r in SeqIO.parse(handle, 'fasta')]

    with DataWriter(args.o, inference) as data:
        data.write_contigs(refs) # Write references to output path

        func = generate_infer if inference else generate_train

        arguments = []
        for n, r in refs:
            for region in generate_regions(r, n):
                a = (args.X, r, region) if inference else (args.X, args.Y, r, region)
                arguments.append(a)

        print(f'Data generation started, number of jobs: {len(arguments)}.')

        with Pool(processes=args.t) as pool:
            finished = 0
            for result in pool.imap(func, arguments):
                if not result:
                    continue
                c, p, x, y = result
                # name, positions, X, Y
                data.store(c, p, x, y)
                finished += 1

                if finished % 10 == 0:
                    print('Writing to disk started')
                    data.write()
                    print('Writing to disk finished')

            print('Writing to disk started')
            data.write()
            print('Writing to disk finished')


if __name__ == '__main__':
    main()

