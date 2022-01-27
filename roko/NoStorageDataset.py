from Bio import SeqIO
from features import generate_regions, is_in_region 
import gen
from labels import *
import math
import numpy as np
import pysam
import torch

'''
NoStorageDataset takes 2 arguments to initialize: 
1. A list of tuples (path to the bam file for calls to draft, path to the draft fasta, path to the bam file for truth to draft)
2. A boolean, True for inference mode and False otherwise.
If inference mode, the third item in the tuples is not required.
'''

EMPTY_DICT = dict()
ENCODED_UNKNOWN = encoding[UNKNOWN]   

def get_string_tuples_from_file(path):
    outer_list = []
    with open(path) as f:
        for line in f.readlines():
            inner_list = line.split()
            outer_list.append(inner_list)
    return outer_list

def load_arguments_from_file(path, inference):
    tuples = get_string_tuples_from_file(path)
    return generate_args_for_feature_gen(tuples, inference)

# Inputs: A list of tuples of lengths 2 or 3. The items are 1. path to calls to draft bam 2. path to draft fasta 3. path to truth to draft bam (optional).
# Tuples should contain 3. if inference is true.
# Ouput:1. A list of of tuples (path to calls to draft bam, path to truth to draft bam, draft sequence, a region in the draft) if inference is True,
# and path to truth to draft bam is None otherwise. 2. 
def generate_args_for_feature_gen(list_of_string_tuples, inference):
    results = []
    contigs = {}
    for t in list_of_string_tuples:
        bam_X = t[0]
        draft_fasta = t[1]
        bam_Y = None
        if not inference:
            bam_Y = t[2]
        with open(draft_fasta, 'r') as handle:
            drafts = [(str(r.id), str(r.seq)) for r in SeqIO.parse(handle, 'fasta')]
        for name, seq in drafts:
            contigs[name] = (seq, len(seq))
            for region in generate_regions(seq, name):
                a = (bam_X, None, seq, region) if inference else (bam_X, bam_Y, seq, region)
                results.append(a)

    return results, contigs

def worker_init_fn(worker_id):
    worker_info = torch.utils.data.get_worker_info()
    dataset = worker_info.dataset  # the dataset copy in this worker process
    overall_start = dataset.start
    overall_end = dataset.end
    # configure the dataset to only process the split workload
    per_worker = int(math.ceil((overall_end - overall_start) / float(worker_info.num_workers)))
    worker_id = worker_info.id
    dataset.start = overall_start + worker_id * per_worker
    dataset.end = min(dataset.start + per_worker, overall_end)


class NoStorageDataset(torch.utils.data.IterableDataset):
    def __init__(self, input_tuples, infer):
        self.list_of_args, self.contigs = generate_args_for_feature_gen(input_tuples, infer) # unchanged between workers
        self.infer = infer
        # By default, each worker operates on the whole list -> To be modified based on work id 
        self.start = 0
        self.end = len(self.list_of_args)
        
    def __iter__(self):
        for i in range(self.start, self.end):
            bam_X, bam_Y, ref, region = self.list_of_args[i]
            if self.infer:  
                region_string = f'{region.name}:{region.start+1}-{region.end}'
                result = gen.generate_features(bam_X, ref, region_string, EMPTY_DICT, 0)
                for P, X, Y, X2 in zip(*result):
                    P = [list(p) for p in P]
                    yield region.name, torch.tensor(P), X, X2.astype(np.int16)        
                print(f'Finished generating features for {region.name}:{region.start+1}-{region.end}.')
                
            else: 
                alignments = get_aligns(bam_Y, ref_name=region.name, start=region.start, end=region.end)
                filtered = filter_aligns(alignments)
                print(f'Finished generating labels for {region.name}:{region.start+1}-{region.end}.')

                if not filtered:
                    print('No alignments.')
                    return iter(())

                for a in filtered:
                    pos_labels = dict()
                    n_pos = set()

                    t_pos, t_labels = get_pos_and_labels(a, ref, region)
                    for p, l in zip(t_pos, t_labels):
                        if l == ENCODED_UNKNOWN:
                            n_pos.add(p)
                        else:
                            pos_labels[p] = l

                    pos_sorted = sorted(list(pos_labels.keys()))
                    region_string = f'{region.name}:{pos_sorted[0][0]+1}-{pos_sorted[-1][0]}'
                    result = gen.generate_features(bam_X, str(ref), region_string, pos_labels, 1)

                    for P, X, Y, X2 in zip(*result):
                        to_yield = True

                        for p in P:
                            assert is_in_region(p[0], filtered)

                            if p in n_pos:
                                to_yield = False
                                break

                        if to_yield:
                            yield X, Y, X2.astype(np.int16)
                print(f'Finished generating examples for {region.name}:{region.start+1}-{region.end}.')