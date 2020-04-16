#!/usr/bin/env python3

import pandas
import pickle

key_file = snakemake.input['key_file']
pickle_file = snakemake.output['pickle']

# read key file
key_data = pandas.read_csv(key_file, delimiter='\t')

# remove spaces from mararoa-downs
key_data['sample'] = key_data['sample'].str.replace('\s', '-', regex=True)
key_data['sample'] = key_data['sample'].str.replace('.', '-', regex=False)

all_individuals = sorted(set(key_data['sample']))

my_individuals = enumerate(all_individuals)
individual_i = {y: x for x, y in my_individuals}
# pickle the individual_i dict for other rules to use
with open(pickle_file, 'wb+') as f:
    pickle.dump(individual_i, f)
