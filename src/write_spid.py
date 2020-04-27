#!/usr/bin/env python3

# read the meraculous config
with open(snakemake.input['spid'], 'rt') as f:
    spid_string = ''.join(f.readlines())

my_conf = spid_string.format(
    snakemake.params['popmap'])

with open(snakemake.output['spid'], 'wt') as f:
    f.write(my_conf)
