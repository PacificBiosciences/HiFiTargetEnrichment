import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

DPI=200

readCsv      = snakemake.input.csv
targetBed    = snakemake.input.bed
targetBuffer = int(snakemake.params.buffer)
outDir       = snakemake.params.odir


targets = pd.read_csv(targetBed,
                      sep='\t',
                      names=['chr','start','stop','target'])\
            .set_index('target')
#Set length of target region
targets['tlength'] = targets.eval('stop - start + 2 * @targetBuffer')

data = pd.read_csv(readCsv)
#replace "." with "off-target"
ot = 'off-target'
data.target = data.target.str.replace('.',ot,regex=False)


####
#mean base coverage
####

pdata = data.query(' target != @ot ')\
            .groupby( 'target' )\
            .length.sum().reset_index()

try:
    pdata['meanBaseCoverage'] = pdata.length / pdata.target.map(targets.tlength)
except pd.errors.InvalidIndexError as e:
    print(f'Error: Is your target name list unique?\n\n{e}')
    sys.exit()

order = sorted(pdata.target)
g = sns.catplot(data=pdata,
                x='target',order=order,
                y='meanBaseCoverage',
                kind='bar',aspect=2)
plt.xticks(rotation=45)
g.savefig(f'{outDir}/mean_base_coverage_by_target.png',dpi=DPI)
plt.clf()

###
#sample readlength
###
pdata = data.query( ' target != @ot ')

g = sns.displot(pdata.length,kind='hist',kde=True)
g.savefig(f'{outDir}/readlength_hist.png',dpi=DPI)
plt.clf()

####
#readlength by target
###

pdata = data

order = sorted(pdata.target.unique())
g = sns.FacetGrid(data=pdata,
                  xlim=(0,10000),
                  col='target',col_wrap=7,
                  col_order=order,
                  hue='sample',
                  sharey=False)

g.map(sns.kdeplot,'length')

g.set_xlabels('Read Length')
g.set_ylabels('HiFi Reads')
g.add_legend()
g.savefig(f'{outDir}/readlength_hist_by_target.png',dpi=DPI)
plt.clf()

####
#readlength by target
####
pdata = data
         
order = sorted(pdata.target.unique())
g = sns.catplot( data=pdata,
                 x='target',
                 y='length',
                 order=order,
                 kind='violin',
                 aspect=2)
plt.xticks(rotation=45);
g.savefig(f'{outDir}/readlength_violins_by_target.png',dpi=DPI)
plt.clf()

####
#dedup count by target
####
pdata = data
            
order = sorted(pdata.target.unique())

g = sns.catplot( data=pdata,
                 x='target',
                 y='duplicates',
                 order=order,
                 kind='box',aspect=2)
plt.xticks(rotation=45);
g.savefig(f'{outDir}/dedup_count_by_target.png',dpi=DPI)
plt.clf()


####
#dedup rate by target
####
def dedup_rate(d):
    dups = d.duplicates.sum()
    return  dups / ( dups + len(d) ) 


pdata = data.groupby('target')\
            .apply(dedup_rate)\
            .rename('Duplication Rate')\
            .reset_index()

order = sorted(pdata.target)
g = sns.catplot(data=pdata,
                x='target',order=order,
                y='Duplication Rate',
                kind='bar',aspect=2)
plt.xticks(rotation=45)
g.savefig(f'{outDir}/dedup_rate_by_target.png',dpi=DPI)
plt.clf()

####
#rq by target
###
pdata = data

order = sorted(pdata.target.unique())
            
g = sns.catplot( data=pdata,
                 x='target',
                 y='qual',
                 order=order,
                 kind='box',aspect=2)
plt.xticks(rotation=45)
g.savefig(f'{outDir}/readqual_by_target.png',dpi=DPI)
plt.clf()

###
#off target length
###
pdata = data.query( ' target == @ot ')

g = sns.displot(pdata.length,kind='hist',kde=True)
g.savefig(f'{outDir}/off-target_length_hist.png',dpi=DPI)
