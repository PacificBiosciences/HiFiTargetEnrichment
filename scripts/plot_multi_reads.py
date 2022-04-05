import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

DPI=400

readCsv      = sys.argv[1]
targetBed    = sys.argv[2]
targetBuffer = int(sys.argv[3])
outDir       = sys.argv[4]

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
print('writing mean base cov')
pdata = data.query(' target != "off-target" ')\
            .groupby( ['target','sample'] )\
            .length.sum().reset_index()

pdata['meanBaseCoverage'] = pdata.length / pdata.target.map(targets.tlength)


order = sorted(pdata.target.unique())
g = sns.catplot(data=pdata,
                x='target',order=order,
                y='meanBaseCoverage',
                kind='box',aspect=2)
plt.xticks(rotation=45)
g.savefig(f'{outDir}/mean_base_coverage.png',dpi=DPI)
plt.clf()

g = sns.catplot(data=pdata,
                x='target',order=order,
                y='meanBaseCoverage',
                hue='sample',
                kind='strip',aspect=2)
plt.xticks(rotation=45)
g.savefig(f'{outDir}/mean_base_coverage_by_sample.png',dpi=DPI)
plt.clf()


####
#dedup rate by target
####
print('writing dedup_rate')
def dedup_rate(d):
    dups = d.duplicates.sum()
    return  dups / ( dups + len(d) ) 


pdata = data.groupby(['target','sample'])\
            .apply(dedup_rate)\
            .rename('Duplication Rate')\
            .reset_index()

order = sorted(pdata.target.unique())
g = sns.catplot(data=pdata,
                x='target',order=order,
                y='Duplication Rate',
                kind='box',aspect=2)
plt.xticks(rotation=45);
g.savefig(f'{outDir}/dedup_rate_by_target.png',dpi=DPI)
plt.clf()

####
#readlength by target
####
#TODO better represent dups by replicating length dup times
print('writing dedup_length')
pdata = data
         
order = sorted(pdata.target.unique())
g = sns.catplot( data=pdata,
                 x='target',
                 y='length',
                 kind='violin',aspect=2)
plt.xticks(rotation=45);
g.savefig(f'{outDir}/dedup_length_by_target.png',dpi=DPI)
plt.clf()

####
#readlength by target
####
print('writing readlength by target')
pdata = data
wrap  = int( len(order) ** 0.5 + 1 )         
g = sns.FacetGrid(data=pdata,
                  hue='sample',
                  col='target',col_wrap=wrap,
                  col_order=order,sharey=False)
order = sorted(pdata.target.unique())



g.map(sns.kdeplot,'length')

g.set_xlabels('Read Length')
g.set_ylabels('HiFi Reads')
g.add_legend()
g.savefig(f'{outDir}/readlength_hist_by_target.png',dpi=DPI)
plt.clf()
