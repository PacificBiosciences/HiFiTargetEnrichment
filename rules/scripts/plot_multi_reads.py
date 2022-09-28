import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

DPI=400

readCsv       = snakemake.input.csv
targetBed     = snakemake.input.bed
targetBuffer  = int(snakemake.params.buffer)
targetPerPlot = int(snakemake.params.targetsPerPanel)
outDir        = snakemake.params.odir

# function for partitioning target sets into readable subplots
def partition(lst, size):
    for i in range(0, len(lst), size):
        yield lst[i:i+size]

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

# assign groups to partition targets
order = sorted(pdata.target.unique())
grps = { tgt:i for i,grp in enumerate( partition( order, targetPerPlot ) ) for tgt in grp }
pdata[ 'plotGroup' ] = pdata.target.map(grps)

g = sns.catplot(data=pdata,
                x='target',sharex=False,
                y='meanBaseCoverage',
                col='plotGroup', col_wrap=2,
                kind='box',aspect=2)
g.set_xticklabels(rotation=45);
# remove plotgroup labels
for ax in g.axes.flatten():
    ax.set_title('')
plt.tight_layout()
g.savefig(f'{outDir}/mean_base_coverage.png',dpi=DPI)
plt.clf()

g = sns.catplot(data=pdata,
                x='target', sharex=False,
                y='meanBaseCoverage',
                col='plotGroup', col_wrap=2,
                hue='sample',
                kind='strip',aspect=2 )
                #facet_kws=dict(legend_out=True))
g.set_xticklabels(rotation=45);
for ax in g.axes.flatten():
    ax.set_title('')
sns.move_legend(g, "upper left", bbox_to_anchor=(1, 0.75), frameon=False)
plt.tight_layout()
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
grps = { tgt:i for i,grp in enumerate( partition( order, targetPerPlot ) ) for tgt in grp }
pdata[ 'plotGroup' ] = pdata.target.map(grps)

g = sns.catplot(data=pdata,
                x='target', sharex=False,
                col='plotGroup',col_wrap=2,
                y='Duplication Rate',
                kind='box',aspect=2)
g.set_xticklabels(rotation=45);
for ax in g.axes.flatten():
    ax.set_title('')
plt.tight_layout()
g.savefig(f'{outDir}/dedup_rate_by_target.png',dpi=DPI)
plt.clf()

####
#readlength by target
####
#TODO better represent dups by replicating length dup times
print('writing dedup_length')
pdata = data
         
order = sorted(pdata.target.unique())
grps = { tgt:i for i,grp in enumerate( partition( order, targetPerPlot ) ) for tgt in grp }
pdata[ 'plotGroup' ] = pdata.target.map(grps)
g = sns.catplot( data=pdata,
                 x='target', sharex=False,
                 col='plotGroup',col_wrap=2,
                 y='length',
                 kind='violin',aspect=2)
g.set_xticklabels(rotation=45);
for ax in g.axes.flatten():
    ax.set_title('')
plt.tight_layout()
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
