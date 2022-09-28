import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

DPI=200

coverage = pd.read_csv( snakemake.input[0] )

order = sorted(coverage.target.unique())
maxCols = int( len(order) ** 0.5 ) + 1

plt.figure(figsize=(40,40))

g = sns.FacetGrid(data=coverage,sharex=False,
                  col='target',col_wrap=maxCols,col_order=order,
                  hue='target')

g.map(plt.plot, 'start','coverage')

g.set_xlabels('chr start pos')
g.set_ylabels('Coverage')

g.savefig(f'{snakemake.params.odir}/coverage_by_target.png', dpi=DPI)
