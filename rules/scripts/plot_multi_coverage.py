import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from math import ceil

DPI=400

data = pd.read_csv( snakemake.input[0] )
#replace "." with "off-target"
ot = 'off-target'
data.target = data.target.str.replace('.',ot,regex=False)

####
#multi-sample coverage per target
####

order = sorted(data.target.unique())
ncols = ceil(len(order)**.5)

g = sns.FacetGrid(data=data,sharex=False,
                  col='target',col_wrap=ncols,col_order=order,
                  hue='sample')

g.map(plt.plot, 'start','coverage')\
 .set_xlabels('chr start pos')\
 .set_ylabels('Coverage')\
 .add_legend()\
 .savefig(f'{snakemake.params.odir}/multi_coverage_by_target.png',dpi=DPI)
