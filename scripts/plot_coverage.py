#! /home/UNIXHOME/jharting/anaconda3/bin/python3.7

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

DPI=400

covCsv = sys.argv[1]
odir   = sys.argv[2]

coverage = pd.read_csv(covCsv)

order = sorted(coverage.target.unique())

g = sns.FacetGrid(data=coverage,sharex=False,
                  col='target',col_wrap=7,col_order=order,
                  hue='target')

g.map(plt.plot, 'start','coverage')

g.set_xlabels('chr start pos')
g.set_ylabels('Coverage')

g.savefig(f'{odir}/coverage_by_target.png', dpi=DPI)
