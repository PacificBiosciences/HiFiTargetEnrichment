import pandas as pd

res = pd.concat( [ pd.read_csv( csv, index_col='readname' )
                   for csv in snakemake.input ],
                  axis=1 )\
        .reset_index()
res.insert( 0, 'sample', snakemake.params.sample )
res.to_csv( snakemake.output[0], index=False )
