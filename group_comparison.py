#!/usr/bin/env python
#Author: Behzad Moumbeini
# this script has these 4 outputs:
# 1- normalized data
# 2- group comparison plot
# 3- heatmap
# 4- scatter plot
# usage: group_comparison.py [-h] -G GENE
'''
put all the .RCC files in the same directory

'''
import seaborn as sns
import os

def files(path):
    folder = os.fsencode(path)
    filenames = []
    for file in os.listdir(folder):
        filename = os.fsdecode(file)
        if filename.endswith( ('.RCC') ):
            filenames.append(filename)
            filenames.sort()
    return filenames



import pandas as pd
from functools import reduce
from itertools import chain

def convert(filenames):
    dataframes = []
    # load all the dataframes in a list (dataframes)
    for filename in filenames:
        with open(filename) as f:
            total = f.readlines()
        skip_value = total.index('CodeClass,Name,Accession,Count\n')
        df = pd.read_csv(filename, skiprows=skip_value, skipfooter=5, sep=',')
        df = df.rename(columns={'Count': filename})
        dataframes.append(df)
    # merge the dataframes
    df_merged = reduce(lambda x,y: pd.merge(x,y, on=['CodeClass', 'Name', 'Accession'], how='outer'), dataframes)
    df_merged.to_csv('result15.csv')


import subprocess

def normalization(script):
    subprocess.call (["/usr/bin/Rscript", "--vanilla", script])
    return


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
def group_comparison(file1, file2, gene_name):
    df1 = pd.read_csv(file1, sep = '\t')
    df2 = pd.read_csv(file2, sep = '\t')
    df4 = df2.drop(['CodeClass', 'Accession'], axis=1)
    df5 = df4.set_index('Name').stack().reset_index().join(df1.set_index('sample'),on='level_1').rename(columns={'level_1': 'sample', 0: 'Value'})
    df5["expression log2"] = np.log(df5["Value"])
    df6 = df5.loc[df5['Name'] == gene_name]
    sns.set(rc={"axes.facecolor":"#e6e6e6",
            "axes.grid":False,
            'axes.labelsize':30,
            'figure.figsize':(20.0, 10.0),
            'xtick.labelsize':25,
            'ytick.labelsize':20})
    p = sns.violinplot(data=df6,
                   x = 'group',
                   y = 'expression log2',
                   notch=True).set_title(gene_name)
    plt.xticks(rotation=45)
    l = plt.xlabel('')
    plt.ylabel('expression log2')
    plt.savefig(f'{gene_name}.pdf')


import os
import pandas as pd
from heatmap import HEATMAP


def heatmap(normalzed_data_file):
    infile = normalzed_data_file
    prefix = os.path.splitext(infile)[0]
    df = pd.read_table(infile, index_col=0)
    df2 = np.log2(np.array(df.iloc[:, 3:]))
    df.iloc[:, 3:] = df2
    final = df
    final2 = final.iloc[:, 3:]
    HEATMAP(final2, prefix=prefix, scale="column")


#this piece generates the scatter plot n*n dimension
import pandas as pd
import numpy
import seaborn as sns
import seaborn
def scatter_plot(Normal_file, out_plot):
    df = pd.read_csv(Normal_file, sep="\t")
    df1 = df.iloc[:, 3:]
    columns = list(df1.columns.values)
    df2 = (numpy.log(df1))
    for i in columns:
        ax = seaborn.pairplot(df2, vars=columns, kind='reg')
        ax.savefig(out_plot)
    
def heatmap2(script):
    subprocess.call (["/usr/bin/Rscript", "--vanilla", script])
    return


def main():
    import argparse
    ap = argparse.ArgumentParser(description="")
    ap.add_argument('-G', '--gene', required=True)
    args = ap.parse_args()
    print(args)    
    comparison = group_comparison("metadata.txt", "Normalized_data.txt", args.gene)
    
    return 

if __name__ == "__main__":
    path = os.path.join(os.getcwd())
    file_list = files(path)
    print(convert(file_list))

    print(normalization("aval.r"))
    print(normalization("normalization.r"))

    comparison = main()
    print(comparison)
    print(heatmap("Normalized_data.txt"))
    print(heatmap2("heatmap.r"))
#    print(scatter_plot("Normalized_data.txt", "mp.pdf"))


