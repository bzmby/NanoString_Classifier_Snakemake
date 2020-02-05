#!/usr/bin/env python
#Author: Behzad Moumbeini
# this script has these  outputs:
# 1- vsc file with 2 QC criteria 
# 2- plots with 2 lines
# usage: python3 QC.py
'''
put all the .RCC files in the same directory

'''


import matplotlib.pyplot as plt
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from numpy import array
import numpy
import csv
import os



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


def QC(file):
    with open(file) as f:
        for line in f.read().split():
            segments = line.split(",")
            if segments[0] == "FovCount":
                FC = segments[1]
            elif segments[0] == "FovCounted":
                FCed = segments[1]
                criteria = int(FC)/int(FCed)
            elif segments[0] == "BindingDensity":
                BD = segments[1]
    return [criteria, BD]


def QC2(fil):
    Di = {}
    Pos = []
    Neg = []
    with open(fil) as fl:
        for line in fl.read().split():
            seg = line.split(",")
            if seg[0] == "Positive":
                value = seg[3]
                Pos.append(value)
            elif seg[0] == "Negative":
                value = seg[3]
                Neg.append(value)
            Di["Positive"] = Pos
            Di["Negative"] = Neg
    return Di


def pre_plot(file):
    num = dict()
    Pos = []
    with open(file, 'r') as opened:
        for row in opened:
            a = ' '.join(row.splitlines()).split(',')
            if a[0] == 'Positive':
                num[a[1][4]] = int(a[3])
    return num



def QC_seg(file):
    main = {}
    for i in file:
        main[i] = QC(i)
    return main


def QC_plot(file):
    main2 = {}
    for i in file:
        main2[i] = QC2(i)
    return main2



def plot(file_list):
    nums = []
    Dis = []
    for i in range(len(file_list)):
        nums.append(pre_plot(file_list[i]))
        Dis.append(QC2(file_list[i]))


    plt.style.use('default')
    sns.set()
    sns.set_style('whitegrid')
    sns.set_palette('Set1')
    # concentration = np.log2(np.array([0.125, 0.5, 2, 8, 32, 128]))
    concentration = np.log2(np.array([128, 32, 8, 2, 0.5, 0.125]))
    Pos = [ [i[e] for e in sorted(i.keys())] for i in nums ]
    Pos_log = numpy.log2(numpy.array(Pos))
    Neg = [np.mean(np.log2(np.array(i['Negative'],np.float64)))+2*np.std(np.log2(np.array(i['Negative'],np.float64))) for i in Dis] 
    fig = plt.figure()
    for i in range(len(Neg)):
        plt.plot(concentration, Pos_log[i], label=file_list[i])
        plt.axhline(y=Neg[i], color='b', linestyle='-')
        plt.legend()
        plt.xlabel("log2 concentration")
        plt.ylabel("log2 raw counts")
        plt.ylim(0, 40)
        fig = plt.gcf()
        # plt.show()
        fig.savefig(f'{i}.pdf')
        plt.close('all')



def Quality_measures(path):
    file = files(path)
    qc = []
    for i in file:
        qc.append(QC(i))


    new_dict = dict(zip(file, qc))
    return new_dict



if __name__ == "__main__":

    path = os.path.join(os.getcwd())

    file_list = files(path)

    figure1 = plot(file_list)
    print(figure1)

    QC_measure = Quality_measures(path)

    with open('quality_control.csv', 'w') as f:
        header = ["file_name", "imaging_QC", "Binding_Density"]
        w = csv.writer(f)
        w.writerow(header)
        for key, lst in QC_measure.items():
            w.writerow([key] + lst)


