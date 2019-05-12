#!/usr/bin/env python

from __future__ import division
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import argparse
from matplotlib.colors import LinearSegmentedColormap

parser = argparse.ArgumentParser()
parser.add_argument("-f", required=False, default = None, type=str)
parser.add_argument("-outfile", required=False, default = "_raw_data", type=str)
parser.add_argument("-labels", required=False, default = None, type=str)
parser.add_argument("-width", required = False, default = 2.48, type = float)
parser.add_argument("-aspect", required = False, default = 0.5, type = float)
parser.add_argument("-order", required = False, default = "horizontal", type = str)
parser.add_argument("-fontsize", required = False, default = 6, type = float)
parser.add_argument("-every", required = False, default = 1, type = int)
parser.add_argument("-xlabelevery", required = False, default = 5, type = int)
parser.add_argument("-row_every", required = False, default = 2, type = int)
parser.add_argument("-plot_type", required = False, default = "heatmap", type = str)
parser.add_argument("-xticks", required = False, default = "yes", type = str)
parser.add_argument("-xticklabels", required = False, default = "yes", type = str)
parser.add_argument("-xlinewidth", required = False, default = 0.2, type = float)
parser.add_argument("-ylinewidth", required = False, default = 0.2, type = float)
parser.add_argument("-outfile_type", required = False, default = "pdf", type = str)
parser.add_argument("-spines", required = False, default = "no", type = str)
parser.add_argument("-xlabel", required = False, default = "Samples", type = str)
args = parser.parse_args()

args.outfile = (args.f).replace(".txt", "_") + args.plot_type + "." + args.outfile_type

print args

def Size(width, height, s=1):
    return s * 34/max(width, height)


font = {'family':'sans-serif',
        'size':args.fontsize,
        'sans-serif':['Arial']}
matplotlib.rc('font', **font)

plt.style.use('seaborn-colorblind')

df = pd.read_csv(args.f, header=None, index_col=None, sep = " ")
if args.row_every == 2:
    df = df.iloc[1::args.row_every,:]
else:
    if args.row_every == 1:
        df = df.iloc[0::args.row_every,:]

labels = (pd.read_csv(args.labels, header=None, index_col=None)).iloc[:,1] \
            if args.labels != None \
            else [str(i+1) for i in range(df.shape[1])]

if args.order == "horizontal":
    df = df.T

cmap_name = "bw"
colors = ["white","black"]
cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=2)

fig, ax = plt.subplots()
fig.set_size_inches(args.width, args.width*args.aspect)

if args.plot_type == "spy":
    im = ax.spy(df.values, precision=0, marker = "o", ms=Size(df.shape[1], df.shape[0]),
                aspect = "auto", markeredgecolor='none')
if args.plot_type == "heatmap":
    im = ax.pcolormesh(df.values, cmap=cm, antialiased=False)
    
if args.order == "horizontal":
    if args.plot_type == "spy":
        ax.set_yticks(np.arange(0.0, df.shape[0], args.every))
        ax.set_yticklabels([(el if i%args.every == 0 else "") for i,el in enumerate(labels)])
        ax.set_xticks([])
        ax.set_xticklabels([])

        ax.xaxis.set_ticks_position('bottom')        
        for i in range(df.shape[0]):
            ax.axhline(i+0.5, xmin=0, xmax=df.shape[1], color = "lightgrey",
                      lw=Size(df.shape[1], df.shape[0], s=0.1))
        for i in range(df.shape[1]):
            ax.axvline(i+0.5, ymin=0, ymax=df.shape[0], color = "lightgrey",
                      lw=Size(df.shape[1], df.shape[0], s=0.1))
    else:
        ax.set_yticks(np.arange(0.5, df.shape[0], args.every))
        ax.set_yticklabels([(el if i%args.every == 0 else "") for i,el in enumerate(labels)])
        if args.xticks == "yes":
            ax.set_xticks(np.arange(0.5, df.shape[1], args.every))
            ax.set_xticklabels(["" for el in ax.get_xticks()])
        if args.xticklabels == "yes":
            a1 = [0.5]
            a2 = list(np.arange(args.xlabelevery-0.5, df.shape[1]+0.5, args.xlabelevery))
            ax.set_xticks(np.array(a1+a2))
            xticklabels = [(str(int(i)-1) \
                            if i>1 \
                            else str(int(i))) for i in (np.arange(1, df.shape[1]+1.1, args.xlabelevery))]
            ax.set_xticklabels(xticklabels)
        ax.xaxis.set_ticks_position('bottom')
        args.gap = 0.1/df.shape[1]
        xgap = args.xlinewidth
        ygap = args.ylinewidth
        for i in range(df.shape[0]+1):
            ax.fill_between([0,df.shape[1]+1], i-ygap, i+ygap, color="white", lw=0)
        for i in range(df.shape[1]+1):
            ax.fill_betweenx([0,df.shape[0]+1], i-xgap, i+xgap, color="white", lw=0)
        ax.set_xlim(-xgap, df.shape[1]+xgap)
        ax.set_ylim(-ygap, df.shape[0]+ygap)
        
    ax.set_ylabel("Features", fontsize=args.fontsize)
    ax.set_xlabel(args.xlabel, fontsize=args.fontsize, labelpad=0)
    ax.invert_yaxis()
    ax.grid(False)

if args.outfile_type == "pdf":    
    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages(args.outfile.replace(".png", ".pdf"))
    pdf.savefig(bbox_inches = 'tight',figure=fig)
    pdf.close()
else:
    fig.savefig(args.outfile, bbox_inches = 'tight', dpi=600, transparent=True)