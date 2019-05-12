#!/usr/bin/env python2

from __future__ import division
import numpy as np
import scipy as sp
import pandas as pd
import argparse
import sys
import matplotlib
import matplotlib.pyplot as plt
import statsmodels.nonparametric.api as smnp

plt.style.use('seaborn-colorblind')

# Options for axis scaling
# 1. f1 only: same across all plots, normalise
# args.sharey = "yes"
# #args.sharey_all = "no"
# #args.sharey_feature = "yes"
# args.normalise = "no" : should make no difference

# 2. f1 only: max of each plot, normalise
# args.sharey = "no"
# #args.sharey_all = "no"
# #args.sharey_feature = "yes"
# args.normalise = "no" : should make no difference

# 3. f1 and f2: like 2, but max across the pair of features
# args.aspect = 0.9
# args.sharey = "no"
# args.sharey_all = "no"
# args.sharey_feature = "yes"
# args.normalise = "no" : means that, given the feature is acquired on trajectory in 2, when?

parser = argparse.ArgumentParser()
parser.add_argument("-f", required=False, default = None, type=str)
parser.add_argument("-f2", required=False, default = None, type=str)
parser.add_argument("-outfile", required=False, default = "_raw_data", type=str)
parser.add_argument("-labels", required=False, default = None, type=str)
parser.add_argument("-width", required = False, default = 3, type = float)
parser.add_argument("-aspect", required = False, default = 0.5, type = float)
parser.add_argument("-ordered", required = False, default = "horizontal", type = str)
parser.add_argument("-fontsize", required = False, default = 6, type = float)
parser.add_argument("-xevery", required = False, default = 2, type = int)
parser.add_argument("-kde_points", required = False, default = 100, type = int)
parser.add_argument("-normalise", required = False, default = "no", type = str)
parser.add_argument("-bar", required = False, default = "yes", type = str)
parser.add_argument("-kde", required = False, default = "yes", type = str)
parser.add_argument("-seperated", required = False, default = "yes", type = str)
parser.add_argument("-bar_alpha", required = False, default = 1, type = float)
parser.add_argument("-kde_alpha", required = False, default = 0.25, type = float)
parser.add_argument("-bw", required = False, default = 0.5, type = float)
parser.add_argument("-cut", required = False, default = 1, type = float)
parser.add_argument("-sharey_type", required = False, default = "all", type = str)
parser.add_argument("-sharey", required = False, default = "yes", type = str)
parser.add_argument("-sharey_all", required = False, default = "yes", type = str)
parser.add_argument("-sharey_feature", required = False, default = "yes", type = str)
parser.add_argument("-transition_data", required = False, default = None, type = str)
parser.add_argument("-transition_data2", required = False, default = None, type = str)
parser.add_argument("-leftover_frequency", required = False, default = "no", type = str)
parser.add_argument("-greyregion", required = False, default = "#efefef", type = str)
parser.add_argument("-haxis", required = False, default = "yes", type = str)
parser.add_argument("-pad", required = False, default = 0.1, type = float)
parser.add_argument("-out_type", required = False, default = "pdf", type = str)
parser.add_argument("-verbose", required = False, default = "yes", type = str)
args = parser.parse_args()

# See above for update
# Axes sharing cases (also matters what normalise is...)
# args.sharey
# args.sharey_all
# args.sharey_feature
#     1. Single data: max y across all axes: sharey == "yes"
#     2. Single data: max y for a given feature: sharey == "no"
#     3. Double data: max y across all axes: sharey_all == "yes"
#     4. Double data: max y across all features for each set: sharey_all == "no", sharey == "yes"
#     5. Double data: max y across both within same feature: sharey_all == "no", sharey == "no", sharey_feature == "yes"
#     6. Double data: max y for each data for each feature: sharey_all == "no", sharey == "no", sharey_feature == "no"
#     7. Double data (normed to 1): max y for each data for each feature: sharey_all == "no", sharey == "no", sharey_feature == "no"


# DEFAULT TO STYLE 3, if f2 != None, and default to same transition data
if args.f2 != None:
    args.aspect = 0.9
    args.normalise = "no"
    args.sharey_type = "feature"
    if args.transition_data2 == None:
        args.transition_data2 = args.transition_data

if args.sharey_type == "all":
    args.sharey = "yes"
    args.sharey_all = "yes"

if args.sharey_type == "feature":
    args.sharey = "no"
    args.sharey_all = "no"
    args.sharey_feature = "yes"


font = {'family':'sans-serif',
        'size':args.fontsize,
        'sans-serif':['Arial']}

matplotlib.rc('font', **font)

df = pd.read_csv(args.f, index_col=None, sep=",").astype(int)
L = np.amax(df.values)+1
if args.f2 != None:
    df2 = pd.read_csv(args.f2, index_col=None, sep=",").astype(int)
    L = max(np.amax(df.values), np.amax(df2.values))+1


# Methods to make PDFs and KDEs
def MakePdf(df, l=L):
    X = np.zeros((l,l))
    for i, el in enumerate(df.values):
        X[el[1], el[0]] += 1
    return X

def Order(pdf, l=L):
    steps = np.array(range(l))
    means = [(np.dot(steps, el), i) for i, el in enumerate(pdf)]
    return [el[1] for el in sorted(means)]

def GetLabels(file_name, l=L):
    if file_name!=None:
        return pd.read_csv(file_name, header=None, index_col=None).iloc[:,1]
    else:
        return [str(i+1) for i in range(L)]

def Kde(data, bw=args.bw, kernel="gau", gridsize=100., cut=args.cut, clip=(-np.inf, np.inf), cumulative=False):
    """Compute a univariate kernel density estimate using statsmodels."""
    fft = kernel == "gau"
    kde = smnp.KDEUnivariate(np.array([float(el) for el in data]))
    kde.fit(kernel, bw, fft, gridsize=gridsize, cut=cut, clip=clip)
    return kde

def MakeKde(df, points = args.kde_points, start = -args.cut*args.bw, l=L):
    end = -start + l - 1
    grid = np.linspace(start, end, points)
    kdes = []
    for i in range(l):
        data = df[df["feature"].astype(int)==i]["step"]
        if len(data) == 0:
            kdes.append(np.array([0 for el in grid]))
            continue
        kde = Kde(data)
        kde_grid = kde.evaluate(grid)
        kdes.append(kde_grid)
    width = (end-start)/points
    adjust = [np.sum(el) * width for el in kdes]
    kdes = [(el1/el2 if el2 > 0 else el1) for el1,el2 in zip(kdes, adjust)]
    return kdes, grid, width

def FindN(df):
    N=1
    for i in range(1, len(df)):
        if df.iloc[i,0] <= df.iloc[i-1,0]:
            N+=1
    return N

def Normalise(X, totals, N, normalise="yes'"):
    if normalise == "yes":
        return np.array([el*totals[i]/N for i, el in enumerate(X)])
    return np.array([el for i, el in enumerate(X)])

def MaxYs(X, Y):
    maxys = []
    for i in range(len(X)):
        maxy = max(np.amax(X[i,:]),np.amax(Y[i,:]))
        if np.isnan(maxy) or maxy == 0:
            maxy = 1e-9
        maxys.append(maxy)
    return maxys

def MakePdfComplete(df, normalise = args.normalise):
    counts = np.array(MakePdf(df))
    N = FindN(df)
    totals = np.array([np.sum(el) for el in counts])
    pdf = np.array([(el/totals[i] if totals[i] != 0 else el) for i, el in enumerate(counts)])
    if normalise == "yes":
        pdf = Normalise(pdf, totals, N, normalise)
    return pdf, totals, N


def MakePdfKdeYs(df, normalise = args.normalise, make_kdes=args.kde, make_bars=args.bar):
    pdf, totals, N = MakePdfComplete(df, normalise=normalise)
    kdes = None
    if make_kdes == "yes":
        kdes, grid, width = MakeKde(df)
        
    # Normalise by counts and get maxys
    pdf = Normalise(pdf, totals, N, normalise)
    if make_kdes == "yes":
        kdes = Normalise(kdes, totals, N, normalise=normalise)
        maxys = MaxYs(pdf, kdes)
    else:
        maxys = MaxYs(pdf, pdf)
    return pdf, kdes, maxys


def EarliestAndLatest(in_file=args.transition_data, l=L):
    out = [[0,l] for i in range(l)]
    if in_file != None:
        out = [[l,0] for i in range(l)]
        df = pd.read_csv(in_file, header=None, index_col=None, sep = ' ').astype(int)
        for i in range(0,len(df),2):
            ones_start = list(df.iloc[i,:]).count(1)
            ones_end   = list(df.iloc[i+1,:]).count(1)
            for j in range(len(df.iloc[i,:])):
                if df.iloc[i+1,j] != df.iloc[i,j]:
                    out[j][0] = ones_start if ones_start < out[j][0] else out[j][0]
                    out[j][1] = ones_end if ones_end+1 > out[j][1] else out[j][1]
    return out

# Make pdfs, kdes, maxys
pdf, kdes, maxys = MakePdfKdeYs(df)
ts = EarliestAndLatest(args.transition_data, l=L)
if args.f2 != None:
    pdf2, kdes2, maxys2 = MakePdfKdeYs(df2)
    ts2 = EarliestAndLatest(args.transition_data2, l=L)

# Make orders
order = Order(pdf) if args.ordered == "yes" else range(L)
if args.f2 != None:
    order2 = Order(pdf) if args.ordered == "yes" else range(L)
labels = GetLabels(file_name=args.labels, l=L)


if args.verbose == "yes":
    print "Fontsize:", args.fontsize
    print [np.sum(el) for el in pdf]
    print order
    print maxys, len(maxys), min(maxys), max(maxys)


# -------------------------------------------
# PLOTTING
# -------------------------------------------
fig, ax = plt.subplots(L, sharex=True)
fig.set_size_inches(args.width, args.aspect*args.width)
xs = np.arange(0, L, 1)
start = -args.bw*args.cut
end = -start + L - 1
if args.verbose == "yes":
    print start, end
grid = np.linspace(start, end, args.kde_points)
if args.f2 == None:
    if args.sharey == "yes":
        ylims = [(0, max(maxys)) for _ in range(L)]
    if args.sharey == "no":
        ylims = [(0, maxys[i]) for i in range(L)]
else:
    if args.sharey_all == "yes":
        max_all = max(max(maxys), max(maxys2))
        ylims = [(-max_all, max_all) for _ in range(L)]
    else:
        if args.sharey_all == "no" and args.sharey == "yes":
            ylims = [(-max(maxys2), max(maxys)) for _ in range(L)]
        else:
            if args.sharey_all == "no" and args.sharey == "no" and args.sharey_feature == "yes":
                ylims = [(-max(maxys[i],maxys2[i]),max(maxys[i],maxys2[i])) for i in range(L)]
            else:
                if args.sharey_all == "no" and args.sharey == "no" and args.sharey_feature == "no":
                    ylims = [(-maxys2[i],maxys[i]) for i in range(L)]

# Pad limits                    
pad_factor = 1/(1-args.pad)
ylims = [(ylims[i][0]*pad_factor, ylims[i][1]*pad_factor) for i in range(len(ylims))]

if args.verbose == "yes":
    for i, lim in enumerate(ylims):
        print lim, max(pdf[i,:])

for i in range(L):
    ax[i].bar(xs, pdf[order[i],:], alpha=args.bar_alpha, color="C0", lw=0)
    if args.kde == "yes":
        ax[i].fill_between(grid, 0, kdes[order[i]], alpha=args.kde_alpha, color="C0", lw=0)
    if args.f2 != None:
        ax[i].bar(xs, -pdf2[order[i],:], alpha=args.bar_alpha, color="C2", lw=0)
        if args.kde == "yes":
            ax[i].fill_between(grid, -kdes2[order[i]], 0, alpha=args.kde_alpha, color="C2", lw=0)
    if args.transition_data != None and args.transition_data2 == None:
        ax[i].axvspan(xmin=-.5, xmax=ts[order[i]][0]-.5, ymin=0, ymax=1,
                      color=args.greyregion, zorder=-1, lw=0)
        ax[i].axvspan(xmin=ts[order[i]][1]-.5, xmax=L-.5, ymin=0, ymax=1,
                      color=args.greyregion, zorder=-1, lw=0)
    if args.transition_data != None and args.transition_data2 != None:
        ax[i].axvspan(xmin=-.5, xmax=ts[order[i]][0]-.5, ymin=0.5, ymax=1,
                      color=args.greyregion, zorder=-1, lw=0)
        ax[i].axvspan(xmin=ts[order[i]][1]-.5, xmax=L-.5, ymin=0.5, ymax=1,
                      color=args.greyregion, zorder=-1, lw=0)
        ax[i].axvspan(xmin=-.5, xmax=ts2[order[i]][0]-.5, ymin=0, ymax=0.5,
                      color=args.greyregion, zorder=-1, lw=0)
        ax[i].axvspan(xmin=ts2[order[i]][1]-.5, xmax=L-.5, ymin=0, ymax=0.5,
                      color=args.greyregion, zorder=-1, lw=0)

# Adjust the subplot properties:
# Remove all axes ticks
# Remove all borders
for i in range(L):
    ax[i].spines['right'].set_visible(False)
    ax[i].spines['top'].set_visible(False)
    ax[i].spines['left'].set_visible(False)
    ax[i].spines['bottom'].set_visible(False)
    ax[i].set_yticks([])
    ax[i].set_yticklabels([])
    if args.f2 == None:
        ax[i].set_ylabel(labels[order[i]], rotation="horizontal", ha="right", va="baseline", y=0.0)
    else:
        ax[i].set_ylabel(labels[order[i]], rotation="horizontal", ha="right", va="center")

    ax[i].set_ylim(ylims[order[i]])
    
    # Horizontal axis
    if args.haxis == "yes":
        ax[i].axhline(0, xmin=0, xmax=L, lw=0.25, color = "dimgrey", zorder=-2)
    else:
        if args.haxis == "scaled":
            ax[i].axhline(0, xmin=0, xmax=L, lw=0.25*20/L, color = "dimgrey", zorder=-2)
            
    if i < L-1:
        ax[i].set_xticks([])
        ax[i].tick_params(axis=u'both', which=u'both',length=0)
    else:
        ax[i].spines['bottom'].set_visible(False)
        ax[i].tick_params(axis=u'both', which=u'both',length=0)
        ax[i].set_xlabel("Number of features acquired")
        ax[i].set_xticks(np.arange(-0.5,L+0.5,args.xevery))
        ax[i].set_xticklabels(range(0,L+1,args.xevery))

if args.seperated == "yes":
    fig.subplots_adjust(hspace=0.)
else:
    fig.subplots_adjust(hspace=-.25)

if args.out_type == "pdf":
    from matplotlib.backends.backend_pdf import PdfPages
    plot = PdfPages(args.outfile + ".pdf")
    plot.savefig(bbox_inches='tight', figure=fig)
    plot.close()
else:
    if args.out_type == "png":
        fig.savefig(args.outfile + ".png", bbox_inches='tight', dpi=600)

