#!/usr/bin/env python

from __future__ import division
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
import argparse
from collections import OrderedDict

def smooth(x, window_len=11, window='hanning'):
    # From: https://scipy-cookbook.readthedocs.io/items/SignalSmooth.html
    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."
    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."
    if window_len<3:
        return x
    
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"

    s=np.r_[2*x[0]-x[window_len:1:-1], x, 2*x[-1]-x[-1:-window_len:-1]]

    if window == 'flat': #moving average                                                                                                    
        w = np.ones(window_len,'d')
    else:
        w = getattr(np, window)(window_len)
    y = np.convolve(w/w.sum(), s, mode='same')
    return y[window_len-1:-window_len+1]

def MinAIC(parameters, variable):
    grouped_by_parameter = OrderedDict()
    for i in range(len(parameters)):
        p = int(parameters[i])
        v = float(variable[i])
        if p in grouped_by_parameter:
            grouped_by_parameter[p].append(v)
        else:
            grouped_by_parameter[p] = [v]

    mins = []
    params= []
    for key, value in grouped_by_parameter.items():
        params.append(int(key))
        minimum = np.min(value)
        mins.append(minimum)
    return mins, params


plt.style.use('seaborn-colorblind')
cwd = os.getcwd()
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--f", type=str, required=True)
parser.add_argument("-f2", type=str, required=False)
parser.add_argument("-outfile", type=str, required=False, default="_test")
parser.add_argument("-loc", type=str, required=False)
parser.add_argument("-width", type=float, required=False, default=3.4)
parser.add_argument("-height", type=float, required=False, default=2)
parser.add_argument("-aspect", type=float, required=False, default=None)
parser.add_argument("-fontsize", type=float, required=False, default=6)
parser.add_argument("-minvalue", type=float, required=False, default=None)
parser.add_argument("-pdf", type=str, required=False, default="yes")
parser.add_argument("-smooth_so", type=int, required=False, default=None)
parser.add_argument("-legend_columns", type=int, required=False, default=2)
args = parser.parse_args()

print args

font = {'family':'sans-serif',
        'size':args.fontsize,
        'sans-serif':['Arial']}
matplotlib.rc('font', **font)

file1 = args.f
df_so = pd.read_csv(file1,index_col = None)
params_so   = df_so['parameters']
log_lik_so  = df_so['log_likelihood_aic']
aic_so      = df_so['aic_minimum']

aic_null = aic_so.iloc[params_so.idxmin()]
ll_null  = -aic_null/2

min_el_so     = aic_so.idxmin()
min_k_so      = params_so.iloc[min_el_so]
min_aic_so    = aic_so.iloc[min_el_so]
min_aic_ll_so = log_lik_so.iloc[min_el_so]
min_params_so = params_so.iloc[min_el_so]

mean_aic_so, params_so = MinAIC(list(params_so), list(aic_so))
if args.smooth_so != None:
    smooth_so = smooth(np.array(mean_aic_so), window_len=args.smooth_so)
else:
    smooth_so = np.array(mean_aic_so)

if args.f2 != None:
    file2 = args.f2
    df_fo   = pd.read_csv(file2, index_col = None)
    params_fo  = df_fo['parameters']
    log_lik_fo = df_fo['log_likelihood_aic']
    aic_fo     = df_fo['aic_minimum']

    min_el_fo     = aic_fo.idxmin()
    min_k_fo      = params_fo.iloc[min_el_fo]
    min_aic_fo    = aic_fo.iloc[min_el_fo]
    min_aic_ll_fo = log_lik_fo.iloc[min_el_fo]
    min_params_fo = params_fo.iloc[min_el_fo]

    mean_aic_fo, params_fo = MinAIC(list(params_fo), list(aic_fo))
    smooth_fo = np.array(mean_aic_fo)

    
print "Zero order", 0, aic_null, ll_null
if args.f2!= None:
    print "First order:", min_k_fo, min_aic_fo, min_aic_ll_fo
print "Second order:", min_k_so, min_aic_so, min_aic_ll_so

fig = plt.figure()
fig.set_size_inches(args.width, args.height if args.aspect == None else args.width * args.aspect)
ax1 = fig.add_subplot(111)
ax1.set_xlabel("Parameters", labelpad=-.1)
ax1.set_ylabel("AIC")

alpha_custom = 0.5
alpha_point = 0.75
ms = 15
zm = 10
po_zo_aic = ax1.scatter(0, aic_null, color='C0', label="Zero order", s=ms, marker="o", zorder=zm, alpha=alpha_point, lw=0)
if args.f2 != None:
    po_fo_aic = ax1.scatter(min_k_fo, min_aic_fo, color='C1', label="First order regularised", s=ms, marker="o", zorder=zm, alpha=alpha_point, lw=0)
po_so_aic = ax1.scatter(min_k_so, min_aic_so, color='C2', label="Second order regularised", s=ms, marker="o", zorder=zm, alpha=alpha_point, lw=0)

if args.f2 != None:
    p_fo_aic = ax1.plot(params_fo, smooth_fo, color="C1", linestyle=":", label='First order regularisation', zorder=0)
p_so_aic = ax1.plot(params_so, smooth_so, color="C2", linestyle=":", label='Second order regularisation', zorder=0)

if args.f2 != None:
    ys = [aic_null, min_aic_fo, min_aic_so, max(smooth_so), min(smooth_so), max(smooth_fo), min(smooth_fo)]
else:
    ys = [aic_null, min_aic_so, max(smooth_so), min(smooth_so)]
    
if args.minvalue != None:
    ys.append(args.minvalue)

max_y1 = max(ys)
min_y1 = min(ys)
extra = 0.2
bbox_extra = 0.125
space = extra*(max_y1-min_y1)
ax1.set_ylim(min_y1 - space, max_y1+space)

sp = fig.subplotpars
bb=[sp.left, sp.top+bbox_extra*(sp.top-sp.bottom), sp.right-sp.left, bbox_extra*(sp.top-sp.bottom) ]

handles, labels = ax1.get_legend_handles_labels()
handles = handles[2:] + handles[:2]
labels  = labels[2:] + labels[:2]

if args.legend_columns == 1:
    leg = ax1.legend(handles, labels, bbox_to_anchor=bb, mode="expand", borderaxespad=0,
                     bbox_transform=fig.transFigure, fontsize=args.fontsize, ncol=1)
else:
    if args.legend_columns == 2:
        leg = ax1.legend(handles, labels, bbox_to_anchor=bb, mode="expand", borderaxespad=0,
                         bbox_transform=fig.transFigure, fontsize="x-small", ncol=2)

out_file = args.outfile
if args.pdf == "yes":
    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages(out_file + ".pdf")
    pdf.savefig(bbox_inches = 'tight', figure=fig)
    pdf.close()
if args.pdf == "no":
    fig.savefig(out_file + ".png", dpi=600, bbox_inches='tight', transparent=True)

