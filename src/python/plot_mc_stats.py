#!/usr/bin/env python

from __future__ import division
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse

plt.style.use('ggplot')
parser = argparse.ArgumentParser()
parser.add_argument('-f', required=False, default="stats.csv", type=str)
parser.add_argument('-b', required=False, type=int, default = 50000)
parser.add_argument('-last', required=False, type=int, default=500)
args = parser.parse_args()

print args

file1 = args.f

# Load and label data
df1 = pd.read_csv(file1,index_col = None)
mc_step1 = df1['mcmc_step'] if 'mcmc_step' in df1 else df1['smc_step']
log_lik1 = df1['log_likelihood']
accepts1 = df1['accepts']
accepts_lower1 = df1['accepts_lower']

start = list(mc_step1).index(args.b)
last = int(args.last/mc_step1[1])

# Calculate acceptance over "last" number of time steps in "stats.csv"
temp_ = []
accepts1_l = list(df1["accepts"])
step = range(0,len(df1)+1)
for i in range(len(df1)):
    if i > last+1:
        temp_.append((accepts1_l[i]*step[i+1] - accepts1_l[i-last+1]*(step[i-last+1]))/(step[i+1] - (step[i-last+1])))
    else:
        temp_.append(np.nan)
accepts3 = temp_

temp_ = []
accepts_l = list(df1["accepts"])
accepts0_l = list(df1["accepts0"])
accepts_lower_l = list(df1["accepts_lower"])
accepts0_lower_l = list(df1["accepts0_lower"])
step = range(0,len(df1)+1)
accepts_3 = []
accepts0_3 = []
accepts_lower_3 = []
accepts0_lower_3 = []
div_3 = []
div0_3 = []
data_lists = [accepts_l, accepts0_l, accepts_lower_l, accepts0_lower_l]
new_data_lists = [accepts_3, accepts0_3, accepts_lower_3, accepts0_lower_3]
for i in range(len(df1)):
    for k, l in enumerate(data_lists):
        if i > last+1:
            new_data_lists[k].append((l[i]*step[i+1] - l[i-last+1]*(step[i-last+1]))/(step[i+1] - (step[i-last+1])))
        else:
            new_data_lists[k].append(np.nan)
    try:
        div_3.append(accepts_lower_3[i]/new_data_lists[0][i])
    except:
        div_3.append(np.nan)
    try:
        div0_3.append(new_data_lists[3][i]/new_data_lists[1][i])
    except:
        div0_3.append(np.nan)
        
def PlotFigure(mc_step1, log_lik1, div_3, div0_3, new_data_lists, start):
    plt.clf()
    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    ax1.plot(mc_step1, log_lik1, label="log likelihood")
    ax2.plot(mc_step1[start:], div_3[start:], label="MH decrease accept: " + str(last), color = 'C1')
    ax2.plot(mc_step1[start:], div0_3[start:], label="MH APM decrease accept: " + str(last), color = 'C4')

    ax2.plot(mc_step1[start:], new_data_lists[0][start:], label="Accept proportion: " + str(last), color='C2')
    ax2.plot(mc_step1[start:], new_data_lists[1][start:], label="MH APM accept proportion: " + str(last), color='C3')
    max_y = max(log_lik1[start:])
    min_y = min(log_lik1[start:])
    tolerance=1.1
    ax1.set_ylim(min_y*tolerance, max_y/tolerance)
    ax2.set_ylim(0,1)
    ax1.set_xlim(min(mc_step1[start:]),max(mc_step1[start:]))
    ax2.set_xlim(min(mc_step1[start:]),max(mc_step1[start:]))
    ax1.set_ylabel("Log likelihood")
    if start == 0:
        ax1.set_xlabel("MC step")
        ax2.set_xlabel("MC step")
    else:
        ax1.set_xlabel("MC step after burn-in")
        ax2.set_xlabel("MC step after burn-in")
    ax2.set_ylabel("Proportion")
    ax2.legend(loc="upper left", fontsize='x-small', frameon=True)
    fig.subplots_adjust(left=None, bottom=None, right=None, top=None,
                    wspace=0.3, hspace=None)
    fig.savefig("stats-" + str(start) +".png", transparent=True, bbox_icnhes='tight')

PlotFigure(mc_step1, log_lik1, div_3, div0_3, new_data_lists, 0)
PlotFigure([el - args.b for el in list(mc_step1)], log_lik1, div_3, div0_3, new_data_lists, start)

print "-"*100
print "Acceptance proportion at end over last " + str(args.last) + " iterations               : " +str(round(new_data_lists[0][-1],3))
print "Increasing likelihood proportion end over last " + str(args.last) + " iterations       : " +str(round(div_3[-1],2))
print "Acceptance proportion at end over last (APM) " + str(args.last) + " iterations         : " +str(round(new_data_lists[1][-1],3))
print "Increasing likelihood proportion end over last (APM) " + str(args.last) + " iterations : " +str(round(div0_3[-1],2))
print "-"*100