#!/usr/bin/env python

from __future__ import division
import matplotlib as mpl
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import networkx as nx
from networkx.drawing.nx_agraph import write_dot, graphviz_layout
import numpy as np
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f", required=False, default = None, type=str)
parser.add_argument("-transitions", required=False, default = None)
parser.add_argument("-bar", required=False, default = "no")
parser.add_argument("-heatmap", required=False, default = "yes")
parser.add_argument("-heatmap_max", required=False, default = None, type=float)
parser.add_argument("-data_type", required=False, default = "match-data", type=str)
parser.add_argument("-prob_type", required=False, default = "joint", type=str)
parser.add_argument("-layout_type", required=False, default = "dot", type=str)
parser.add_argument("-labels", required=False, default = "labels.csv", type=str)
parser.add_argument("-include_self_edges", required=False, default = None)
parser.add_argument("-fontsize", required=False, type = float, default = 4.)
parser.add_argument("-labelsfontsize", required=False, type = float, default = 2.)
parser.add_argument("-xtick_rotation", required=False, type = float, default = 90.)
parser.add_argument("-outfile_type", required=False, default = "pdf", type=str)
parser.add_argument("-width", required=False, default = 3., type=float)
parser.add_argument("-aspect", required=False, default = 0.5, type=float)
parser.add_argument("-s", required=False, default = 100, type=float)
parser.add_argument("-yevery", required=False, default = 1, type=int)
parser.add_argument("-outfile", required=False, default = None)
parser.add_argument("-node_size", required=False, default = 100., type=float)
parser.add_argument("-edge_width", required=False, default = 2.5, type=float)
parser.add_argument("-edge_width_exponent", required=False, default = 1, type=float)
parser.add_argument("-graph_fontsize", required=False, default = 8., type=float)
parser.add_argument("-connection_style", required=False, default = "arc3,rad=-0.2", type=str)
parser.add_argument("-any_time", required=False, default = 0, type=int)
parser.add_argument("-end", required=False, default = "no", type=str)
args = parser.parse_args()

print args

def LoadPi(in_file, N=100, rand=0):
    pis = []
    with open(in_file, 'r') as f:
        counter = 0
        for line_raw in f:
            if counter >= N:
                break
            line_raw2 = line_raw.rstrip("\r\n").split("\t")[1:]
            pi = []
            for line in line_raw2:
                pi.append([float(el) for el in line.split(" ")])
            pis.append(pi)
            counter += 1
    return pis
    
def MakeGraph(df, df2, L, data_type="zero-one"):
    G = nx.DiGraph()
    nt = np.sum(df["weight"])/L
    for i in range(len(df)):
        G.add_edge(df.iloc[i, 0],
                   df.iloc[i, 1],
                   weight=float(df.iloc[i, 2]),
                   inv_lweight=-np.log(float(df.iloc[i, 2]/nt)))
    
    for node in G.nodes():
        tweight = 0
        for edge in G.out_edges(node):
            tweight += G.edges[edge]["weight"]
        for edge in G.out_edges(node):
            G.edges[edge]["cprob"] = G.edges[edge]["weight"]/tweight
            G.edges[edge]["inv_lcprob"] = -np.log(G.edges[edge]["weight"]/tweight)

    for i in range(len(df2)):
        G.nodes()[df2.iloc[i,0]]["weight"] = df2.iloc[i,1]
    return G, nt

def BinFromInt(s, L):
    return list(reversed([int(x) for x in ('{:0' + str(L) + 'b}').format(s)]))

def Diff(s1, s2):
    for i, el in enumerate(s1):
        if el != s2[i]:
            return i
        
def MakeFeatureAdjacencyJoint(df, L, any_time=args.any_time, verbose=1):
    feature_transitions = [[0 for j in range(L+2)] for i in range(L+2)]
    count = 0
    for i in range(len(df)):
        traj = df.iloc[i,0]
        traj = [el for el in traj.split("-")[1:]]
        previous = 0
        if any_time == 0:
            for j, el in enumerate(traj):
                try:
                    el1 = int(el)+1
                except:
                    if str(el) != "end":
                        print el
                    el1 = L+1
                if previous == L+1:
                    print traj
                feature_transitions[previous][el1] += 1
                previous = el1
                count += 1
        else:
            temp = [0]
            for j, el in enumerate(traj):
                try:
                    el1 = int(el)+1
                except:
                    if str(el) != "end":
                        print el, type(el)
                    el1 = L+1
                for elp in temp:
                    feature_transitions[elp][el1] += 1
                temp.append(el1)
            count += 1
    for i, eli in enumerate(feature_transitions):
        for j, elj in enumerate(eli):
            feature_transitions[i][j] /= count
    return feature_transitions


def MakeColors(weights, minima_=-20, maxima_=20, colormap="bwr_r"):
    import matplotlib.cm as cm
    minima = minima_
    maxima = maxima_
    norm = mpl.colors.Normalize(vmin=minima, vmax=maxima)
    mapper = cm.ScalarMappable(norm=norm, cmap=plt.get_cmap(colormap))
    colors = []
    for v in weights:
        colors.append(mapper.to_rgba(v))
    return colors

def RemoveAxes(ax):
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.get_yaxis().set_ticks([])
    ax.get_yaxis().set_visible(False)
    ax.get_xaxis().set_ticks([])
    ax.get_xaxis().set_visible(False)


def PlotGraph(G, pos, ax, constyle, node_size, edge_width, font_size,
              edge_width_exponent=args.edge_width_exponent, colors=None,
              special=0):
    color = "C5"
    sign = lambda x: x and (1, -1)[x<0]
    widths = [float(G.edges()[edge]["weight"]) for edge in G.edges()]
    widths = [(abs(el)/max([abs(eli) for eli in widths]))**edge_width_exponent * edge_width * sign(el) for el in widths]
    weights = [(abs(el)/max([abs(eli) for eli in widths]))**edge_width_exponent * sign(el) for el in widths]
    nodes  = nx.draw_networkx_nodes(G, pos=pos, nodelist=G.nodes(),
                                    ax=ax,
                                    node_color = color,
                                    node_size = node_size)

    if colors == None:
        edges = nx.draw_networkx_edges(G, pos=pos, ax=ax,
                                       edgelist=G.edges(),
                                       edge_width=widths,
                                       arrowstyle="->",
                                       node_size = node_size,
                                       arrowsize= edge_width/0.25)

        for i, edge in enumerate(edges):
            edge.set_connectionstyle(constyle)
            edge.set_linewidth(abs(widths[i]*edge_width/max(widths)))
            edge.set_alpha(abs(widths[i]/max(widths)))
    else:
        colors = MakeColors([el for el in weights], minima_=-1, maxima_=1, colormap="RdBu")
        edges = nx.draw_networkx_edges(G, pos=pos, ax=ax,
                                       edgelist=G.edges(),
                                       edge_width=widths,
                                       arrowstyle="->",
                                       node_size = node_size,
                                       arrowsize= edge_width/0.25)
        
        for i, edge in enumerate(edges):
            edge.set_connectionstyle(constyle)
            edge.set_linewidth(abs(widths[i]*edge_width/max(widths)))
            edge.set_color(colors[i])
            edge.set_alpha((widths[i]+1)/2)
            
    labels = nx.draw_networkx_labels(G, pos=pos, ax=ax, font_size=font_size)
    RemoveAxes(ax)
    return

prefix = ""
args.outfile = prefix + "forwards-feature_graph_" + args.data_type + "_" + args.prob_type + "_" + args.layout_type
args.transitions = prefix + "forwards_list-pord-trajectory-" + args.data_type + ".txt"

font = {'family':'sans-serif',
        'size':args.fontsize,
        'sans-serif':['Arial']}    
matplotlib.rc('font', **font)

state_labels = [str(el) for el in list(pd.read_csv(args.labels, header=None, index_col = None).iloc[:,1])]
state_labels = ["root"] + state_labels
if args.end == "yes":
    state_labels = state_labels + ["end"]
    
n_ = 0
pis = LoadPi(args.f, N=100)
L = len(pis[0][0])

df = pd.read_csv(args.transitions, index_col=None, header = None)

feature_transitions = MakeFeatureAdjacencyJoint(df, L)
if args.end != "yes":
    feature_transitions = [el[:-1] for i, el in enumerate(feature_transitions[:-1])]
if args.prob_type == "conditional":
    feature_transitions = [(el/np.sum(el) if np.sum(el) > 0 else el) for el in feature_transitions]

G = nx.DiGraph()
node_sizes_d = {}
for i, label in enumerate(state_labels):
    node_sizes_d[label] = 1
for i, label in enumerate(state_labels):
     G.add_node(label, ns=node_sizes_d[label])

for i, el1 in enumerate(feature_transitions):
    for j, el2 in enumerate(el1):
        if abs(el2-0.0) > 1e-10:
            G.add_edge(state_labels[i], state_labels[j], weight=el2, invweight=abs(-np.log(el2)))

s = args.s
connection_style = args.connection_style
node_sizes = [G.nodes[node]["ns"] for node in G.nodes()]
node_sizes = [s*el/max(node_sizes) for el in node_sizes]
weights = [G.edges[edge]['weight'] for edge in G.edges()]

max_weight = max([el for el in weights])
colors = MakeColors(weights, minima_=0, maxima_=1, colormap="Blues")
widths = [((el/max_weight)**args.edge_width_exponent * args.edge_width if el > 0 else 0) for el in weights]

if args.layout_type == "circular":
    pos = {}
    for node in G.nodes():
        try:
            el = int(node)
            angle = 2*np.pi * el/(L+1)
        except:
            if node == "root":
                angle = 0
            if node == "end":
                pos[node] = (0, 0)
                continue
            if node != "root" and node != "end":
                pos = nx.circular_layout(G)
                break
        pos[node] = ((L+1)*np.sin(angle), (L+1)*np.cos(angle))
    
if args.layout_type == "dot":
    pos = graphviz_layout(G, prog='dot')

plt.clf()
plt.style.use('seaborn-colorblind')
fig = plt.figure()
fig.set_size_inches(args.width,args.width)
ax2 = fig.add_subplot(111)
ax2.set_aspect(1)
PlotGraph(G=G,
          ax=ax2,
          constyle=connection_style,
          pos=pos,
          node_size=args.node_size,
          edge_width=args.edge_width,
          font_size=args.graph_fontsize)    
if args.outfile_type == "pdf":
    pdf = PdfPages(args.outfile + "_graph.pdf")
    pdf.savefig(bbox_inches = 'tight',figure=fig)
    pdf.close()
if args.outfile_type == "png":
    fig.savefig(args.outfile + "_graph.png", bbox_inches='tight', dpi=600, transparent=True)

    
plt.clf()
plt.style.use('seaborn-colorblind')
fig = plt.figure(frameon=False)
fig.set_size_inches(args.width,args.width)
gs = gridspec.GridSpec(1, 1)

ax = fig.add_subplot(gs[0, 0])

if args.heatmap == "yes" and args.bar != "yes":
    cmap = plt.get_cmap('viridis')
    if args.heatmap_max == None:
        _vmax = np.amax(np.array(feature_transitions))
    else:
        _vmax = args.heatmap_max
    im = ax.pcolormesh(np.array(feature_transitions), cmap=cmap,
                       vmin=0, vmax=_vmax)
    ax.invert_yaxis()
    fig.colorbar(im, ax=ax)
    ax.xaxis.set_ticks_position("top")
    ax.xaxis.set_label_position("top")
    ax.set_xticks(np.arange(0.5,len(state_labels)+0.5,1))
    ax.set_xticklabels([label for label in state_labels], rotation=args.xtick_rotation)
    ax.set_yticks(np.arange(0.5,len(state_labels)+0.5,1))
    ax.set_yticklabels([label for label in state_labels])
    ax.tick_params(axis='x', which='both', bottom=False, top=False)
    ax.tick_params(axis='y', which='both', left=False, right=False)
    
    if args.any_time == 0:
        if args.prob_type == "joint" and args.any_time == 0:
            ax.set_xlabel("Next acquisition")
        if args.prob_type == "conditional" and args.any_time == 0:
            ax.set_xlabel("Next conditional acquisition")
        ax.set_ylabel("Previous acquisition")
    else:
        ax.set_xlabel("Acquisition of Y given X previously acquired")
        ax.set_ylabel("Previous acquisition of X")    
    ax.set_aspect("equal")
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    
if args.outfile_type == "pdf":
    pdf = PdfPages(args.outfile + "_adjacency.pdf")
    pdf.savefig(bbox_inches = 'tight',figure=fig)
    pdf.close()
if args.outfile_type == "png":
    fig.savefig(args.outfile + "_adjacency.png", bbox_inches='tight', dpi=600, transparent=True)

try:
    write_dot(G, args.outfile + "_graph.gv")
except:
    pass
