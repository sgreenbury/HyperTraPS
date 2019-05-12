#!/usr/bin/env python

from __future__ import division
import matplotlib as mpl
import matplotlib
from matplotlib.patches import Ellipse
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import numpy as np
import pandas as pd
import argparse
import os
import operator as op
from collections import OrderedDict
np.random.seed(0)

parser = argparse.ArgumentParser()
parser.add_argument("-f", required=True, default = None)
parser.add_argument("-svg", required=False, default="no")
parser.add_argument("-edge_scale_exponent", required=False, default=0.5, type = float)
parser.add_argument("-edge_scale", required=False, default=2.5, type = float)
parser.add_argument("-edge_color", required=False, default="dimgrey", type = str)
parser.add_argument("-edge_alpha", required=False, default=0.5, type = float)
parser.add_argument("-show_edge_alpha", required=False, default="yes", type = str)
parser.add_argument("-node_scale", required=False, default=1.0, type = float)
parser.add_argument("-node_scale_exponent", required=False, default=1.0, type = float)
parser.add_argument("-node_alpha", required = False, default = 1, type = float)
parser.add_argument("-write_transitions", required=False, default="no")
parser.add_argument("-xaxis", required=False, default="yes")
parser.add_argument("-xevery", required=False, default=2, type = int)
parser.add_argument("-direction_weight", required=False, default=1, type = float)
parser.add_argument("-tlimit", required=False, default=int(1e8), type = int)
parser.add_argument("-labels", required=False, default=None)
parser.add_argument("-extra_labels", required=False, default=None)
parser.add_argument("-extra_labels_routes", required=False, default=None)
parser.add_argument("-fontsize", required=False, default=6, type=float)
parser.add_argument("-labels_fontsize", required=False, default=5, type=float)
parser.add_argument("-outfile_graph", required=False, default=None, type=str)
parser.add_argument("-edge_labels", required=False, default="False", type=str)
parser.add_argument("-gamma", required=False, default=1.0, type=float)
parser.add_argument("-gap_type", required=False, default="hamming", type=str)
parser.add_argument("-probable_paths", required = False, default = None, type = float)
parser.add_argument("-node_normalise", required = False, default = 1.0, type = float)
parser.add_argument("-edge_normalise", required = False, default = 1.0, type = float)
parser.add_argument("-layout_type", required = False, default = "noverlap", type = str)
parser.add_argument("-noverlap_space", required = False, default = 0.1, type = float)
parser.add_argument("-bbox", required = False, default = "tight", type = str)
parser.add_argument("-aspect", required = False, default = 0.75, type = float)
parser.add_argument("-width", required = False, default = 3, type = float)
parser.add_argument("-axlim", required = False, default = 1.3, type = float)
parser.add_argument("-max_width", required = False, default = 3, type = float)
parser.add_argument("-grid", required = False, default = "no", type = str)
parser.add_argument("-node_type", required = False, default = "int", type = str)
parser.add_argument("-out_type", required = False, default = "pdf", type = str)
parser.add_argument("-label_type", required = False, default = "greedy_data", type = str)
parser.add_argument("-xlabel", required = False, default = "Number of features acquired", type = str)
parser.add_argument("-show_se", required = False, default = "bracket", type = str)
parser.add_argument("-label_bbox", required = False, default = "yes", type = str)
parser.add_argument("-alternate_label", required = False, default = "no", type = str)
parser.add_argument("-radial_label_distance", required = False, default = -1, type = float)
parser.add_argument("-bbox_alpha", required = False, default = 0.75, type = float)
parser.add_argument("-extra", required = False, default = "no", type = str)
parser.add_argument("-label_colour_type", required = False, default = "gradient", type = str)
parser.add_argument("-colormap", required = False, default = "plasma", type = str)
parser.add_argument("-transition_data", required = False, default = None, type = str)
parser.add_argument("-arb", required = False, default = None, type = str)
args = parser.parse_args()

print args

# Methods

def BinFromInt(i, L):
    return list(reversed([int(x) for x in ('{:0' + str(L) + 'b}').format(i)]))

def ConvertToInt(v):
    count = 0
    for i, el in enumerate(v):
        count += 2**i if el == 1 else 0
    return count

def StateToString(state):
    return ",".join([str(el) for el in state])
def StringToState(state):
    return [int(el) for el in state.split(",")]

def MakeRoutes(df):
    routes = []
    for i in range(len(df)):
        if i > 0:
            if df.iloc[i,0] == 0 or (df.iloc[i,0] <= df.iloc[i-1,0]):
                if i > 0:
                    routes.append(route)
                    route = []
                if len(routes) > args.tlimit:
                    break
        else:
            route = []
        route.append(df.iloc[i,:]["feature"])
    return routes

def MakeTransitions(routes, uts=None, node_type = "ints"):
    transitions = {}
    if uts == None:
        uts = range(len(routes))
    for i in uts:
        route = routes[i]
        state = [0]*L
        for el in route:
            if node_type == "int":
                source = ConvertToInt(state)
            else:
                source = StateToString(state)
            state[el] = 1
            if node_type == "int":
                target = ConvertToInt(state)
            else:
                target = StateToString(state)
            pair = tuple([source, target])
            if pair in transitions:
                transitions[pair] += 1
            else:
                transitions[pair] = 1
    return transitions

def WriteOutTranstions(trans, ofile):
    out = []
    for key, value in trans.iteritems():
        el = [key[0], key[1], value]
        out.append(el)    
    df_out = pd.DataFrame(out, columns=["from", "to", "weight"], index=None)
    df_out.to_csv(ofile, index=None, sep=" ")

def WriteOutRoutes(rout, ofile):
    df_out_routes = pd.DataFrame(rout, columns=None, index=None)
    df_out_routes.to_csv(outfile.replace("edge-list"+node_tag, "routes"), index = None, header=None, sep = " ")

def ncr(n, r):
    r = min(r, n-r)
    numer = reduce(op.mul, xrange(n, n-r, -1), 1)
    denom = reduce(op.mul, xrange(1, r+1), 1)
    return numer//denom

def MakeGraph(df, L):
    # ONLY CORRECT FOR ZERO-ONE CURRENTLY
    G = nx.DiGraph()
    nt = np.sum(df["weight"])/L
    for i in range(len(df)):
        G.add_edge(df.iloc[i, 0],
                   df.iloc[i, 1],
                   weight=float(df.iloc[i, 2]),
                   inv_lweight=-np.log(float(df.iloc[i, 2]/nt)))
    
    # Calculate the probability of taking a particular edge given
    # at a particular node
    for node in G.nodes():
        tweight = 0
        for edge in G.out_edges(node):
            tweight += G.edges[edge]["weight"]
        for edge in G.out_edges(node):
            G.edges[edge]["cprob"] = G.edges[edge]["weight"]/tweight
            G.edges[edge]["inv_lcprob"] = -np.log(G.edges[edge]["weight"]/tweight)
    return G, nt

def AddNodeWeight(G, df, L):
    for i in range(len(df)):
        G.nodes()[df.iloc[i, 0]]["weight"] = df.iloc[i, 1]
    return G

def MinimumArborescence(G):
    B = nx.minimum_spanning_arborescence(G, attr="inv_lcprob")
    for edge in B.edges():
        if edge in G.edges():
            for key in G.edges()[edge]:
                B.edges()[edge][key] = G.edges()[edge][key]
    for node in B.nodes():
        if node in G.nodes():
            for key in G.nodes()[node]:
                B.nodes()[node][key] = G.nodes()[node][key]
    return B

def CheckGraph(G):
    PROBLEM = 0
    for node in G.nodes():
        in_weight=0
        out_weight=0
        node_weight = G.nodes()[node]["weight"]
        for edge in G.out_edges(node):
            out_weight += G.edges()[edge]["weight"]
        for edge in G.in_edges(node):
            in_weight += G.edges()[edge]["weight"]

        if in_weight > node_weight or out_weight > node_weight:
            print "Node", "Node weight", "In weight", "Out weight"
            print node, node_weight, in_weight, out_weight
            print "NODE PROBLEM..:", node
            PROBLEM = 1
    if PROBLEM == 0:
        print "All weights correct..."
    return
    

def ZeroOrderProb(state_i, L, node_type = args.node_type):
    if node_type == "int":
        ones = (BinFromInt(state_i, L)).count(1)
    else:
        ones = StringToState(state_i).count(1)
    return 1/ncr(L,ones) * 1/(L-ones)

def ProbOfPath(path, G):
    prob = 1
    for i in range(len(path)-1):
        prob*=G.edges[(path[i],path[i+1])]["cprob"]
    return prob

def CalculateProbOfPaths(paths, G):
    return np.sum([ProbOfPath(path, G) for path in paths])

def k_shortest_paths(G, source, target, k, weight=None):
    paths = []
    for path in nx.shortest_simple_paths(G, source, target, weight=weight):
        paths.append(path)
        if len(paths) > k:
            return paths
def k_prob_shortest_paths(G, source, target, p=0.95, weight=None):
    paths = []
    cp = 0
    count = 0
    d=[]
    for path in nx.shortest_simple_paths(G, source, target, weight=weight):
        d.append([count, cp])
        paths.append(path)
        cp += ProbOfPath(path, G)
        count += 1
        if cp > p:
            df = pd.DataFrame(d, columns=["count","cdf"])
            return paths, df

def MakeNewGraph(G, paths):
    G2 = nx.DiGraph()
    for path in paths:
        for i in range(len(path)-1):
            from_node = path[i]
            to_node   = path[i+1]
            weight = G.edges[(from_node, to_node)]["weight"]
            G2.add_edge(from_node, to_node, weight=weight)
    return G2        

def Hamming(v1, v2):
    count = 0
    for i, el in enumerate(v1):
        if v1[i] != v2[i]:
            count += 1
    return count

def CalculateAngle(v, G, pos, direction_weight=args.direction_weight, node_type=args.node_type):
    ones = v.count(1)
    state = StateToString(v)
    if ones == 0:
        return 0
    if ones == 1:
        el = v.index(1)
        ac = np.arccos(1-2*(el/(len(v)-1)))/(np.pi) * 180
        return ac
    if node_type == "int":
        K = ConvertToInt(v)
    if node_type == "string":
        K = state
    in_edges = G.in_edges(K)
    weights = 0
    weighted_angle = 0
    
    for edge in in_edges:
        w = G.edges[edge]['weight']
        f = edge[0]
        weights += w**direction_weight
        weighted_angle += pos[f][4] * w**direction_weight
    ew = weighted_angle/weights
    if ew > 181:
        print ew, "Over-expectation"
    return ew

def CoordsSphere(r, theta, phi):
    x = -r * np.cos(theta/180 * np.pi)
    y = r * np.sin(theta/180 * np.pi) * np.cos(phi/180 * np.pi)
    return (x, y, r, theta, phi)

def AddPosition(pos, k, G, L, node_type = args.node_type):
    if node_type == "int":
        state = BinFromInt(k, L)
    if node_type == "string":
        state = StringToState(k)
    theta = np.arccos(1-2*state.count(1)/L) /np.pi * 180
    phi = CalculateAngle(state, G, pos)
    pos[k] = CoordsSphere(1, theta, phi)
    return

def PrintPos(pos):
    for key in pos:
        print pos[key]

def AEQ(a,b,ep=1e-8):
    return True if abs(a-b) < ep else False
    
def GroupPos(pos, L, node_type = args.node_type):
    pos2 = {}
    for key in pos:
        els = int(round(L*(pos[key][0]+1))/2)
        if node_type == "int":
            state = tuple([tuple(BinFromInt(key, L))])
        else:
            state = tuple([tuple(StringToState(key))])
        try:
            pos2[els].append(pos[key] + state + tuple([key]))
        except KeyError:
            pos2[els] = [pos[key] + state + tuple([key])]
    return pos2

def Hamming(s1, s2):
    count = 0
    for i, el in enumerate(s1):
        if s1[i]!=s2[i]:
            count += 1
    return count

def MaxHam(n, L):
    if n <= np.floor(L/2):
        return 2*n
    return 2*(L-n)

def MakeFixedPosition(G, L, node_type = args.node_type):
    fixed_position = {}
    state_dic = {}
    for node in G.nodes():
        if node_type == "int":
            state = BinFromInt(node, L)
        else:
            state = StringToState(node)
        ones = state.count(1)
        state_dic[node] = ones
    state_dic = sorted(state_dic.items(), key=op.itemgetter(1))
    
    for node in state_dic:
        AddPosition(fixed_position, node[0], G, L)
    return fixed_position


def AdjustGroupPos(pos, L, gpos, gamma=0.0, gap_type="uniform"):
    _pos = {}
    for key in gpos:
        el1 = gpos[key]
        el2 = sorted(el1, key=lambda tup: tup[1])
        amin_d = -np.sin(np.arccos((1-2*key/L)))
        amax_d = -amin_d
        span = abs(amax_d-amin_d)

        if len(el2)>1:
            dists = [el2[i][1]-el2[i-1][1] for i in range(1,len(el2))]
            hammings = [Hamming(el2[i+1][-2], el2[i][-2]) for i in range(len(dists))]
            hamming_max = max(hammings)/MaxHam(key, L)
            if AEQ(el2[0][1], el2[-1][1]):
                min_d = amin_d
                max_d = amax_d
            else:
                min_d = el2[0][1]
                max_d = el2[-1][1]

            if gap_type == "uniform":
                gap   = (max_d - min_d)/len(dists)
                dists2 = np.array([(1-gamma)*el + gamma*gap for el in dists])
                dists2 = dists2/np.sum(dists2) * (max_d - min_d)
            else:
                if gap_type == "hamming":
                    dists2 = []
                    hamming_sum = np.sum(hammings)

                    for i, el in enumerate(dists):
                        _el = (1-gamma)*el + gamma*Hamming(el2[i+1][-2],
                                      el2[i][-2]) * (max_d-min_d)/hamming_sum
                        dists2.append(_el)
                    dists2 = np.array(dists2)/np.sum(dists2) * (max_d - min_d)

            cdists = [min_d]
            for i in range(len(dists2)):
                cdists.append(cdists[i]+dists2[i])
            for i, el in enumerate(el2):
                _pos[el[-1]] = tuple([el[0],cdists[i]])
        else:
            _pos[el2[0][-1]] = tuple([el2[0][0],el2[0][1]])
    return _pos

def CalculateNodeSizes(G,
                       nt,
                       node_normalise=args.node_normalise,
                       node_scale_exponent=args.node_scale_exponent,
                       node_scale=args.node_scale,
                       node_type=args.node_type):
    node_sizes = []
    for node in G.nodes():
        node_sizes.append(G.nodes()[node]["weight"])
    m = max(node_sizes)
    return [((el/m)**node_scale_exponent)*node_scale for el in node_sizes]

def Noverlap(G, pos, gpos, node_sizes, space_to_diam=0.5, node_scale_exponent=0.5, aspect=args.aspect):
    _pos = {}
    adjustments = {}
    for key in gpos:
        el1 = gpos[key]
        el2 = sorted(el1, key=lambda tup: tup[1])
        amin_d = el2[0][1]
        amax_d = el2[0][-1]
        amin_d = -np.sin(np.arccos((1-2*key/L)))
        amax_d = -amin_d
        diams = []
        for el in el2:
            node = el[-1]
            diams.append(node_sizes[node]**0.5)
        if len(el2) < 2:
            continue
        adjustments[key] = (amax_d - amin_d)*(1-space_to_diam)/np.sum(diams)
        
    min_adjustment = min(adjustments.values())
    node_sizes_adjust = {}
    for node in node_sizes:
        node_sizes_adjust[node] = (node_sizes[node]**0.5)*min_adjustment/2
    max_m = max([node_sizes_adjust[node] for node in node_sizes])
    max_r = (2/L)/args.max_width
    min_adj = max(1, max_m/max_r)

    for node in node_sizes:
        node_sizes_adjust[node] = node_sizes_adjust[node]/min_adj
    
    for key in gpos:
        el1 = gpos[key]
        el2 = sorted(el1, key=lambda tup: tup[1])
        amin_d = -np.sin(np.arccos((1-2*key/L)))
        amax_d = -amin_d
    for key in gpos:
        el1 = gpos[key]
        el2 = sorted(el1, key=lambda tup: tup[1])
        if len(el2) > 1:
            amin_d = -np.sin(np.arccos((1-2*key/L)))
            amax_d = -amin_d
            dists = []
            diams = [node_sizes_adjust[el[-1]]*2 for el in el2]
            gap = ((amax_d - amin_d) - np.sum(diams))/(len(el2)-1)
            for i in range(1,len(diams)):
                dist1 = diams[i-1]/2
                dist2 = diams[i]/2
                distt = dist1 + dist2 + gap
                dists.append(distt)
            dists = np.array(dists)
            cdists = [amin_d]
            for i, el in enumerate(dists):
                cdists.append(cdists[i]+el)
            for i, el in enumerate(el2):
                _pos[el[-1]] = tuple([el[0],cdists[i]])
        else:
            _pos[el2[0][-1]] = tuple([el2[0][0],el2[0][1]])
    return _pos, node_sizes_adjust




def Adjust(w, factor, mult):
    w2 = []
    m = max(w)
    for el in w:
        w2.append(((el/m)**factor)*mult)
    return w2


def MakeColors(weights, scale=args.edge_scale, colormap="Greys"):
    import matplotlib as mpl
    import matplotlib.cm as cm
    minima = 0
    maxima = 1
    norm = mpl.colors.Normalize(vmin=minima, vmax=maxima)
    mapper = cm.ScalarMappable(norm=norm, cmap=plt.get_cmap(colormap))
    colors = []
    scaler = scale if scale > 0 else 1
    for v in weights:
        colors.append(mapper.to_rgba(v/scaler))
    return colors

def MakeTopSizesByOnes(G, node_sizes, node_type = args.node_type):
    import operator
    top_sizes = {}
    for i, node in enumerate(G.nodes()):
        if node_type == "int":
            state = BinFromInt(node, L)
        else:
            state = StringToState(node)
        ones = state.count(1)
        ns = node_sizes[i]
        if ones not in top_sizes:
            top_sizes[ones] = {}
        top_sizes[ones][node] = ns
    for key in top_sizes:
        top_sizes[key] = sorted(top_sizes[key].items(), key=operator.itemgetter(1), reverse=True)
    return top_sizes

def DrawNodesAsCircles(G, pos, ax, radii, aspect = args.aspect, **kwargs):
    Nodes = []
    for i, node in enumerate(G.nodes()):
        circle = Ellipse(pos[node], width = 2*radii[i]*aspect, height=2*radii[i], **kwargs)
        circle.set_linewidth(0)
        Nodes.append(circle)
        ax.add_patch(circle)
    return Nodes

def SortedOutEdges(node, G):
    edgews = {}
    for edge in G.out_edges(node):
        edgews[edge] = G.edges[edge]['weight']
    edgews = sorted(edgews.items(), key=op.itemgetter(1), reverse=True)
    return edgews
    
def Diff(node1, node2, L, node_type = args.node_type):
    if args.node_type == "int":
        b1 = BinFromInt(node1, L)
        b2 = BinFromInt(node2, L)
    else:
        b1 = StringToState(node1)
        b2 = StringToState(node2)
    for i,el in enumerate(b1):
        if b2[i] == 1 and b1[i] ==0:
            return i

def MakeGreedyLabelSet(state_labels,
                       colormap,
                       labels_definition = args.extra_labels_routes,
                       label_colour_type = args.label_colour_type,
                       node_type = args.node_type,
                       probable_paths = args.probable_paths):
    lr = [int(el) for el in args.extra_labels_routes.split(",")]
    limit = lr[0]
    n_r = lr[1]
    labels = {}
    big_labels = []
    for r in range(n_r):
        if args.node_type == "int":
            start = 0
        else:
            start = StateToString([0]*L)
        current_node = start
        small_labels = OrderedDict()
        for i in range(limit):
            edges = SortedOutEdges(current_node, G)
            updated = 0
            for j in range(len(edges)):
                node_from = edges[j][0][0]
                node_to = edges[j][0][1]
                if probable_paths != None:
                    if node_to in G_2.nodes():
                        if node_to not in labels:
                            el = Diff(node_from, node_to, L)
                            labels[node_to] = state_labels[el]
                            small_labels[node_to] = state_labels[el]
                            updated = 1
                            current_node = node_to
                            break
                else:
                    if node_to not in labels:
                        el = Diff(node_from, node_to, L)
                        # ** adjust with plus for now
                        labels[node_to] = state_labels[el]
                        small_labels[node_to] = str("+") + state_labels[el]
                        updated = 1
                        current_node = node_to
                        break
            if updated == 0:
                big_labels.append(small_labels)
                break
            if i == limit-1:
                big_labels.append(small_labels)
    if label_colour_type == "standard":
        colours = [pal[CycleColors(i)] for i in range(n_r)]
    else:
        colours = MakeColors(weights=range(n_r), scale=(n_r-1)*1.5, colormap=colormap)
    return big_labels, colours

def MakeProbableLabelSet(G,
                         state_labels,
                         colormap,
                         labels_definition = args.extra_labels_routes,
                         label_colour_type = args.label_colour_type,
                         node_type = args.node_type,
                         probable_paths = args.probable_paths):
    lr = [int(el) for el in args.extra_labels_routes.split(",")]
    limit = lr[0]
    n_r = lr[1]
    # Find n_r shortest paths
    if node_type == "int":
        source = 0
        target = 2**L -1 
    if node_type == "string":
        source = StateToString([0]*L)
        target = StateToString([1]*L)
    paths = k_shortest_paths(G=G,
                             source=source,
                             target=target,
                             k=n_r-1,
                             weight="inv_lcprob")
    labels = {}
    big_labels = []
    for path in paths:
        small_labels = OrderedDict()
        for i in range(len(path)-1):
            node_from = path[i]
            node_to = path[i+1]
            if probable_paths != None:
                if node_to in G_2.nodes():
                    if node_to not in labels:
                        if node_type == "int":
                            el = Diff(node_from,
                                      node_to,
                                      L)
                        else:
                            el = Diff(StringToState(node_from),
                                      StringToState(node_to),
                                      L)
                        labels[node_to] = state_labels[el]
                        small_labels[node_to] = state_labels[el]
            else:
                if node_to not in labels:
                    el = Diff(node_from, node_to, L)
                    labels[node_to] = state_labels[el]
                    small_labels[node_to] = state_labels[el]
        big_labels.append(small_labels)
    
    if label_colour_type == "standard":
        colours = [pal[CycleColors(i)] for i in range(n_r)]
    else:
        colours = MakeColors(weights=range(n_r), scale=(n_r-1)*1.5, colormap=colormap)
    return big_labels, colours, paths

def MakeProbableLabelSetData(G,
                             state_labels,
                             colormap,
                             L,
                             transition_data = args.transition_data,
                             labels_definition = args.extra_labels_routes,
                             label_colour_type = args.label_colour_type,
                             node_type = args.node_type,
                             probable_paths = args.probable_paths):
    df = pd.read_csv(transition_data, sep = " ", index_col = None, header=None)
    paths = []
    for i in range(int(len(df)/2)):
        source = list(df.iloc[2*i,:])
        target = list(df.iloc[2*i+1,:])
        if node_type == "int":
            source = ConvertToInt(source)
            target = ConvertToInt(target)
        if node_type == "string":
            source = StateToString(source)
            target = StateToString(target)
        if source not in G.nodes() or target not in G.nodes():
            continue
        try:
            path = k_shortest_paths(G=G,
                                    source=source,
                                    target=target,
                                    k=0,
                                    weight="inv_lcprob")
        except:
            continue
        paths.append(path[0])

    labels = {}
    big_labels = []
    for path in paths:
        small_labels =OrderedDict()
        for i in range(len(path)-1):
            node_from = path[i]
            node_to = path[i+1]
            if probable_paths != None:
                if node_to in G_2.nodes():
                    if node_to not in labels:
                        if node_type == "int":
                            el = Diff(node_from,
                                      node_to,
                                      L)
                        else:
                            el = Diff(StringToState(node_from),
                                      StringToState(node_to),
                                      L)
                        labels[node_to] = state_labels[el]
                        small_labels[node_to] = state_labels[el]
            else:
                if node_to not in labels:
                    el = Diff(node_from, node_to, L)
                    labels[node_to] = state_labels[el]
                    small_labels[node_to] = state_labels[el]
        big_labels.append(small_labels)
    
    if label_colour_type == "standard":
        colours = [pal[CycleColors(i)] for i in range(len(paths))]
    else:
        colours = MakeColors(weights=range(len(paths)), scale=len(paths)-1, colormap=colormap)
    return big_labels, colours, paths


# MAIN BLOCK
font = {'family':'sans-serif',
        'size':args.fontsize,
        'sans-serif':['Arial']}
matplotlib.rc('font', **font)

prefix = ""
file_name = args.f

datafile = prefix + file_name

if file_name.find("match-data") != -1 or file_name.find("zero-one") != -1:
    node_tag = "-edge-list-long"
    node_tag2= "-state-list-long"
    args.write_transitions = "no"
else:
    node_tag = "-edge-list-int" if args.node_type == "int" else "-edge-list-string"
outfile = datafile.replace("list-pord", "list-pord"+node_tag).replace(".csv", ".txt")
outfile2= datafile.replace("list-pord", "list-pord"+node_tag2).replace(".csv", ".txt")

outfile_graph = os.path.dirname(args.f)
if len(outfile_graph) > 0:
    outfile_graph += "/"
outfile_graph = outfile if args.outfile_graph == None else outfile_graph + args.outfile_graph


df = pd.read_csv(datafile, index_col=None)
L = max(df["feature"]+1)
if args.radial_label_distance == -1:
    args.radial_label_distance = 1/L

args.extra_labels_routes = args.extra_labels_routes if args.extra_labels_routes != None else "1," + str(L-1) + ",0"
if args.write_transitions != "no":
    routes = MakeRoutes(df)
    uts = np.random.shuffle(range(min(len(routes), int(args.tlimit))))
    transitions = MakeTransitions(routes, uts=uts, node_type = args.node_type)
    WriteOutTranstions(transitions, outfile)

df = pd.read_csv(outfile,  index_col = None, sep=" ")
df2= pd.read_csv(outfile2, index_col = None, sep=" ")

G, nt = MakeGraph(df=df, L=L)
G = AddNodeWeight(G, df=df2, L=L)
if args.arb != None:
    B = MinimumArborescence(G)
    G = B
CheckGraph(G)    




if args.probable_paths != None:
    if args.node_type == "int":
        source = 0
        target = 2**L -1 
    if args.node_type == "string":
        source = StateToString([0]*L)
        target = StateToString([1]*L)
    paths, df_d = k_prob_shortest_paths(G=G,
                                        source=source,
                                        target=target,
                                        p=args.probable_paths,
                                        weight="inv_lcprob")
    G_2 = MakeNewGraph(G, paths)
    print G_2.edges()
    fig, ax = plt.subplots(1,1)
    if L < 9:
        n_p = len(list(nx.all_simple_paths(G,
                                           source=source,
                                           target=target)))
        ax.plot(df_d.iloc[:,0]/n_p, df_d.loc[:,"cdf"])
        ax.set_xlabel("Proportion of paths found")
    else:
        ax.plot([el for el in df_d.loc[:,"count"]], df_d.loc[:,"cdf"])
        ax.set_xlabel("Paths found")
    ax.set_ylim(0,1)
    ax.set_xscale('log')
    ax.set_ylabel("Probability of first x paths")
    out = outfile_graph.replace("edge-list"+node_tag,"count-cdf")
    out = out + ".png" if out.find(".txt") == -1 else out.replace(".txt", ".png")
    plt.savefig(out, dpi=300)



for node in G.nodes():
    if len(G.in_edges(node)) < 1:
        print node

fixed_position = MakeFixedPosition(G, L)

if args.layout_type == "noverlap":
    node_sizes = CalculateNodeSizes(G, node_scale_exponent=args.node_scale_exponent, node_scale=args.node_scale, nt=nt)
    node_sizes_d = dict(zip(G.nodes(), node_sizes))
    gps = GroupPos(fixed_position, L=L)
    fixed_positions, node_sizes_d = Noverlap(G=G, pos=fixed_position, gpos=gps,
                                             node_sizes=node_sizes_d,
                                             node_scale_exponent=0.5,
                                             space_to_diam=args.noverlap_space)
    node_sizes = [node_sizes_d[node] for node in G.nodes()]

if args.layout_type == "spring":
    gps = GroupPos(fixed_position, L)
    fixed_positions = AdjustGroupPos(fixed_position, L, gps, gamma=args.gamma, gap_type=args.gap_type)
    node_sizes = CalculateNodeSizes(G=G,
                                    node_normalise=args.node_normalise,
                                    node_scale_exponent=args.node_scale_exponent*0.5,
                                    node_scale=args.node_scale*(2/L)*(1/args.max_width), nt=nt)

    node_sizes_d = dict(zip(G.nodes(), node_sizes))

if args.layout_type == "cube":
    a = 1
    b = 1/3
    fixed_positions = {0:(-1,0), 1:(-b,2*b), 2:(-b,0), 4:(-b,-2*b), 3:(b,2*b), 5:(b,0), 6:(b,-2*b), 7:(1,0)}
    node_sizes = CalculateNodeSizes(G=G,
                                    node_normalise=args.node_normalise,
                                    node_scale_exponent=args.node_scale_exponent*0.5,
                                    node_scale=args.node_scale*(2/L)*(1/args.max_width), nt=nt)

    node_sizes_d = dict(zip(G.nodes(), node_sizes))




fixed_nodes = fixed_positions.keys()
nodes = G.nodes()
edges = G.edges()
if args.edge_normalise == 1.0:
    weights = [G[u][v]['weight'] for u,v in edges]
    if args.probable_paths != None:
        weights_2 = [G_2[u][v]['weight'] for u,v in G_2.edges()]
else:
    weights = []
    for u,v in edges:
        if args.node_type == "int":
            ones = BinFromInt(u, L).count(1)
        else:
            ones = StateToString(u).count(1)
        weights.append((G[u][v]['weight']/nt) / (1/ncr(L, ones) * 1/(L-ones)) * nt)
    if args.probable_paths != None:
            weights_2 = []
            for u,v in G_2.edges():
                if args.node_type == "int":
                    ones = BinFromInt(u, L).count(1)
                else:
                    ones = StateToString(u).count(1)
                weights_2.append((G_2[u][v]['weight']/nt) / (1/ncr(L, ones) * 1/(L-ones)) * nt)


weights_original = Adjust(weights, float(args.edge_alpha), float(args.edge_scale))
weights = Adjust(weights, float(args.edge_scale_exponent), float(args.edge_scale))
if args.probable_paths != None:
    weights_original_2 = Adjust(weights_2, float(args.edge_alpha), float(args.edge_scale))
    weights_2 = Adjust(weights_2, float(args.edge_scale_exponent), float(args.edge_scale))    



# Plotting
old_colorblind = ["#0072b2", "#009e73", "#d55e00", "#cc79a7", "#f0e442", "#56e4b9"]
plt.clf()
plt.style.use('seaborn-colorblind')
ncolors_old = 6
ncolors = ncolors_old
pal = sns.color_palette("colorblind", ncolors_old)
pal = sns.color_palette(old_colorblind)

fig = plt.figure()
fig.set_size_inches(args.width,args.width*args.aspect)
ax = fig.add_subplot(111)
pos = fixed_positions
colors = MakeColors(weights)
if args.node_normalise == 1.0 or args.edge_normalise == 1.0:
    colorn = pal[0]
else:
    colorn = pal[0]

def CycleColors(i, ncolors=ncolors, taken=2):
    return (i%(ncolors-taken)+int(taken/2))
    
def CycleColors(i, ncolors=ncolors, taken=[0,2]):
    els = list(set(range(ncolors)) - set(taken))
    return els[i]

def CycleColors(i, ncolors=ncolors, taken=[0,2]):
    return i
    
mw = max(node_sizes)/((2/L)*(1/args.max_width))*8/L
if args.probable_paths != None:
    node_sizes_2 = [node_sizes_d[el] for el in G_2.nodes()]
    
        
if args.probable_paths == None:
    widths = [el*mw for el in weights]
    Edges = nx.draw_networkx_edges(G, pos, edges=edges, width=widths,
                                   edge_color = args.edge_color,
                                   arrowstyle = mpl.patches.ArrowStyle("-"),
                                   ax = ax,
                                   node_size=0)

    Nodes = DrawNodesAsCircles(G=G, pos=pos, radii=node_sizes, aspect=args.aspect, alpha=args.node_alpha,
                                    color=colorn,
                                    ax=ax)
    if args.show_edge_alpha == "yes":
        for i, edge in enumerate(Edges):
            edge.set_alpha(weights_original[i]/args.edge_scale)
        
else:
    widths_2 = [el*mw for el in weights_2]
    Edges = nx.draw_networkx_edges(G_2, pos, edges=G_2.edges(), width=widths_2,
                                   edge_color = args.edge_color,
                                   ax = ax,
                                   arrowstyle = mpl.patches.ArrowStyle("-"), node_size=0)

    Nodes = DrawNodesAsCircles(G=G_2, pos=pos, radii=node_sizes_2,
                               aspect=args.aspect,
                               alpha=args.node_alpha,
                               color=colorn, ax=ax)
    if args.show_edge_alpha == "yes":
        for i, edge in enumerate(Edges):
            edge.set_alpha(weights_original_2[i]/args.edge_scale)


if args.show_se != "no":
    labels_ends = {}
    if args.node_type == "int":
        start = 0
        end = 2**L-1
    else:
        start = StateToString([0]*L)
        end = StateToString([1]*L)
    if args.show_se == "bracket":
        labels_ends[start] = "{}"
    else:
        labels_ends[start] = "no\nfeatures"
        labels_ends[end] = "all\nfeatures"
    if args.show_se != "bracket":
        nx.draw_networkx_labels(G, pos, labels_ends, font_size=args.labels_fontsize, font_color='black')
    else:
        if args.show_se.lower() != "none":
            nx.draw_networkx_labels(G, pos, labels_ends, font_size=args.labels_fontsize, font_color='white')

labels = {}
top_sizes = MakeTopSizesByOnes(G, node_sizes)
if args.extra_labels != None:
    if args.labels != None:
        state_labels = [str(el) for el in list(pd.read_csv(args.labels, header=None, index_col = None).iloc[:,1])]
    else:
        state_labels = range(L)
    ns = [int(el) for el in (args.extra_labels).split(",")]
    lim = ns[-1]
    ns = ns[:-1]
    for i, n in enumerate(ns):
        count = 0
        for j, el in enumerate(top_sizes[n]):
            if j >= lim:
                break
            node = el[0]
            if args.node_type == "int":
                state = BinFromInt(node, L)
            else:
                state = StringToState(node)
            label = []
            for k, el in enumerate(state):
                if el == 1:
                    label.append(str(state_labels[k]))
            label =  "{" + (",".join(label)) + "}"
            labels[node] = label

if args.extra_labels_routes != None:
    if int(args.extra_labels_routes.split(",")[1]) == 0:
        args.extra_labels_routes = None

print args.extra_labels_routes

if args.extra_labels_routes != None:
    if args.labels != None:
        state_labels = [str(el) for el in list(pd.read_csv(args.labels, header=None, index_col = None).iloc[:,1])]
    else:
        state_labels = [str(i+1) for i in range(L)]
    if args.label_type == "greedy":
        big_labels, label_colours = MakeGreedyLabelSet(state_labels=state_labels,
                                                       colormap=args.colormap)
    else:
        if args.label_type == "probable":
            big_labels, label_colours, paths = MakeProbableLabelSet(G=G,
                                                                    state_labels=state_labels,
                                                                    colormap=args.colormap)
        else:
            if (args.label_type == "probable_data" or args.label_type == "greedy_data" or args.label_type == "None") and args.transition_data != None:
                big_labels, label_colours, paths = MakeProbableLabelSetData(G=G,
                                                                            state_labels=state_labels,
                                                                            colormap=args.colormap,
                                                                            transition_data=args.transition_data, L=L)
            if args.label_type == "greedy_data":
                big_labels, label_colours = MakeGreedyLabelSet(state_labels=state_labels,
                                                               colormap=args.colormap)

            if (args.label_type == "probable_data" or args.label_type == "greedy_data" or args.label_type == "None") and args.transition_data != None:
                G_3=nx.DiGraph()
                for path in paths:
                    el1 = path[0]
                    el2 = path[-1]
                    if el1 in G.nodes():
                        G_3.add_node(el1)
                    if el2 in G.nodes():
                        G_3.add_node(el2)
                node_sizes_3 = [node_sizes_d[el] for el in G_3.nodes()]
                print "Data N:", len(G_3.nodes())
                Nodes_3 = DrawNodesAsCircles(G=G_3, pos=pos, radii=node_sizes_3,
                                             aspect=args.aspect, alpha=args.node_alpha,
                                             color=pal[2],
                                             ax=ax)
                
                
    if args.label_type != "None":
        if args.extra_labels_routes != None:
            keep_els = [int(el) for el in args.extra_labels_routes.split(",")]
            keep_els = keep_els[2:] if len(keep_els) > 2 else []

        for i, el in enumerate(big_labels):
            bbox2 = dict(facecolor='lightgrey', edgecolor='none', boxstyle='round')
            if i not in keep_els and len(keep_els)>1:
                continue
            for node_i, node in enumerate(el):
                lab = el[node]
                pos_lab = pos[node]
                theta=np.arccos(-pos_lab[0])
                lab_abs = args.radial_label_distance
                if args.alternate_label == "no":
                    ysign = 1 if pos_lab[1] > 0 else -1
                    xsign = 1
                else:
                    if args.alternate_label == "yes":
                        ysign = 1 if node_i % 2 == (0 if pos_lab[1] > 0 else 1)  else -1
                        xsign = 1 if node_i % 2 == 0 else -1
                    else:
                        if args.alternate_label == "1":
                            ysign = 1 if pos_lab[1] > 0 and i != 1 else -1
                            xsign = 1

                if args.extra == "yes":
                    extra = (1-np.abs(pos_lab[0]))**4
                else:
                    extra = 0
                text_lab = (pos_lab[0] - xsign*np.cos(theta)*lab_abs,
                            pos_lab[1] + ysign*(np.sin(theta)*lab_abs*(1+extra)))

                if args.label_bbox == "yes":
                    ax.annotate(s=lab,
                                xy=pos_lab,
                                xycoords='data',
                                xytext=text_lab,
                                textcoords='data',
                                color=label_colours[i],
                                fontsize=args.labels_fontsize,
                                bbox=dict(boxstyle="round", facecolor="#efefef", alpha=args.bbox_alpha, edgecolor='none', zorder=-100),
                                horizontalalignment='center',
                                verticalalignment='center')

                    ax.annotate(s="",
                                xy=pos_lab,
                                xycoords='data',
                                xytext=text_lab,
                                textcoords='data',
                                color=label_colours[i],
                                fontsize=args.labels_fontsize,
                                bbox=dict(boxstyle="round", facecolor="none", alpha=0.5, edgecolor='none', zorder=-100),
                                arrowprops=dict(arrowstyle="-", color=label_colours[i], lw=0.2, connectionstyle="arc3,rad=0.",
                                                shrinkA=0, shrinkB=0, alpha=0.5),
                                horizontalalignment='center',
                                verticalalignment='center',
                                zorder=-200)
                else:
                    if args.label_bbox == "no":
                        ax.annotate(s=lab,
                                    xy=pos_lab,
                                    xycoords='data',
                                    xytext=pos_lab,
                                    textcoords='data',
                                    color=label_colours[i],
                                    fontsize=args.labels_fontsize,
                                    horizontalalignment='center',
                                    verticalalignment='center')


                
if args.grid == "no":
    if args.xaxis.lower() == "yes":
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.axes.get_xaxis().set_visible(True)
        xticks = [-1 + 2*i/(L) for i in range(0, L+1, args.xevery)]
        xticklabels = range(0, L+1, args.xevery)
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels, fontsize=args.fontsize)
        ax.xaxis.set_ticks_position('bottom')
        ax.get_yaxis().set_ticks([])
        ax.get_yaxis().set_visible(False)
        ax.set_xlabel(args.xlabel, fontsize=args.fontsize)
    else:
        plt.grid('off')

axlim = args.axlim
axlim_y = axlim/args.aspect
ax.set_xlim((-axlim,axlim))
ax.set_ylim((-axlim_y,axlim_y))
if args.layout_type != "cube":
    ax.set_aspect(aspect = args.aspect, adjustable='datalim')


def MakeLabel(X, app):
    if X.find(".txt") != -1:
        return X.replace(".txt", app)
    else:
        return X + app

if args.out_type == "pdf":
    pdf = PdfPages(MakeLabel(outfile_graph, ".pdf"))
    if args.bbox == "tight":
        pdf.savefig(bbox_inches = 'tight',figure=fig)
    else:
        pdf.savefig(figure=fig)
    pdf.close()
if args.out_type == "png":
    if args.bbox == "tight":
        fig.savefig(MakeLabel(outfile_graph, ".png"), dpi=600, bbox_inches = "tight")
    else:
        fig.savefig(MakeLabel(outfile_graph, ".png"), dpi=600)
if args.svg.lower() == "yes":
    fig.savefig(MakeLabel(outfile_graph, ".svg"))
