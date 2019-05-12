#!/usr/bin/env python

# _____________________________________________________________________________
# Description:
#     This python file takes an input sample dataframe (and phylogeny) and
#     outputs the transition dataset that can be parsed by HyperTraPS for
#     performing inference

# Python Package requirements (all available on conda):
#     ete3
#     numpy
#     pandas
#     networkx
#     colleactions
#     argparse
#     Optional:
#        pygraphviz

# Inputs:
#     -data:
#         A comma seperated dataframe file of form:
        
#         Sample ID    , feature1 , feature2, ... ,feature L
#         1            ,       1_1,      1_2, ... , 1_L
#         ...
#         N            ,       N_1,      N_2, ... , N_L
    
#         This can be cross-sectional, longitudinal (must be ordered by collection time)
#         or phylogenetic data (must include a tree)

#     -input_type:
#         "cross-sectional", "longitudinal", "sample-tree" or "phylogenetic"

#     -phylogeny:
#         A phylogeny where leaves match the "Sample Index" passed to -data
#         Must be in Newick format and is imported with ete3
#         If unrooted, sets root as the outgroup of the midpoint (arbitrary)

#     -phylogeny_format:
#         The format type for the ete3 import:
#             see:
#                 http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#reading-and-writing-newick-trees
#     -constant_transitions:
#         0 or 1
#         Whether or not to include transitions in dataset where no change occurs.
#         Adds a constant amount to liklihood so is set to 0 to switch off as default

#     -labels:
#         The header with the feature names are output to the file passed here
#         Default: labels.csv

#     -outfile:
#         The transition dataset that is passed to HyperTraPS is passed here.
#         Default: transitions.txt

#     -outfile_cs:
#         The cross-section (cross-sectional and leaves for phylogeny) or end for longitudinal
        
# Outputs:
#     - Transition dataset for input to HyperTraPS
#     - Raw data of cross-section
#     - Labels file for use in post-HyperTraPS analysis procedures
    
# _____________________________________________________________________________

from __future__ import division
from ete3 import Tree
import numpy as np
import pandas as pd
import networkx as nx
from collections import OrderedDict
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-data", type=str, required=False)
parser.add_argument("-input_type", type=str, required=False, default="phylogenetic")
parser.add_argument("-phylogeny", type=str, required=False, default=None)
parser.add_argument("-phylogeny_format", type=int, required=False, default=3)
parser.add_argument("-sample_tree", type=str, required=False, default=None)
parser.add_argument("-constant_transitions", type=int, required=False, default=1)
parser.add_argument("-labels", type=str, required=False, default="labels.csv")
parser.add_argument("-graphviz", type=str, required=False, default="no")
parser.add_argument("-outfile", type=str, required=False, default="transitions.txt")
parser.add_argument("-outfile_cs", type=str, required=False, default="raw-data.txt")
args, unknown = parser.parse_known_args()

print args

def LoadData(file_name, keep_nan=0, index_col=0):
    df = pd.read_csv(file_name, sep=",", index_col=index_col).astype(int)
    df = df.replace(r'\s+', np.nan, regex=True)
    if keep_nan == 1:
        df = df.replace(np.nan, 2)
    else:
        df = df.dropna()
    return df

def LoadStringsToDataDict(df):
    index_labels = OrderedDict()
    for index, row in df.iterrows():
        index_labels[str(index)] = list(row)
    return index_labels

def WriteTransitions(obs, outfile):
    obs_flat = []
    for ob in obs:
        obs_flat.append(ob[0])
        obs_flat.append(ob[1])
    df = pd.DataFrame(obs_flat)
    df.to_csv(outfile, sep=" ", header=None, index=None)

def WriteCS(df_in, outfile):
    obs_flat = []
    for i in range(len(df_in)):
        obs_flat.append(df_in.iloc[i,:])
    df = pd.DataFrame(obs_flat)
    df.to_csv(outfile, sep=" ", header=None, index=None)
    
def WriteLabels(df, simpler=[], outfile="labels.csv", start=0):
    labels = []
    for i, el in enumerate(df.columns[start:]):
        if len(simpler) != 0:
            labels.append([el, simpler[i]])
        else:
            labels.append([str(i), str(el)])
    df_labels = pd.DataFrame(labels)
    df_labels.to_csv(outfile, header=None, index=None)

    
def MakeCSObservations(df):
    obs = []
    L = len(df.columns)
    G = nx.DiGraph()
    for index, row in df.iterrows():
        obs.append([[0]*L, list(row)])
        n_c = len(G.nodes())
        G.add_node(n_c, state=tuple([0]*L))
        G.add_node(n_c+1, state=tuple(list(row)))
        G.add_edge(n_c, n_c+1)
    return obs, G

def MakeLongitudinalObservations(df, constant_transitions=0, verbose=0):
    obs = []
    cs_obs = []
    start = 1
    L = len(df.iloc[0,start:])
    ids = list(set(df.iloc[:,0]))
    if verbose == 1 or True:
        print ids
    G = nx.DiGraph()
    for idx in ids:
        df_t = df[df.iloc[:,0]==idx]
        count = 0
        temp_obs = []
        for index, row in df_t.iterrows():
            row_l = list(row[start:])
            if count == 0 and row_l.count(0) != L:
                temp_obs.append(tuple([0]*L))
            temp_obs.append(tuple(row_l))
            count += 1
        def SetOrder(seq, constant_transitions):
            seen = set()
            seen_add = seen.add
            if constant_transitions == 1:
                return [x for x in seq]
            return [x for x in seq if not (x in seen or seen_add(x))]

        temp_obs = SetOrder(temp_obs, constant_transitions=constant_transitions)
        print temp_obs
        for i, el in enumerate(temp_obs):
            if i == len(temp_obs)-1:
                break
            obs.append([list(temp_obs[i]),list(temp_obs[i+1])])
            n_c = len(G.nodes())
            if i == 0:
                G.add_node(n_c, state=tuple(temp_obs[i]))
                G.add_node(n_c+1, state=tuple(temp_obs[i+1]))
                G.add_edge(n_c, n_c+1)
            else:
                G.add_node(n_c, state=tuple(temp_obs[i+1]))
                G.add_edge(n_c-1, n_c)
        cs_obs.append(list(temp_obs[-1]))
    return obs, pd.DataFrame(cs_obs), G


def BitAnd(a, b):
    _p = []
    for i, el in enumerate(a):
        if a[i] == b[i]:
            pk = a[i]
        else:
            if (a[i] == 2 and b[i] != 0) or (a[i] != 0 and b[i] == 2):
                pk = 2
            else:
                pk = 0
        _p.append(pk)
    return _p

def BitEqual(a, b):
    for i, el in enumerate(a):
        if a[i] != b[i]:
            return False
    return True

def MakeInitialLeaves(tree, node_dict, verbose=0):
    leaves = []
    children = {}
    parent_counts = {}
    for leaf in tree:
        if leaf.up in children:
            # If parent already has recorded children, store tuple:
            # (childA, childB, parent)
            leaves.append((leaf, children[leaf.up], leaf.up))
            parent_counts[leaf.up] +=1
        else:
            children[leaf.up] = leaf
            parent_counts[leaf.up] = 1

    # Check whether there are any parents that have more than 2 leaves
    for key in parent_counts:
        if parent_counts[key] > 2:
            print key, parent_counts[key]
            raise ValueError('Parent has more than two children...')
    if verbose == 1:
        for leaf in leaves:
            print leaf, parent_counts[leaf[2]]
    return leaves

def MakeObservations(tree, node_dict, label_dict, constant_transitions = 0, verbose = 0, check_root = 0, resolve_single_leaves = 0):
    # Algorithm
    # 1. Pick two leaves at random that share a branch point: a and b
    # 2. Calculate a AND b = x
    # 3. if x /= a, then add x -> a to observations
    # 4. if x /= b, then add x -> b to observations
    # 5. replace a and b with x at their root
    # 6. go to 1 while there remain leaves
    G = nx.DiGraph()
    leaves = MakeInitialLeaves(tree, node_dict)
    obs = []
    while len(leaves) != 0:
        if verbose == 1:
            print leaves
        leaf_tuple = leaves[0]
        
        # Work out parent state with BIT AND 
        a_s = node_dict[leaf_tuple[0]]
        b_s = node_dict[leaf_tuple[1]]
        p_s = BitAnd(a_s, b_s)

        if verbose == 1:
            print a_s, type(a_s)
            print b_s, type(b_s)
            print p_s, type(p_s)
        
        # Add parent node to node dict and label dict
        node_dict[leaf_tuple[2]] = p_s
        label_dict[leaf_tuple[2]] = len(label_dict)
        
        G.add_node(label_dict[leaf_tuple[2]], state=tuple(p_s), n=p_s.count(1))
        G.add_node(label_dict[leaf_tuple[0]], state=tuple(a_s), n=a_s.count(1))
        G.add_node(label_dict[leaf_tuple[1]], state=tuple(b_s), n=b_s.count(1))
        G.add_edge(label_dict[leaf_tuple[2]], label_dict[leaf_tuple[0]])
        G.add_edge(label_dict[leaf_tuple[2]], label_dict[leaf_tuple[1]])
        
        # Add observations
        if constant_transitions == 0:
            if BitEqual(a_s, p_s) == False:
                obs.append([p_s, a_s])
            if BitEqual(b_s, p_s) == False:
                obs.append([p_s, b_s])
        if constant_transitions == 1:
            obs.append([p_s, a_s])
            obs.append([p_s, b_s])

        if verbose == 1:
            print leaf_tuple

        # Delete children nodes
        parent = leaf_tuple[2]
        parent.remove_child(leaf_tuple[0])
        parent.remove_child(leaf_tuple[1])
        leaves.remove(leaf_tuple)
        
        for leaf in tree:
            # Check whether p_s shares any parents
            if leaf.up == parent.up and leaf != parent:
                leaves.append((leaf, parent, parent.up))
                break

        if resolve_single_leaves == 1:
            tree, node_dict, label_dict = CollapseSingleLeaves(tree, node_dict, label_dict)
                
        if resolve_single_leaves == 1 and len(leaves) == 0:
            for leaf in tree:
                if leaf.is_root():
                    continue
                parent = leaf.up
                if len(parent.children) == 2:
                    childA = parent.children[0]
                    childB = parent.children[1]
                    if childA.is_leaf() and childB.is_leaf():
                        leaves.append((parent.children[0], parent.children[1], parent))
                        break
                    
        if parent.is_root() == True:
            try:
                p_s = node_dict[parent]
            except:
                if verbose == 1:
                    print "Can't print parent state currently", parent
            if verbose == 1:
                print len(parent.children)
            if p_s.count(1) != 0 or constant_transitions == 1:
                obs.append([[0]*len(p_s), p_s])
    if verbose == 1:
        for leaf in tree:
            print leaf
            try:
                print node_dict[leaf]
            except:
                continue
    return obs, G

    
def CheckLeaves(tree):
    for leaf in tree:
        parent = leaf.up
        if len(parent.children)!= 2:
            print parent, parent.children


def MakeNodeDict(tree, leaf_dict, verbose = 0, resolve_single_leaves=0):
    node_dict = {}
    label_dict = {}
    for key in leaf_dict:
        try:
            key_p = tree.get_leaves_by_name(name=key)[0]
            node_dict[key_p] = leaf_dict[key]
            label_dict[key_p] = len(label_dict)
        except IndexError:
            if verbose == 1:
                print "INDEX ERROR: ", key
            continue
            
    if resolve_single_leaves == 1:
        tree, node_dict, label_dict = CollapseSingleLeaves(tree, node_dict, label_dict)
    
    if verbose == 1:
        print "Number of leaves in tree:", len(tree)
        print "Number of keys in node_dict:", len(node_dict)
    return node_dict, label_dict

def CollapseSingleLeaves(tree, node_dict, label_dict, verbose=0):
    single_leaves = True
    while single_leaves:
        single_leaves = False
        for leaf in tree:
            parent = leaf.up
            if parent != None:
                if len(parent.children) == 1:
                    try:
                        node_dict[parent] = node_dict[leaf]
                        del node_dict[leaf]
                    except:
                        print "Leaf not in node dict..."
                    try:
                        label_dict[parent] = label_dict[leaf]
                        del label_dict[leaf]
                    except:
                        print "Lead not in label dict"
                    if verbose == 1 and False:
                        print node_dict[parent]
                    leaf.detach()
                    single_leaves = True
                    break
    return tree, node_dict, label_dict

def MatchLabelsToLeaves(tree, df, resolve_single_leaves=0, verbose=0):
    leaf_dict = LoadStringsToDataDict(df)
    REMOVE = True
    while REMOVE:
        REMOVE = False
        for leaf in tree:
            if leaf.name not in leaf_dict:
                leaf.detach()
                REMOVE = True
                break
    if verbose == 1:
        print "Number of leaves in tree:", len(tree)
        print "Number of labels in leaf_dict:", len(leaf_dict)
        
    return tree, leaf_dict
    
def LoadTree(phylogeny, format_type, root_if_unrooted):
    t = Tree(args.phylogeny, format=format_type)
    # Check if tree is rooted and root to midpoint if not
    if len((t.get_tree_root()).children) > 2 and root_if_unrooted == 1:
        R = t.get_midpoint_outgroup()
        t.set_outgroup(R)
    return t

def PrintGraphProperties(G, verbose=0):
    if verbose == 1:
        for edge in G.edges():
            print edge, G.nodes[edge[0]], G.nodes[edge[1]]
    H = G.to_undirected()
    print "Connected components: ", nx.number_connected_components(H)
    print "Number of nodes:      ", len(G.nodes())
    print "Number of edges:      ", len(G.edges())

def WriteGraph(G, outfile=args.outfile, graphviz=args.graphviz):
    if graphviz == "yes":
        from networkx.drawing.nx_agraph import write_dot
        try:
            write_dot(G, outfile[:-4] + ".gv")
        except:
            pass
    nx.write_edgelist(G, outfile[:-4] + ".edgelist", delimiter="\t", data=False)
    nodes = []
    for node in G.nodes():
        nodes.append([node, " ".join(map(str,G.nodes()[node]["state"]))])
    df_nodes = pd.DataFrame(nodes)
    df_nodes.to_csv(outfile[:-4] + ".nodelist", sep="\t", header=None, index=None)

if args.input_type == "cross-sectional":
    data        = LoadData(args.data, keep_nan = 1)
    obs, G      = MakeCSObservations(data)
    WriteLabels(data, outfile=args.labels, start=0)
    WriteTransitions(obs, args.outfile)
    WriteCS(data, args.outfile_cs)
    WriteGraph(G, args.outfile)
    
if args.input_type == "longitudinal":
    data           = LoadData(args.data, keep_nan = 0, index_col=None)
    obs, cs_obs, G = MakeLongitudinalObservations(data, constant_transitions=args.constant_transitions)
    WriteLabels(data, outfile=args.labels, start=1)
    WriteTransitions(obs, args.outfile)
    WriteCS(cs_obs, args.outfile_cs)
    WriteGraph(G, args.outfile)
    
if args.input_type == "sample-tree":
    G = nx.read_edgelist(args.sample_tree, create_using=nx.DiGraph())
    df = pd.read_csv(args.data, index_col=0).astype(int)
    transitions = []
    df.index = df.index.astype(str)
    df.loc["0"] = [0]*df.shape[1]
    for edge in G.edges():
        s = edge[0]
        t = edge[1]
        if t in df.index:
            if s in df.index:
                s_l = tuple(df.loc[s])
            else:
                s_l = tuple([0]*len(df.iloc[0]))
            t_l = tuple(df.loc[t])
            if s_l != t_l or args.constant_transitions == 1:
                transitions.append(s_l)
                transitions.append(t_l)
            G.nodes()[s]["state"] = s_l
            G.nodes()[t]["state"] = t_l
    df_out = pd.DataFrame(transitions)
    df_out.to_csv(args.outfile, sep=" ", header=None, index=None)
    WriteLabels(df)
    WriteCS(df, args.outfile_cs)
    WriteGraph(G, args.outfile)
    
if args.input_type == "phylogenetic":
    data          = LoadData(args.data, keep_nan = 0)
    args.resolve_single_leaves = 1
    t             = LoadTree(args.phylogeny,
                             format_type=args.phylogeny_format,
                             root_if_unrooted=1)
    t, leaf_dict          = MatchLabelsToLeaves(t, data)
    node_dict, label_dict = MakeNodeDict(t, leaf_dict, resolve_single_leaves = args.resolve_single_leaves)
    CheckLeaves(t)
    obs, G = MakeObservations(t,
                              node_dict,
                              label_dict,
                              constant_transitions=args.constant_transitions,
                              resolve_single_leaves = args.resolve_single_leaves)
    WriteLabels(data, outfile=args.labels, start=0)
    WriteTransitions(obs, args.outfile)
    WriteCS(data, args.outfile_cs)
    PrintGraphProperties(G)
    WriteGraph(G, args.outfile)
