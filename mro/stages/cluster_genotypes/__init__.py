#!/software/python-2.7.6/bin/python
import subprocess
import glob
import gzip
import vcf
import pyfasta
import json
import random
import numpy as np
from sklearn.cluster import SpectralClustering
import pandas as pd

def split(args):
    chunks = []
    for vcf in args.vcfs:        
        chunks.append({'vcf':vcf, '__mem_gb':8})
    return {'chunks':chunks,'join':{'__mem_gb':8}}

def getGC(fasta, chrom, pos, window):
    seq = fasta[chrom][max(0,pos-window/2):min(len(fasta[chrom]),pos+window/2)]
    gc = 0
    for c in seq.upper():
        if c == 'G' or c == 'C':
            gc += 1
    return gc/float(len(seq))


def main(args, outs):
    args.coerce_strings()
    outs.coerce_strings()
    num_calls_used = 0
    fasta = pyfasta.Fasta(args.fasta)
    
    num_calls = 0
    metrics = {}
    outs.metrics = metrics
    with open(args.vcf,'r') as vcf_file:
        vcf_reader = vcf.Reader(vcf_file, compressed=True)
        with open(outs.calls[:-3], 'w') as calls:
            for rec in vcf_reader:
                if not rec.is_snp:
                    continue
                num_calls +=  1
                if str(rec.ALT[0]) == "<*>":
                    continue                
                if rec.QUAL < args.qual_filter:
                    continue
                if getGC(fasta, rec.CHROM, rec.POS, 30) <= 0.15:
                    continue
                assert(fasta[rec.CHROM][rec.POS-1] == rec.REF)
                if fasta[rec.CHROM][rec.POS-2] == rec.ALT[0] and fasta[rec.CHROM][rec.POS] == rec.ALT[0]:
                    continue
                refs = rec.get_hom_refs()#[call.sample for call in rec.get_hom_refs()]
                alts = rec.get_hom_alts()#[call.sample for call in rec.get_hom_alts()]
                hets = rec.get_hets()#[call.sample for call in rec.get_hets()]
                ref = len(refs)
                alt = len(alts)
                if args.ploidy > 1:
                    ref += len(hets)
                    alt += len(hets)
                alts_split = [[] for x in range(len(rec.ALT))]
                for alt_rec in alts:
                    alts_split[int(alt_rec.data.GT[0])-1].append(alt_rec)
                for het_rec in hets:
                    alts_split[int(het_rec.data.GT[-1])-1].append(het_rec)
                    
                if ref >= args.ref_alleles_required and alt >= args.alt_alleles_required:
                    num_calls_used += 1
                    alt_string = ""
                    for index,alt_s in enumerate(rec.ALT):
                        alt_string += str(alt_s)
                        if index < len(rec.ALT) -1:
                            alt_string += ","
                    ref_cells = [call.sample for call in refs]
                    for het in hets:
                        ref_cells.append(het.sample)
                    if len(alt_string) > 1:
                        continue
                    calls.write(rec.CHROM+"\t"+str(rec.POS)+"\t.\t"+rec.REF+"\t"+alt_string+"\t"+str(rec.QUAL)+"\t.\t")
                    calls.write(",".join(ref_cells))
                    for alt_recs in alts_split:
                        alt_cells = [call.sample for call in alt_recs]
                        calls.write(";"+",".join(alt_cells))
                    calls.write("\n")
    metrics["num_calls"] = num_calls
    metrics["num_calls_used"] = num_calls_used
    metrics["alt_alleles_required"] = args.alt_alleles_required
    metrics["ref_alleles_required"] = args.ref_alleles_required
    outs.metrics = metrics
    make_graph(outs.calls[:-3], outs.graph, outs.graph[:-4]+".abc", args.reads_per_cell, outs.graph+"alleles.csv", args.downsample_graph, args.ground_truth)
    #make_graph(outs.calls, outs.graph, args.reads_per_cell, outs.graph+"alleles.csv", args.downsample_graph, None, args.transcript_clusters)

def make_graph(calls_fn, dot_out, abc_out, reads_per_cell, allele_counts_fn, downsample = None, ground_truth = None, transcript_clusters = None):
    reads_per_cell_counts = {}
    graph_stats = {}
    all_nodes = set()
    with open(reads_per_cell) as rpc:
        for line in rpc:
            tokens = line.strip().split(',')
            reads_per_cell_counts[tokens[0]] = int(tokens[1])
            all_nodes.add(tokens[0])

    if downsample == None:
        nodes_to_keep = all_nodes
    else:
        random.seed(downsample)
        nodes_to_keep = set(random.sample(list(all_nodes), min(len(all_nodes),downsample)))

    graph = {} # map from cell1:cell2:count
    bad_graph = {} # map of discordant alleles from cell1:cell2:count
    with open(calls_fn) as calls: 
        with open(allele_counts_fn,'w') as allele_counts_file:
            for line in calls:
                tokens = line.strip().split()
                alleles = tokens[7]
                bxs = alleles.split(";")#[1:] # SKIPPING REF ALLELES
                cells = [] # cells is [[cell1_allele1, cell2_allele2,...],[cell1_allele2,cell2_allele2],...]
                for allele in bxs:
                    cells.append(allele.split(",")) 
                indel = "indel" if not(len(tokens[3]) == len(tokens[4])) else "snp"
                if len(cells) == 2:
                    allele_counts_file.write(indel + "," + str(len(cells[0]))+","+str(len(cells[1]))+"\n")
                min_allele_frequency = 1.0
                for (allele_index, allele_cell_list) in enumerate(cells):
                    allele_frequency = len(allele_cell_list)/float(len(cells[1-allele_index]))
                    allele_frequency = max(0.1, allele_frequency)
                    allele_frequency = min(0.9, allele_frequency)
                    min_allele_frequency = min(min_allele_frequency, allele_frequency)
                    for index1, cell1 in enumerate(allele_cell_list):
                        for index2 in range(index1+1,len(allele_cell_list)):
                            cell2 = allele_cell_list[index2]
                            if cell1 < cell2:
                                cella = cell1
                                cellb = cell2
                            else:
                                cella = cell2
                                cellb = cell1
                            edges = graph.setdefault(cella, {})
                            count = edges.setdefault(cellb, 0)
                            #if allele_index == 0:
                            #    graph[cella][cellb] += 1.0/6.0 
                            #else:
                            #    graph[cella][cellb] += 1
                            graph[cella][cellb] += 1.0 - allele_frequency
                for index in range(len(cells)):
                    for index2 in range(index+1,len(cells)):
                        if index == index2:
                            continue
                        for cell1 in cells[index]:
                            for cell2 in cells[index2]:
                                if cell1 < cell2:
                                    cella = cell1
                                    cellb = cell2
                                else:
                                    cella = cell2
                                    cellb = cell1
                                #if cella in graph and cellb in graph[cella]:
                                #    graph[cella][cellb] -= 1.0 - min_allele_frequency
        
    iqr = sorted(reads_per_cell_counts.values())
    p25 = iqr[int(0.25*len(reads_per_cell_counts))]
    p50 = iqr[int(0.5*len(reads_per_cell_counts))]
    p75 = iqr[int(0.75*len(reads_per_cell_counts))]
    
    metrics = {}
    metrics["p25_reads_per_cell"] = p25
    metrics["p50_reads_per_cell"] = p50
    metrics["p75_reads_per_cell"] = p75
    graph_stats['reads_per_cell_iqr'] = metrics
    graph_stats['nodes'] = len(iqr)
    edges = 0
    edge_weights = [0]


    colors = ["#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#a65628","#f781bf","#f781bf"]
    if ground_truth:
        ground_truth_mapping = {}
        samples_so_far = {}
        with open(ground_truth) as gt:
            gt.readline()
            for line in gt:
                tokens = line.strip().split()
                assignment = tokens[5].split("-")
                if not assignment[1] in samples_so_far:
                    samples_so_far[assignment[1]] = len(samples_so_far)        
                if assignment[0] == "SNG":
                    ground_truth_mapping[tokens[0]] = "color=\""+colors[samples_so_far[assignment[1]]]+"\""
                else:
                    if not assignment[2] in samples_so_far:
                        samples_so_far[assignment[2]] = len(samples_so_far)
                    #ground_truth_mapping[tokens[0]] = "style=filled fillcolor=\""+colors[samples_so_far[assignment[1]]]+"\" color=\""+colors[samples_so_far[assignment[2]]]+"\""
                    ground_truth_mapping[tokens[0]] = "color=\""+colors[samples_so_far[assignment[1]]]+"\""
                    
                
    with open(abc_out,'w') as abc:
        with open(dot_out,'w') as dot:
            dot.write("graph cells {\n")
            dot.write("scale=8.0\n")
            nodes = set()
            for node1, node2s in graph.iteritems():
                nodes.add(node1)
                for node2, count in node2s.iteritems():
                    nodes.add(node2)
            for node in nodes:
                if not node in nodes_to_keep:
                    continue
                color = "color=red"
                count = reads_per_cell_counts[node]
                if count > p25:
                    color = "color=orange"
                if count > p50: 
                    color = "color=yellow"
                if count > p75:
                    color = "color=green"

                if ground_truth:
                    color = ground_truth_mapping[node]
                dot.write(node.replace('-','.')+" [shape=point, "+color+"]\n")
                #dot.write("node [shape=point];\n")

            for node1, node2s in graph.iteritems():
                for node2, count in node2s.iteritems():
                    if count < 2:
                        continue
                    if not node1 in nodes_to_keep:
                        continue
                    if not node2 in nodes_to_keep:
                        continue
                    dot.write(node1.replace('-','.')+" -- "+node2.replace('-','.')+" [color=black, weight="+str(count)+"];\n")
                    abc.write(node1+"\t"+node2+"\t"+str(count)+"\n")
                    if node1 in all_nodes:
                        all_nodes.remove(node1)
                    if node2 in all_nodes:
                        all_nodes.remove(node2)
                    edges += 1
                    edge_weights.append(count)
            dot.write("}\n")
            graph_stats['singletons'] = len(all_nodes)
            graph_stats['edges'] = edges
            edge_weights = sorted(edge_weights)
            edge_weight_dist = {}
            edge_weight_dist['p25'] = edge_weights[int(0.25*len(edge_weights))]
            edge_weight_dist['p50'] = edge_weights[int(0.5*len(edge_weights))]
            edge_weight_dist['p75'] = edge_weights[int(0.75*len(edge_weights))]
            graph_stats['edge_weight_dist'] = edge_weight_dist
            edge_weight_counts = {}
            for edge_weight in edge_weights:
                if edge_weight < 5000:
                    edge_weight_counts.setdefault(int(edge_weight),0)
                    edge_weight_counts[int(edge_weight)] += 1
            graph_stats['edge_weight_histogram'] = edge_weight_counts
    return graph_stats

from plotnine import *
import matplotlib
    
                    

def join(args, outs, chunk_defs, chunk_outs):
    metrics = {}
    usable_calls = []
    for chunk_out in chunk_outs:
        chunk_out.coerce_strings()
        usable_calls.append(chunk_out.calls[:-3])
        for key, val in chunk_out.metrics.iteritems():
            metrics.setdefault(key, 0)
            metrics[key] += val
    metrics["alt_alleles_required"] = args.alt_alleles_required
    metrics["ref_alleles_required"] = args.ref_alleles_required
    with open(outs.calls[:-3]+"unsorted.tsv",'w') as calls:
        command = ['cat']
        command.extend(usable_calls)
        #print command
        subprocess.check_call(command, stdout = calls)
        #for calls_file in usable_calls:
        #    with open(calls_file) as cfile:
        #        for line in cfile:
        #            calls.write(line)
    with open(outs.calls[:-3],'w') as calls:
        subprocess.check_call(['sort', '-k1,1', '-k2,2n', outs.calls[:-3]+"unsorted.tsv"],stdout=calls)
    metrics["graph_stats"] = make_graph(outs.calls[:-3], outs.graph, outs.graph[:-4]+".abc", args.reads_per_cell, outs.graph+"alleles.csv", args.downsample_graph, args.ground_truth)
    #subprocess.check_call(['mcl',outs.graph[:-4]+".abc",'--abc','-o',outs.clusters, '-tf', '#ceilnb(300),#knn(200)'])
    adjacency_matrix, edges, edges_by_node, node_list, node_index = compute_adjacency_matrix(outs.graph[:-4]+".abc")
    clustering_labels = None
    best_score = -10000
    if args.num_clusters == None:
        min_clusters = args.min_clusters
        max_clusters = args.max_clusters
        if min_clusters == None:
            min_clusters = 2
        if max_clusters == None:
            max_clusters = 20
        for num_clusters in range(min_clusters, max_clusters+1):
            spectral_clustering = SpectralClustering(num_clusters, affinity='precomputed', n_init=1000, assign_labels='discretize')
            spectral_clustering.fit(adjacency_matrix)
            #print "clustering k = "+str(i)
            score = clustering_score(edges,spectral_clustering.labels_,node_index)
            if score >= best_score:
                best_score = score
                clustering_labels = spectral_clustering.labels_  
            #score2 = sil_score(edges_by_node, sc.labels_,node_index)
            #print "\t"+str(score)+"\t"+str(score2)
    else:
        spectral_clustering = SpectralClustering(args.num_clusters, affinity='precomputed', n_init=1000, assign_labels='discretize')
        spectral_clustering.fit(adjacency_matrix)
        clustering_labels = spectral_clustering.labels_

    create_cluster_output(outs.clusters, clustering_labels, node_index)
    outs.metrics = metrics
    subprocess.check_call(['bgzip',outs.calls[:-3]])
    subprocess.check_call(['tabix','-p', 'vcf',outs.calls])    
    max_clusters = 10
    clusters, cluster_names, cluster_total, cell_clusters = load_clustering(outs.clusters, max_clusters)
    #clusters, cluster_names, cell_clusters = load_clustering(outs.clusters, max_clusters)
    if not args.ground_truth == None: 
        cluster_concordance(args.ground_truth, outs.concordance, outs.concordance_bargraph, outs.concordance_heatmap, max_clusters, clusters, cluster_names, cell_clusters, cluster_total)
    #if not args.transcription_pca is None:
    #    write_pca_colored_by_cluster(args.transcription_pca, outs.pca_vs_clusters, outs.pca_vs_clusters_data, cell_clusters)
    with open(outs.metrics_json,'w') as mets:
        json.dump(metrics, mets)


import operator 
def create_cluster_output(out, clustering_labels, node_index):
    clusters_to_nodes = {}
    for node, index in node_index.iteritems():
        cluster = clustering_labels[index]
        clusters_to_nodes.setdefault(cluster, [])
        clusters_to_nodes[cluster].append(node)
    with open(out,'w') as outwrite:
        for cluster, node_list in sorted(clusters_to_nodes.items(), key= lambda kv: len(kv[1]),reverse=True):
            outwrite.write("\t".join(node_list)+"\n")

def clustering_score(edges,labels,node_index):
    weight_col = []
    cluster1_col = []
    cluster2_col = []
    clusters = set()
    cluster_num_nodes = {}
    for cluster in labels:
        cluster_num_nodes.setdefault(cluster,0)
        cluster_num_nodes[cluster] += 1
    for (node1, node2), weight in edges.iteritems():
        node1_cluster = labels[node_index[node1]]
        node2_cluster = labels[node_index[node2]]
        cluster1_col.append(node1_cluster)
        cluster2_col.append(node2_cluster)
        weight_col.append(weight)
        clusters.add(node1_cluster)
        clusters.add(node2_cluster)
    clusters = [c for c in clusters]
    df = pd.DataFrame({'weight':weight_col,'c1':cluster1_col,'c2':cluster2_col})
    intra_cluster_weights = {}
    for cluster in clusters:
        df2 = df[np.logical_and(df.c1 == cluster, df.c2 == cluster)]
        total = np.sum(df2.weight)
        denom = cluster_num_nodes[cluster]
        denom *= denom-1
        denom /= 2
        #print "\tintra cluster cluster, total, denom, avg, cluster_size"
        #print "\t"+str(cluster)+"\t"+str(total)+"\t"+str(denom)+"\t"+str(total/float(max(denom,1.0)))+"\t"+str(cluster_num_nodes[cluster])
        mean_weight = total/float(max(denom,1.0))
        intra_cluster_weights[cluster] = mean_weight
    inter_cluster_weights = {}
    for index1, cluster1 in enumerate(clusters):
        for cluster2 in clusters[index1+1:]:
            df2 = df[np.logical_or(np.logical_and(df.c1 == cluster1, df.c2 == cluster2),
                                np.logical_and(df.c1 == cluster2, df.c2 == cluster1))]
            #df2 = df[np.logical_or(np.logical_and(df.c1 == cluster1, np.logical_not(df.c2 == cluster1)),
            #                   np.logical_and(np.logical_not(df.c1 == cluster1), df.c2 == cluster1))]
            total = np.sum(df2.weight)
            #denom = cluster_num_nodes[cluster1]*np.sum([cluster_num_nodes[c] for c in clusters if not c == cluster1])
            denom = cluster_num_nodes[cluster1]*cluster_num_nodes[cluster2]
            #print "\t\tinter cluster cluster1, cluster2, total, denom, avg, cluster_size"
            #print "\t\t"+str(cluster1)+"\t"+str(cluster2)+"\t"+str(total)+"\t"+str(denom)+"\t"+str(total/float(max(denom,1.0)))+"\t"+str(cluster_num_nodes[cluster2])
            mean_weight = total/float(max(denom,1.0))
            inter_cluster_weights.setdefault(cluster1,[])
            inter_cluster_weights.setdefault(cluster2,[])
            inter_cluster_weights[cluster1].append(mean_weight)
            inter_cluster_weights[cluster2].append(mean_weight)
    cluster_scores = {}
    for cluster in clusters:
        cluster_scores[cluster] = (intra_cluster_weights[cluster]-np.max(inter_cluster_weights[cluster])) / max(1.0,intra_cluster_weights[cluster])
    #print cluster_scores
    total_score = 0.0#1000000
    denom = 0.0
    for cluster, score in cluster_scores.iteritems():
        total_score += score*cluster_num_nodes[cluster]#score#min(score,total_score) #* cluster_num_nodes[cluster]
        #print (score, cluster_num_nodes[cluster], score*cluster_num_nodes[cluster])
        denom += cluster_num_nodes[cluster]
    total_score /= denom
    #score = np.mean([x for _, x in cluster_scores.iteritems()])
    #print score
    return total_score

def compute_adjacency_matrix(abc_graph):
    edges = {}
    edges_by_node = {}
    nodes = set()
    with open(abc_graph) as graph:
        for line in graph:
            tokens = line.strip().split("\t")
            weight = float(tokens[2])
            node1 = tokens[0]
            node2 = tokens[1]
            edges[(node1,node2)] = weight
            edges_by_node.setdefault(node1,[])
            edges_by_node.setdefault(node2,[])
            edges_by_node[node1].append((node2,weight))
            edges_by_node[node2].append((node1,weight))
            nodes.add(tokens[0])
            nodes.add(tokens[1])
    node_index = {node:index for (index, node) in enumerate(nodes)}
    node_list = [node for node in nodes]
    adj_mat = np.zeros((len(node_index), len(node_index)))
    for (node1, node2), weight in edges.iteritems():
        node1_index = node_index[node1]
        node2_index = node_index[node2]
        adj_mat[node1_index][node2_index] = weight
        adj_mat[node2_index][node1_index] = weight
    return adj_mat, edges, edges_by_node, node_list, node_index

def write_pca_colored_by_cluster(pca, pca_vs_clusters, pca_vs_clusters_data, cell_clusters):
    with open(pca) as projection:
        with open(pca_vs_clusters_data,'w') as out:
            out.write(projection.readline().strip()+",cluster\n")
            for line in projection:
                cell = line.strip().split(",")[0]
                if cell not in cell_clusters: 
                    continue
                cluster = cell_clusters[cell]
                out.write(line.strip()+","+str(cluster)+"\n")
    df = pd.DataFrame(pca_vs_clusters_data)
    #pca_graph = ggplot(df)+geom_point(aes(x='PC-1',y='PC-2',color='cluster'))
    #pca_graph.write(pca_vs_clusters)

def cluster_concordance(ground_truth, concordance_out, barchart, heatmap, max_clusters, clusters, cluster_names, cell_clusters, cluster_total):
    cluster_assignment_counts, assignment_names, assignment_total = \
    compare_ground_truth(ground_truth, cluster_names, cell_clusters, max_clusters)
    with open(concordance_out, 'w') as concordance:
        concordance.write("ground_truth,cluster,count,assignment_total,cluster_total\n")
        for cluster, assignment_counts in cluster_assignment_counts.iteritems():
            for assignment, count in assignment_counts.iteritems():
                concordance.write(str(assignment)+","+str(cluster)+","+str(count)+","+str(assignment_total[assignment])+","+str(cluster_total[cluster])+"\n")
    df = pd.read_csv(concordance_out)
    #barchart_obj = ggplot(df)+geom_bar(aes(x='ground_truth',y='count',color='cluster'),position='dodge',stat='identity')
    #barchart_obj.save(barchart)
    df['count/cluster_total'] = df['count']/df.cluster_total
    #heatmap_obj = ggplot(df)+geom_tile(aes(x='ground_truth',y='cluster',fill='count/cluster_total'))
    #heatmap_obj.save(heatmap)

def compare_ground_truth(ground_truth, cluster_names, cell_clusters, max_clusters):
    assignment_names = set()
    assignment_total = {}
    cluster_assignment_counts = {}
    with open(ground_truth) as demux:
        demux.readline() # get rid of header
        for line in demux:
            tokens = line.strip().split()
            assignment = tokens[6]
            cell = tokens[0]
            if not cell in cell_clusters:
                continue
            cluster = cell_clusters[cell]
            if cluster > max_clusters:
                cluster = "unknown_cluster"
            cluster_assignment_counts.setdefault(cluster,{})
            cluster_assignment_counts[cluster].setdefault(assignment, 0)
            cluster_assignment_counts[cluster][assignment] += 1
            assignment_names.add(assignment)
            assignment_total.setdefault(assignment,0)
            assignment_total[assignment] += 1
        for cluster in cluster_names:
            for assignment in assignment_names:
                cluster_assignment_counts.setdefault(cluster,{})
                cluster_assignment_counts[cluster].setdefault(assignment,0)
    return cluster_assignment_counts, assignment_names, assignment_total
            

def load_clustering(mcl_file, max_clusters):
    cluster_names = set()
    clusters = {}
    cluster_total = {}
    cell_clusters = {}
    with open(mcl_file) as mcl:
        for index, line in enumerate(mcl):
            clusters[index] = set()
            for cell in line.strip().split():
                clusters[index].add(cell)
                if index <= max_clusters:
                    cluster_names.add(index)
                    cluster = index
                else:
                    cluster = "unknown_cluster"
                    cluster_names.add(cluster)
                cell_clusters[cell] = cluster
                cluster_total.setdefault(cluster,0)
                cluster_total[cluster] += 1
    return clusters, cluster_names, cluster_total, cell_clusters
