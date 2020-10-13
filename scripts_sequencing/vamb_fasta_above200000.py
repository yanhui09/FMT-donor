#!/usr/bin/env python3
import sys
sys.path.append('/home/yanhui/Data/depot_tool/vamb_2.0.1/vamb')
import vamb
def splitclusters(clusters, sep):
    # First split by sample
    split = dict()
    for n, (clustername, contigs) in enumerate(clusters.items()):
        for contig in contigs:
            samplename = contig.partition(sep)[0]
            newclustername = 'c' + str(n) + '_' + 's' + samplename
            if newclustername in split:
                split[newclustername].add(contig)
            else:
                split[newclustername] = {contig}
    return split

def filterclusters(clusters, lengthof):
    # Now filter away the small bins
    filtered_bins = dict()
    for medoid, contigs in clusters.items():
        binsize = sum(lengthof[contig] for contig in contigs)
    
        if binsize >= 200000:
            filtered_bins[medoid] = contigs
    
    return filtered_bins

#load contig names, might take some time
with open('contigs/long_scaffolds.fasta', 'rb') as contigfile:
    tnfs, contignames, lengths = vamb.parsecontigs.read_contigs(contigfile)

# This read cluster.tsv file
with open('output_v2/clusters.tsv') as file:
    cluster_iterator= vamb.cluster.read_clusters(file)
clusters = dict(cluster_iterator)


lengthof = dict(zip(contignames, lengths))
filtered_bins = filterclusters(splitclusters(clusters, '_'), lengthof)
print('Number of bins before splitting and filtering:', len(clusters))
print('Number of bins after splitting and filtering:', len(filtered_bins))

#####personal input file

# This writes a .tsv file with the clusters and corresponding sequences
with open('output_v2/clusters_filtered.tsv', 'w') as file:
    vamb.cluster.write_clusters(file, filtered_bins)

# Only keep contigs in any filtered bin in memory
keptcontigs = set.union(*filtered_bins.values())

with open('contigs/long_scaffolds.fasta', 'rb') as file:
    fastadict = vamb.vambtools.loadfasta(file, keep=keptcontigs)

bindir = 'output_v2/vamb_bins'
vamb.vambtools.write_bins(bindir, filtered_bins, fastadict, maxbins=1100)

