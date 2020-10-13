#!/bin/bash
#pathway defination
curl http://rest.kegg.jp/list/pathway > pathway_anno.tsv
#KO defination
curl http://rest.kegg.jp/list/ko > ko_anno.tsv
#Module defination
curl http://rest.kegg.jp/list/module > module_anno.tsv
# KO2module
curl http://rest.kegg.jp/link/ko/module > ko2module.tsv
# KO2pathway
curl http://rest.kegg.jp/link/ko/pathway > ko2pathway.tsv

