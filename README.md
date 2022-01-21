# FMT-donor
This is a script repository for data analysis in the  FMT-donor paper, "Donor specific microbiome affects fecal microbiota transplantation efficacy on necrotizing enterocolitis in preterm pigs". [Preprint link](https://www.researchsquare.com/article/rs-117422/v1)

# Directory structure
```  
├── data  
│   ├── BIN_REFINEMENT  
│   │   ├── coASS  
│   │   └── singleASSBIN  
│   ├── coverm_genomewide  
│   ├── desep2_kaiju  
│   ├── panphlan  
│   └── strainphlan  
├── figure  
├── scripts_sequencing  
├── scripts_vis  
└── table  
```

**data**: metadata, processed output from *script_sequencing*, etc, are used for data visualization, generated tables  
**figure** generated figures in the paper, *png* file  
**scripts_sequencing**: processing code for sequencing data, and the working path is set accordingly on our sever  
**scripts_vis**: codes for figures and tables  
**table**: generated tables in the paper, *tsv* file  

# Detailed description
* Fig. 1 ([Inserted barplot](/figure/figure_1_barplot.png) generated by [scripts_vis/figure_1_barplot.R](/scripts_vis/figure_1_barplot.R))  
![Figure 1](/figure/figure_1.png)  

* Fig. 2 ([Figure 2](/figure/figure_2.png) generated by [scripts_vis/figure_2.R](/scripts_vis/figure_2.R))  
![Figure 2](/figure/figure_2.png)  

* Fig. 3 ([Figure 3](/figure/figure_3.png) generated by [scripts_vis/figure_3.R](/scripts_vis/figure_3.R))  
![Figure 3](/figure/figure_3.png)  

* Fig. 4 ([Figure 4](/figure/figure_4.png) generated by [scripts_vis/figure_4.R](/scripts_vis/figure_4.R))  
![Figure 4](/figure/figure_4.png)  

* Fig. 5 ([Figure 5](/figure/figure_5.png) generated by [scripts_vis/figure_5.R](/scripts_vis/figure_5.R))
![Figure 5](/figure/figure_5.png)  

* Fig. 6 ([Figure 6](/figure/figure_6.png) generated by [scripts_vis/figure_6.R](/scripts_vis/figure_6.R))  
![Figure 6](/figure/figure_6.png)  

* Fig. 7  Multi-panel figures was combined manually.  
[Fig. 7A](/figure/figure_7A.png) generated by [scripts_vis/figure_7A.R](/scripts_vis/figure_7A.R)  
[Fig. 7B](/figure/figure_7B.png) generated by [scripts_vis/figure_7B.R](/scripts_vis/figure_7B.R)  
[Fig. 7C](/figure/figure_7C.png) generated by [scripts_vis/figure_7C_table_s3.R](/scripts_vis/figure_7C_table_s3.R)  
[Fig. 7D](/figure/figure_7D.png) generated by [scripts_vis/virome.R](/scripts_vis/virome.R)  

* Fig. S1 ([Fig. S1](/figure/figure_s1.png) generated by [scripts_vis/figure_s1_s2.R](/scripts_vis/figure_s1.R))  
![Fig. S1](/figure/figure_s1.png)  

* Fig. S2 ([Fig. S2](/figure/figure_s2.png) generated by [scripts_vis/figure_s1_s2.R](/scripts_vis/figure_s1_s2.R))  
![Fig. S2](/figure/figure_s2.png)  

* Fig. S3 ([Fig. S3](/figure/figure_s3.png) generated by [scripts_vis/figure_s3.R](/scripts_vis/figure_s3.R))  
![Fig. S3](/figure/figure_s3.png)  

* Fig. S4 ([Fig. S4](/figure/figure_s4.png) generated by [scripts_vis/figure_s4.R](/scripts_vis/figure_s4.R))  
![Fig. S4](/figure/figure_s4.png)  

* Fig. S5 ([Fig. S5](/figure/figure_s5.png) generated by [scripts_vis/figure_s5.R](/scripts_vis/figure_s5.R))  
![Fig. S5](/figure/figure_s5.png)  

* Fig. S6 ([Fig. S6](/figure/figure_s6.png) generated by [scripts_vis/figure_s6.R](/scripts_vis/figure_s6.R))  
![Fig. S6](/figure/figure_s6.png)  

* Fig. S7 ([Fig. S7](/figure/figure_s7.png) generated by [MAGpy](https://github.com/WatsonLab/MAGpy) utilities [scripts_sequencing/MAGpy.sh](scripts_sequencing/MAGpy.sh)) 
![Fig. S7](/figure/figure_s7.png)  

* Fig. S8 generated with the help of [scripts_vis/virome.R](scripts_vis/virome.R) and cytoscape 

* Fig. S9 composed of several histological photos 

* Fig. S10 ([Fig. S10](/figure/figure_s10.png) was drawn manually)  
![Fig. S10](/figure/figure_s10.png)  

* Table 1 ([Table 1](/table/table_1.tsv) generated by [scripts_vis/table_1.R](/scripts_vis/table_1.R))  

* Table S1 and S2 ([Table S1](/table/table_s1.tsv) and [Table s2](/table/table_s1.tsv) generated by [scripts_vis/table_s1_s2.R](/scripts_vis/table_s1_s2.R))  
* Table S3 ([Table S3](/table/table_s3.tsv) generated by [scripts_vis/figure_7D_table_s3_s4.R](/scripts_vis/figure_7D_table_s3_s4.R))  
* Table S4 ([Table S4](/table/table_s4.tsv) generated by [scripts_vis/figure_7D_table_s3_s4.R](/scripts_vis/figure_7D_table_s3_s4.R))  











