# FMT-donor
This is a repository of scripts used for FMT-donor paper "Donor specific microbiome differences affects fecal microbiota transplantation efficacy on necrotizing enterocolitis in preterm pigs". [bioaxiv link](XXXX)

# Directory structure
.  
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

**data** metadata, processed output from *script_sequencing*, etc, are used for data visualization, generated tables  
**figure** generated figures in the paper, *png* file  
**scripts_sequencing** processing code for sequencing data  
**scripts_vis** codes for figures and tables  
**table** generated tables in the paper, *tsv* file  

# Detailed description
* Figure 1 ([Inserted barplot](/figure/figure_1_barplot.png) generated by [scripts_vis/figure_1_barplot.R](/scripts_vis/figure_1_barplot.R))  
![Figure 1](/figure/figure_1.png)  

* Figure 2 ([Figure 2](/figure/figure_2.png) generated by [scripts_vis/figure_2.R](/scripts_vis/figure_2.R))  
![Figure 2](/figure/figure_2.png))  

* Figure 3 ([Figure 3](/figure/figure_3.png) generated by [scripts_vis/figure_3.R](/scripts_vis/figure_3.R))  
![Figure 3](/figure/figure_3.png))  

* Figure 4 ([Figure 4](/figure/figure_4.png) generated by [scripts_vis/figure_4.R](/scripts_vis/figure_4.R))  
![Figure 4](/figure/figure_4.png))  

* Figure 5 ([Figure 5](/figure/figure_5.png) generated by [scripts_vis/figure_5.R](/scripts_vis/figure_5.R))
![Figure 5](/figure/figure_5.png))  

* Figure 6 ([Figure 6](/figure/figure_6.png) generated by [scripts_vis/figure_6.R](/scripts_vis/figure_6.R))  
![Figure 6](/figure/figure_6.png))  

* Figure 7  Multi-panel figures was combined manually
[Figure 7A](/figure/figure_7A.png) generated by [scripts_vis/figure_7A.R](/scripts_vis/figure_7A.R))  
[Figure 7B](/figure/figure_7B.png) generated by [scripts_vis/figure_7B.R](/scripts_vis/figure_7B.R))  
[Figure 7C](/figure/figure_7C.png) generated by [MAGpy](https://github.com/WatsonLab/MAGpy) utilities [scripts_squencing/MAGpy.sh](/scripts_squencing/MAGpy.sh))  
[Figure 7D](/figure/figure_7D.png) generated by [scripts_vis/figure_7D_table_s3.R](/scripts_vis/figure_7D_table_s3.R))  

* Figure S1 ([Figure S1](/figure/figure_s1.png) was drawn manually)  
![Figure S1](/figure/figure_s1.png))  

* Figure S2 ([Figure S2](/figure/figure_s2.png) generated by [scripts_vis/figure_s2.R](/scripts_vis/figure_s2.R))  
![Figure S2](/figure/figure_s2.png))  

* Figure S3 ([Figure S3](/figure/figure_s3.png) generated by [scripts_vis/figure_s3.R](/scripts_vis/figure_s3.R))  
![Figure S3](/figure/figure_s3.png))  

* Figure S4 ([Figure S4](/figure/figure_s4.png) generated by [scripts_vis/figure_s4.R](/scripts_vis/figure_s4.R))  
![Figure S4](/figure/figure_s4.png))  

* Figure S5 ([Figure S5](/figure/figure_s5.png) generated by [scripts_vis/figure_s5.R](/scripts_vis/figure_s5.R))  
![Figure S5](/figure/figure_s5.png))  

* Figure S6 ([Figure S6](/figure/figure_s6.png) generated by [scripts_vis/figure_s6.R](/scripts_vis/figure_s6.R))  
![Figure S6](/figure/figure_s6.png))  

* Table 1  

* [Table S1](/table/table_s1.tsv) and [Table s2](/table/table_s1.tsv) generated by [scripts_vis/table_s1_s2.R](/scripts_vis/table_s1_s2.R))  

* [Table S3](/table/table_s3.tsv) generated by [scripts_vis/figure_7D_table_s3.R](/scripts_vis/figure_7D_table_s3.R))  











