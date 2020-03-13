# Diversifying the pipeline for identifying bulk RNA-seq derived biomarkers of cancer within single cell populatins

Hackathon team: Sara Grimm, Jason Wang, Miko Liu, Matt Bernstein

## Background and Objective
Previous work by Matt Bernstein (and others) explored scRNA-seq data of 8 high-grade glioma tumor samples from "Single-cell transcriptome analysis of lineage diversity in high-grade glioma" by Yuan et al (PMID: 30041684). See https://github.com/NCBI-Codeathons/Identifying-bulk-RNA-seq-derived-biomarkers-of-cancer-risk-within-single-cell-populations for their results.  

Our primary objective is to extend their analysis by developing a method to stratify cells in a given scRNA-seq dataset according to malignancy status. For this stratification we rely primarily on a CNV (copy number variation) metric.  

We also hoped to (a) identify gene markers or gene sets correlated with malignancy status, and (b) associate specific cell types from clustered scRNA-seq data with malignancy status.  Unfortunately, despite several attempted methods, we were not successful in these endeavors.


## Workflow  
![image](https://user-images.githubusercontent.com/46359281/76649230-9de6a800-6536-11ea-9458-55a6e5440f0c.png)


## screenshot examples and output
![PJ016 UMAP-byRankICNV](https://user-images.githubusercontent.com/46359281/76648676-81963b80-6535-11ea-808e-295a022e9360.png)![UMAP_proliferation_PJ016](https://user-images.githubusercontent.com/46359281/76647682-a8ec0900-6533-11ea-9ecc-8ffc142e61a1.png)  

![PJ017 UMAP-byRankICNV](https://user-images.githubusercontent.com/46359281/76648677-81963b80-6535-11ea-878e-51725a7b53f4.png)![UMAP_proliferation_PJ017](https://user-images.githubusercontent.com/46359281/76648202-a211c600-6534-11ea-9481-16b5151bae1f.png)  

![UMAP_proliferation_PJ018](https://user-images.githubusercontent.com/46359281/76648204-a211c600-6534-11ea-899e-a069f6c5455c.png)![PJ018 UMAP-byRankICNV](https://user-images.githubusercontent.com/46359281/76648679-822ed200-6535-11ea-8cd2-fd8a74857478.png)  

![PJ025 UMAP-byRankICNV](https://user-images.githubusercontent.com/46359281/76648681-822ed200-6535-11ea-8e6e-5907a6940e08.png)![UMAP_proliferation_PJ025](https://user-images.githubusercontent.com/46359281/76648205-a211c600-6534-11ea-90e8-cadd1a8d6d59.png)  

![PJ030 UMAP-byRankICNV](https://user-images.githubusercontent.com/46359281/76648671-80fda500-6535-11ea-878a-5d09d8f85636.png)![UMAP_proliferation_PJ030](https://user-images.githubusercontent.com/46359281/76648207-a211c600-6534-11ea-8cd8-bf939087f68d.png)  

![UMAP_proliferation_PJ032](https://user-images.githubusercontent.com/46359281/76648209-a211c600-6534-11ea-909e-3f8ecb604f17.png)![PJ032 UMAP-byRankICNV](https://user-images.githubusercontent.com/46359281/76648672-80fda500-6535-11ea-8ac2-e973d972d973.png)  

![UMAP_proliferation_PJ035](https://user-images.githubusercontent.com/46359281/76648210-a2aa5c80-6534-11ea-866c-150849da5fe5.png)  ![PJ035 UMAP-byRankICNV](https://user-images.githubusercontent.com/46359281/76648674-81963b80-6535-11ea-8630-e27917f0fb12.png)

![UMAP_proliferation_PJ048](https://user-images.githubusercontent.com/46359281/76648211-a2aa5c80-6534-11ea-8484-35f661fc2699.png)![PJ048 UMAP-byRankICNV](https://user-images.githubusercontent.com/46359281/76648675-81963b80-6535-11ea-8b7e-41eb5411fe57.png)







## future directions

## Dependencies

## input and output





building on : https://github.com/NCBI-Codeathons/Identifying-bulk-RNA-seq-derived-biomarkers-of-cancer-risk-within-single-cell-populations
