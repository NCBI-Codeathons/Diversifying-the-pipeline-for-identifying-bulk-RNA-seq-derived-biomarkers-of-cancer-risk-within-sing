# Diversifying the pipeline for identifying bulk RNA-seq derived biomarkers of cancer within single cell populations

Hackathon team: Sara Grimm, Jason Wang, Miko Liu, Matt Bernstein

## Background and Objective
Previous work by Matt Bernstein (and others) explored scRNA-seq data of 8 high-grade glioma tumor samples from "Single-cell transcriptome analysis of lineage diversity in high-grade glioma" by Yuan et al (PMID: 30041684). See [paper](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-018-0567-9) and [repository](https://github.com/NCBI-Codeathons/Identifying-bulk-RNA-seq-derived-biomarkers-of-cancer-risk-within-single-cell-populations) for their results.  

Our primary objective is to extend their analysis by developing a method to stratify cells in a given scRNA-seq dataset according to malignancy status. For this stratification we rely primarily on a CNV (copy number variation) metric.  

We also hoped to (a) identify gene markers or gene sets correlated with malignancy status, and (b) associate specific cell types from clustered scRNA-seq data with malignancy status.  Unfortunately, despite several attempted methods, we were not successful in these endeavors.


## Workflow

![image](https://user-images.githubusercontent.com/46359281/76649230-9de6a800-6536-11ea-9458-55a6e5440f0c.png)

## Workflow Steps and Code Bits
Step 1:  Seurat

For each tumor sample, run Seurat to extract the following:  normalized read counts (SCT) per gene per cell, assigned cluster per cell, UMAP coordinates, average expression per gene per cluster, and marker genes per cluster.  The R workflow using Seurat can be found at working_data/seurat_simple/{sampleID}/protocol-{sampleID}.Rtxt

Step 2:  Process inferCNV data.

(Thanks to Matt B for providing inferCNV outputs for each of the 8 tumor samples using [InferCNV](https://github.com/broadinstitute/inferCNV/wiki).) From the ~.observations.txt and ~.references.txt files output from inferCNV, the summed and average |locusScore-medianScore| were calculated over all available loci to provide aggregate CNV metrics per cell.
```
  id = samples[S];
  txtfile=paste(id, ".aggr_medianDelta_per_cell.txt", sep="");
  out=paste("cellID", "cellGroup", "aggregate", "locusCt", "median", sep="\t");
  write.table(out, file=txtfile, sep="\t", quote=F, append=FALSE, row.names=F, col.names=F);
  infileR=paste(id, ".infercnv.references.txt", sep="");
  dataR=read.delim(infileR, header=TRUE, row.names=1);
  nC=ncol(dataR);  # cells
  nL=nrow(dataR);  # loci
  medianR=median(unlist(dataR));
  for (N in 1:nC) {
    dd=dataR[,N];
    delta=abs(dd-medianR);
    deltaSum=sum(delta);
    thiscell=colnames(dataR)[N];
    out=paste(thiscell, "reference", deltaSum, nL, medianR, sep="\t");
    write.table(out, file=txtfile, sep="\t", quote=F, append=TRUE, row.names=F, col.names=F);
  }
  infileT=paste(id, ".infercnv.observations.txt", sep="");
  dataT=read.delim(infileT, header=TRUE, row.names=1);
  nC=ncol(dataT);  # cells
  nL=nrow(dataT);  # loci
  medianT=median(unlist(dataT));
  for (N in 1:nC) {
    dd=dataT[,N];
    delta=abs(dd-medianT);
    deltaSum=sum(delta);
    thiscell=colnames(dataT)[N];
    out=paste(thiscell, "sample", deltaSum, nL, medianT, sep="\t");
    write.table(out, file=txtfile, sep="\t", quote=F, append=TRUE, row.names=F, col.names=F);
  }
```

Then the distribution of the average |locusScore-medianScore| results were plotted to see if the reference and tumor samples were different. (They are!) After repeating over the 8 tumors, we noted that the distribution of these scores in the reference cells are consistently below 0.02, so we are using this as the threshold for what we are confident(???) are non-malignant cells.  These plots can be viewed in working_data/inferCNV_distr.
```
  id = samples[S];
  infile=paste(id, ".aggr_medianDelta_per_cell.txt", sep="");
  data=read.delim(infile, header=TRUE);
  ddR=subset(data, (data$cellGroup == "reference"));
  ddS=subset(data, (data$cellGroup == "sample"));
  loci=ddR[1,4];
  xmin=0; xmax=0.08;
  pngfile=paste(id, ".avg_medianDelta_per_cell.distr.png", sep="");
  png(pngfile, h=500, w=500, res=120);
  maintitle=paste(id, ": average |score-median|\n(", loci, " inferCNV loci)", sep="");
  plot(density(ddR$aggregate/loci), lwd=2, col="black", main=maintitle, xlab="", xlim=c(xmin,xmax), ylim=c(0,yy[S]));
  lines(density(ddS$aggregate/loci), lwd=2, col="limegreen");
  legend("topright", c("reference","sample"), col=c("black","limegreen"), pch=20);
  dev.off();
```

Step 3:  Combine Seurat cell cluster and inferCNV data.

We then integrated the cell cluster assingments from Seurat with the aggregate inferCNV score we calculated, then visualized via UMAP to look at relative location of low/high CNV cells compared to cells with high expression of proliferation marker genes.

## Screenshot Examples and Output  

*yellow: non-malignant*  
*blue: malignant*  

![PJ016 UMAP-byRankICNV](https://user-images.githubusercontent.com/46359281/76648676-81963b80-6535-11ea-808e-295a022e9360.png)![UMAP_proliferation_PJ016](https://user-images.githubusercontent.com/46359281/76647682-a8ec0900-6533-11ea-9ecc-8ffc142e61a1.png)  

![PJ017 UMAP-byRankICNV](https://user-images.githubusercontent.com/46359281/76648677-81963b80-6535-11ea-878e-51725a7b53f4.png)![UMAP_proliferation_PJ017](https://user-images.githubusercontent.com/46359281/76648202-a211c600-6534-11ea-9481-16b5151bae1f.png)  

![UMAP_proliferation_PJ018](https://user-images.githubusercontent.com/46359281/76648204-a211c600-6534-11ea-899e-a069f6c5455c.png)![PJ018 UMAP-byRankICNV](https://user-images.githubusercontent.com/46359281/76648679-822ed200-6535-11ea-8cd2-fd8a74857478.png)  

![PJ025 UMAP-byRankICNV](https://user-images.githubusercontent.com/46359281/76648681-822ed200-6535-11ea-8e6e-5907a6940e08.png)![UMAP_proliferation_PJ025](https://user-images.githubusercontent.com/46359281/76648205-a211c600-6534-11ea-90e8-cadd1a8d6d59.png)  

![PJ030 UMAP-byRankICNV](https://user-images.githubusercontent.com/46359281/76648671-80fda500-6535-11ea-878a-5d09d8f85636.png)![UMAP_proliferation_PJ030](https://user-images.githubusercontent.com/46359281/76648207-a211c600-6534-11ea-8cd8-bf939087f68d.png)  

![PJ032 UMAP-byRankICNV](https://user-images.githubusercontent.com/46359281/76648672-80fda500-6535-11ea-8ac2-e973d972d973.png)![UMAP_proliferation_PJ032](https://user-images.githubusercontent.com/46359281/76648209-a211c600-6534-11ea-909e-3f8ecb604f17.png)  

![PJ035 UMAP-byRankICNV](https://user-images.githubusercontent.com/46359281/76648674-81963b80-6535-11ea-8630-e27917f0fb12.png)![UMAP_proliferation_PJ035](https://user-images.githubusercontent.com/46359281/76648210-a2aa5c80-6534-11ea-866c-150849da5fe5.png)  

![PJ048 UMAP-byRankICNV](https://user-images.githubusercontent.com/46359281/76648675-81963b80-6535-11ea-8b7e-41eb5411fe57.png)![UMAP_proliferation_PJ048](https://user-images.githubusercontent.com/46359281/76648211-a2aa5c80-6534-11ea-8484-35f661fc2699.png)

## Additional Visualizations  

## Annotation of cell clusters using PlangaoDB
![image](https://user-images.githubusercontent.com/46359281/76651763-8d84fc00-653b-11ea-9300-b06bd366bc3c.png)  
## Cells with aggregate inferCNV score in range of reference (assumed non-malignant (red)) cells.
![image](https://user-images.githubusercontent.com/46359281/76653076-5ebc5500-653e-11ea-84d8-ca3fc856294e.png)  
## Expression of selected markers
![biomarker_expression-SOX2](https://user-images.githubusercontent.com/46359281/76650803-b60bf680-6539-11ea-9990-ad2af685d789.png) ![biomarker_expression-OLIG2](https://user-images.githubusercontent.com/46359281/76652567-3a13ad80-653d-11ea-87a8-d4871adf57f0.png)



## Future Directions

- Moving forward, we would like to add gene set enrichment analysis to use enriched gene sets consistent with the caner phenotype, including proliferation and developmental pathways. This would provide a more rigorous statistical score for each cluster as opposed to the proliferation score currently used.
- Furthermore, testing other tools for annotating gene clusters with their respective cell types would aid in identifying non-malignant clusters.
- Examine specific genes for abnormal CNV from the inferred CNV data (TP53, MDM2, RTK, RAS, PI3K, RB, CDK4)
- Use identified cancer cells on the scRNA-Seq dataset integrating all 8 donors to identify an improved boundary (such as using supporter vector machines) rather than the principal components 1 & 2 defined in the paper.

## Dependencies
- Seurat

## Input data  

Find the scRNA-seq data and inferred copy number variant data in input/data. Within working_data, several visualizations have been generated from the scRNA-Seq UMAPs mapping inferred CNV onto the sRNA-Seq UMAPs (umap_color_by_inferCNV_rank, umap_color_by_inferCNV_assumedNonMalignant), proliferation score, and biomarker expression.


