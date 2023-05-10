# rnaseq-salmon


SRA: SRP418305
--
BioProject: PRJNA925533
--
GEO: GSE223292
--

Abstract: from the Bioproject

Retromer controls cellular homeostasis through regulating integral membrane protein sorting and transport and by 
controlling maturation of the endo-lysosomal network. Retromer dysfunction, which is linked to neurodegenerative 
disorders including Parkinson's and Alzheimer's diseases, manifests in complex cellular phenotypes, though the 
precise nature of this dysfunction, and its relation to neurodegeneration, remain unclear.
Here, we perform the first integrated multiomics approach to provide precise insight into the impact of 
Retromer dysfunction on endo-lysosomal health and homeostasis within a human neuroglioma cell model. 
We quantify profound changes to the lysosomal proteome, indicative of broad lysosomal dysfunction and inefficient 
autophagic lysosome reformation, coupled with a reconfigured cell surface proteome and secretome reflective of 
increased lysosomal exocytosis. Through this global proteomic approach and parallel transcriptomic analysis, 
we provide an unprecedented integrated view of Retromer function in regulating lysosomal homeostasis and
emphasise its role in neuroprotection.

Overall design: Comparative gene expression profiling analysis of RNA-seq data for wild-type H4 cells, 
three VPS35 KO clones, and three corresponding VPS35 KO clones rescued with VPS35-GFP expression

I first found out the genes that are differntially expressed between wild type H4 cells and VPS knockout cells.
Then found out those between VPS35 KO cells and the cells which were rescued with VPS35-GFP expression.
Then I checked their correlation. 

All the scripts for quantification, trimming, QC are added to the repo along with the plots of DGE and Gene enrichment.

-----------------------------------------------------------------------------------------------------------------------
Trimmomatic for adapter trimming

Multiqc for QC 

Salmon for  quantification of reads

DESeq2 for differntial gene expression
