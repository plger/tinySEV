# tinySEV

Simple shiny [SummarizedExperiment](https://bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html) Viewer.

While the examples below assume a transcriptomic dataset, this can equally well be applied to proteomics or epigenomics dataset -- or anything that can be squeezed in a SummarizedExperiment (SE).

<br/><br/>

## Preparing SummarizedExperiment objects

If you're not familiar with the SummarizedExperiment object structure, 
[first consult that documentation](https://bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html).

<br/><br/>

General information about the dataset can be saved in the metadata slots 'title', 'name', 'source', or 'description' 
(the last should contain a general description of the object, such as a brief summary of the experimental design and analysis pipeline).
For example:

```
metadata(se)$name <- "LSD RNAseq"
metadata(se)$title <- "RNAseq of rat cultures 2h following LSD"
metadata(se)$description <- "polyA-RNAseq of 15-day rat cortical cultures after 2h exposure to LSD.
  Salmon-based quantification on the Ensembl 107 transcriptome."
```

<br/><br/>

### Including differential expression analyses (DEA)

The results of differential expression analyses should be saved within one column of the `rowData` of the object,
with a column named prefixed with "DEA.". For example, assuming that `dea` is the data.frame containing your DEA results (with
feature names as row names), you can simply do:

```
rowData(se)$DEA.myAnalysis <- dea[row.names(se),]
```

Note that a whole data.frame is thus squeezed into a column of the `rowData`, meaning that you can store multiple such analyses despite 
them having the same column names. In the shiny app, that DEA analysis would now bear the name "myAnalysis".
If the name is not sufficiently clear for others to understand what kind of comparison this is and how it should be interpreted, consider 
adding a description to it. This can be done in the following way:

```
attr(rowData(se)$DEA.myAnalysis, "description") <- "Treatment versus control, correcting for sex differences"
```

DEAs are accepted in formats of common differential expression packages (e.g. edgeR, limma, DESeq2); just make sure that you have 
converted the results to a data.frame.

<br/><br/>

### Specifying default groups and annotation colors

It is possible to specify, in the object, which is the assay and the `colData` columns to use by default when plotting data from the object.
This is not constraining, i.e. one can still change it, but it prevents users from having to select them everytime.
This information is saved in the `default_view` slot of the object's metadata.
For example, suppose that the object has an assay "logcpm" and a colData variable "condition" that I want to be used by default when plotting,
this can be specified in the following way:

```
metadata(se)$default_view <- list(
  assay="logcpm",
  groupvar="condition",
  colvar="condition"
 )
```

This would result in the plots using the logCPM values, splitting samples by "condition", and coloring by conditions.
The following 'default_view' elements are recognized: assay, groupvar, colvar, gridvar (determines the faceting when 
plotting single genes, or the gaps when plotting heatmaps), top_annotation, left_annotation, etc...

#### Colors

When displaying the data along with annotations (e.g. experimental group and other covariates), the tinySEV will be assigning colors to these.
The colors might not be particularly appropriate, and they will not be constant across plots. It is therefore preferable to specify in the object
the annotation colors. This can be done simply by adding a named list of (named) color vectors to `anno_colors` slot of the object's metadata:

```
metadata(se)$anno_colors <- list(genotype=c(WT="grey", mutant="red"))
```

The colors of the heatmap itself is different, and must be passed as its own metadata slot:

```
metadata(SE)$hmcols <- c("darkred","white","darkblue")
```

See the [sechm documentation](https://bioconductor.org/packages/release/bioc/vignettes/sechm/inst/doc/sechm) for more details on this.

<br/><br/>

### Including enrichment analyses
        
Enrichment analysis results can be stored as a named list of data.frames in the metadata as follows:"

```
metadata(se)$EA <- list( analysis1=df1, analysis2=df2 )
```

The dataframes will be displayed as such and require no special formatting. However, if a column named 'genes' is present 
(containing comma-separated lists of genes), it will be used to support extra functionalities (e.g. selecting members of a geneset for plotting).


<br/><br/>

### Including pre-defined lists of genes

It is sometimes convenient to have pre-defined list of genes (e.g. differentially-expressed genes, genes from some publication, etc.) stored into the object, so that these can conveniently be plotted. Provided that the app was launched with the appropriate setting (`feature.listsTab=TRUE`), this can be done by storing a named list of genes in the `feature.lists` slot of the object's metadata, e.g.:

```
metadata(se)$feature.lists <- list(
  "My genes of interest"=c("Ezh2", "Yy1", "Kdm1a"),
  "Immediate early genes"=c("Fos","Fosb","Jun","Junb","Egr1")
)
```

Naturally, the specified gene names must correspond to row.names of the object.
