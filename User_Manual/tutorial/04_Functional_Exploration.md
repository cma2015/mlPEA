<div align='center' >
<p><font size='70'><strong>mlPEA User Manual</strong></font></p>
<font size='100'>(version 1.0)</font>
</div>

- mlPEA is a user-friendly and multi-functionality platform specifically tailored to the needs of streamlined processing of m<sup>6</sup>A-Seq data in a reference genome-free manner. By taking advantage of machine learning (ML) algorithms, mlPEA enhanced the m<sup>6</sup>A-Seq data analysis by constructing robust computational models for identifying high-quality transcripts and high-confidence m<sup>6</sup>A-modified regions.
- mlPEA comprises four functional modules: **Data Preprocessing, Transcriptome Construction, m<sup>6 </sup>A Calling, and Functional Exploration**.
- mlPEA was powered with an advanced packaging technology, which enables compatibility and portability.
- The mlPEA project is hosted on [https://github.com/cma2015/mlPEA](https://github.com/cma2015/mlPEA)
- mlPEA Docker image is available at [https://hub.docker.com/r/malab/mlpea](https://hub.docker.com/r/malab/mlpea)
- The following part shows installation of mlPEA docker image and detailed documentation for each function in mlPEA.

## Functional Exploration

This module provided five functions to perform functional exploration of m<sup>6</sup>A-Seq data

| **Functions**                      |                           **Description**                           |                                    **Input**                                    |                                               **Output**                                               | Time  (test data) |                                                                                                                                                                                                                        **Reference**                                                                                                                                                                                                                        |
| ---------------------------------------- | :-----------------------------------------------------------------------: | :------------------------------------------------------------------------------------: | :----------------------------------------------------------------------------------------------------------: | ----------------- | :---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------: |
| **ML-based transcript annotation** |            Predict the coding region of assembled transcripts            |                         Assembled transcripts in FASTA format                         | Prediction scores of translation initiation and termination sites of assembled transcripts in txt format | ～1 min           |                                                                                                                                                                                                                         In-house scripts                                                                                                                                                                                                                         |
| **Differential methylation**       | Identify differential methylation modifications under multiple conditions | Assembled transcripts in FASTA format and alignment file in bam format |         RNA differential modifications in BED format         | ~5s               |                       In-house scripts                       |
| **Enrichment analysis** | Perform GO or KEGG enrichment analysis for any species |      Transcript list and annotation results      |                                The enriched GO/KEGG terms                                | ~4s            |                                                                                                                                                                                                                         <a href="https://academic.oup.com/bioinformatics/article/27/12/1653/257754" target="_blank">Timothy <I>et al </I>., 2011, Bioinformatics </a>,<a href="https://academic.oup.com/bioinformatics/article/27/12/1696/255896" target="_blank">Philip <I>et al </I>., 2011, Bioinformatics </a>,<a href="https://doi.org/10.1016/j.molcel.2010.05.004" target="_blank">Heinz <I>et al </I>., 2010, Molecular Cell </a>                                                                                                                                                                                                                         |
| **Motif discovery** |          Integrate MEME-ChIP and HOMER to performed de-novo motif discovery          |                        RNA modifications in BED format and assembled transcriptome sequences in FASTA format                        |                                          Discovered motifs in HTML format                                          | ~4s~            | <a href="https://academic.oup.com/bioinformatics/article/27/12/1696/255896" target="_blank">Philip <I>et al</I>., 2011, Bioinformatics </a> |



## **Transcript Annotation**

In this function, we utilized TranslationAI (Fan *et al*, 2023), a deep neural network to directly predict and analyze translation initiation (TIS) and termination sites (TSS) from transfrags. Functional annotation was performed useing eggNOG-mapper (v2.1.9 (Huerta-Cepas *et al*, 2017); database v5.0.2 (Huerta-Cepas *et al*, 2018)) to identify potential functions based on homology.

#### Input

- **Assembly transcript in FASTA format**

#### Output

- **TranslationAI prediction results in TXT format**
- **Transcripts annotation results in TXT format by eggnog**

#### How to use this function

- The following screenshot shows us how to use this function.

![4.1.transcript_anno.png](../img/4.1.transcript_anno.png)

## Enrichment Analysis

This function is designed to perform GO or KEGG enrichment analysis for any species through R package "clusterProfiler".

#### Input

- **Species name (Latin species name)**
- **RNA modifications gene list (a matrix seperated by TAB with one column)**

#### Output

- The enriched GO/KEGG terms
- A PDF focument of top enriched GO/KEGG terms

#### How to use this function

- The following screenshot shows us how to use this function.



## Motif **Discovery**

This function integrates MEME-ChIP and DREME to perform *de-novo* motif discovery.

#### Input

- **RNA modifications (peak regions or single nucleotide resolution) in BED format**
- **Reference genome sequences in FASTA format**

#### Output

- An HTML report generated by DREME or MEME-ChIP

#### How to use this function

- The following screenshot shows us how to use this function.

![4.3.motif.png](../img/4.3.motif.png)
