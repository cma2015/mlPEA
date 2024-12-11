## PanGraphRNA: an efficient and flexible bioinformatics platform for graph pangenome-based RNA-seq data analysis

<a href="https://hub.docker.com/r/malab/pangraphrna" target="_blank"><img src="https://img.shields.io/badge/Docker_image-ready-red.svg" target="_blank"></a><a href="https://hub.docker.com/r/malab/pangraphrna" target="_blank"><img src="https://img.shields.io/docker/pulls/malab/pangraphrna"></a><a href="https://github.com/cma2015/PanGraphRNA" target="_blank"><img src="https://img.shields.io/badge/Source%20codes-support-blue"></a>

## Introduction
- PanGraphRNA is an efficient, flexible and web-based Galaxy platform that can be easily used to construct graph pangenomes from genetic variations at individual, subpopulation, and population levels. It can assist researchers to select appropriate graph pangenomes using various performance metrics for both real and simulation experiments. Currently, PanGraphRNA is composed of four graph pangenome-related functional modules: **Preparation Module, Construction Module, Evaluation Module, and Application Module**. 
- The PanGraphRNA project is hosted on https://github.com/cma2015/PanGraphRNA.
- The PanGraphRNA Docker image can be obtained from https://hub.docker.com/r/malab/pangraphrna.
![PanGraphRNA](https://github.com/cma2015/PanGraphRNA/blob/main/Tutorials/PanGraphRNA_framework.png)
## How to use PanGraphRNA
- Test data and tutorial of PanGraphRNA are presented at https://github.com/cma2015/PanGraphRNA/tree/main/Tutorials
- [The installation of PanGraphRNA](https://github.com/cma2015/PanGraphRNA/blob/main/Tutorials/00_Installation.md)
- [The Graph Pangenome Preparation Module](https://github.com/cma2015/PanGraphRNA/blob/main/Tutorials/01_Graph_Pangenome_Preparation_Module.md)
- [The Graph Pangenome Construction Module](https://github.com/cma2015/PanGraphRNA/blob/main/Tutorials/02_Graph_Pangenome_Construction_Module_and_Alignment.md)
- [The Graph Pangenome Evaluation Module](https://github.com/cma2015/PanGraphRNA/blob/main/Tutorials/03_Graph_Pangenome_Evaluation_Module.md)
- [The Graph Pangenome Application Module](https://github.com/cma2015/PanGraphRNA/blob/main/Tutorials/04_Graph_Pangenome_Application_Module.md)

## How to cite PanGraphRNA
- Bu, Y., Qiu, Z., Sun, W., Han, Y., Liu, Y., Yang, J., Song, M., Li, Z., Zhang, Y., & Ma, C. (2024). PanGraphRNA: an efficient and flexible bioinformatics platform for graph pangenome-based RNA-seq data analysis. (Submitted)

## How to access help
* Comments/suggestions/bugs/issues are welcome reported [here](https://github.com/cma2015/PanGraphRNA/issues) or contact: Yifan Bu b761307648@163.com, Minggui Song smg@nwafu.edu.cn or Chuang Ma chuangma2006@gmail.com

## Change log
- 2024.10 Release PanGraphRNA v1.0
- 2022.09 we launched PanGraphRNA project

## Quick start

- **Step 1**: PanGraphRNA installation from Docker Hub

```
# pull latest PanGraphRNA Docker image from docker hub
$ docker pull malab/pangraphrna
```

- **Step 2**: Launch PanGraphRNA local server

```
$ docker run -it -p 880:8080 malab/pangraphrna bash
$ bash /home/galaxy/run.sh
```

Then, PanGraphRNA local server can be accessed via [http://localhost:8080](http://localhost:8080/)
![quick_start](https://github.com/cma2015/PanGraphRNA/blob/main/Tutorials/Figure/Figure0.jpg)

- **Step 3**: Upload RNA-seq, reference genome and variation data


- **Step 4**: Construct graph pangenome (e.g. individual level graph pangenome) and perform read-genome alignment

![quick_start2](https://github.com/cma2015/PanGraphRNA/blob/main/Tutorials/Figure/Figure0_1.jpg)
