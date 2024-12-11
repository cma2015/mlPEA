## mlPEA: a user-friendly and multi-functionality platform for reference genome-free m<sup>6</sup>A-Seq data analysis

<a href="https://hub.docker.com/r/malab/mlpea" target="_blank"><img src="https://img.shields.io/badge/Docker_image-ready-red.svg" target="_blank"></a><a href="https://hub.docker.com/r/malab/mlpea" target="_blank"><img src="https://img.shields.io/docker/pulls/malab/mlpea"></a><a href="https://github.com/cma2015/mlPEA" target="_blank"><img src="https://img.shields.io/badge/Source%20codes-support-blue"></a>

## Introduction
- mlPEA is a user-friendly and multi-functionality platform specifically tailored to the needs of streamlined processing of m<sup>6</sup>A-Seq data in a reference genome-free manner. By taking advantage of machine learning (ML) algorithms, mlPEA enhanced the m<sup>6</sup>A-Seq data analysis by constructing robust computational models for identifying high-quality transcripts and high-confidence m<sup>6</sup>A-modified regions. mlPEA comprises four functional modules: **Data Preprocessing, Transcriptome Construction, m<sup>6</sup>A Calling, and Functional Exploration**. 
- The mlPEA project is hosted on https://github.com/cma2015/mlPEA.
- The mlPEA Docker image can be obtained from https://hub.docker.com/r/malab/mlpea.

## How to use mlPEA
- Test data and tutorial of mlPEA are presented at https://github.com/cma2015/mlPEA/tree/main/Tutorials
- [The installation of mlPEA](https://github.com/cma2015/mlPEA/blob/main/Tutorials/00_Installation.md)
- [The Data Preprocessing Module](https://github.com/cma2015/mlPEA/blob/main/Tutorials/01_Data_Preprocessing_Module.md)
- [The Transcriptome Construction Module](https://github.com/cma2015/mlPEA/blob/main/Tutorials/02_Transcriptome_Construction_Module.md)
- [The m<sup>6</sup>A Calling Module](https://github.com/cma2015/mlPEA/blob/main/Tutorials/03_m<sup>6</sup>A_Calling_Module.md)
- [The Functional Exploration Module](https://github.com/cma2015/mlPEA/blob/main/Tutorials/04_Functional_Exploration_Module.md)

## How to cite mlPEA
- Yang, J., Song, M., Bu, Y., Liu, C., Zhao, H., Zhang, T., Zhang, C., & Ma, C. (2024). mlPEA:  m<sup>6</sup>A-Seq data analysis without a reference genome. (Submitted)

## How to access help
* Comments/suggestions/bugs/issues are welcome reported [here](https://github.com/cma2015/mlPEA/issues) or contact: Jing Yang sxxayj@nwafu.edu.cn, Minggui Song smg@nwafu.edu.cn or Chuang Ma chuangma2006@gmail.com

## Quick start

- **Step 1**: mlPEA installation from Docker Hub

```
# pull latest mlPEA Docker image from docker hub
$ docker pull malab/mlpea
```

- **Step 2**: Launch mlPEA local server

```
$ docker run -it -p 880:8080 malab/mlpea bash
$ bash /home/galaxy/run.sh
```

Then, mlPEA local server can be accessed via [http://localhost:8080](http://localhost:8080/)

- **Step 3**: Upload m<sup>6</sup>A-Seq data


- **Step 4**: Construct transcriptome and perform ML-based transcript screening

- **Step 5**: Read mapping and m<sup>6</sup>A calling

- **Step 5**: Functional Exploration


## Change log
- 2024.12 Release mlPEA v1.0
- 2024.10 mlPEA web server online
- 2022.09 we launched mlPEA project
