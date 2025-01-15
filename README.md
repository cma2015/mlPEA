## mlPEA: a user-friendly and multi-functionality platform for reference genome-free m<sup>6</sup>A-Seq data analysis

<a href="https://hub.docker.com/r/malab/mlpea" target="_blank"><img src="https://img.shields.io/badge/Docker_image-ready-red.svg" target="_blank"></a><a href="https://hub.docker.com/r/malab/mlpea" target="_blank"><img src="https://img.shields.io/docker/pulls/malab/mlpea"></a><a href="https://github.com/cma2015/mlPEA" target="_blank"><img src="https://img.shields.io/badge/Source%20codes-support-blue"></a>

## Introduction
- mlPEA is a user-friendly and multi-functionality platform specifically tailored to the needs of streamlined processing of m<sup>6</sup>A-Seq data in a reference genome-free manner. By taking advantage of machine learning (ML) algorithms, mlPEA enhanced the m<sup>6</sup>A-Seq data analysis by constructing robust computational models for identifying high-quality transcripts and high-confidence m<sup>6</sup>A-modified regions. mlPEA comprises four functional modules: **Data Preprocessing, Transcriptome Construction, m<sup>6</sup>A Calling, and Functional Exploration**. 
- mlPEA project is hosted on https://github.com/cma2015/mlPEA.
- mlPEA Docker image can be obtained from https://hub.docker.com/r/malab/mlpea.

## How to use mlPEA
- Test data and tutorial of mlPEA are presented at [https://github.com/cma2015/mlPEA/tree/main/User_Manual](https://github.com/cma2015/mlPEA/blob/main/User_Manual/)
- [The installation of mlPEA](https://github.com/cma2015/mlPEA/blob/main/User_Manual/tutorial/00_Installation.md)
- [The Data Preprocessing Module](https://github.com/cma2015/mlPEA/blob/main/User_Manual/01_Data_Preprocessing_Module.md)
- [The Transcriptome Construction Module](https://github.com/cma2015/mlPEA/blob/main/User_Manual//02_Transcriptome_Construction_Module.md)
- [The m<sup>6</sup>A Calling Module](https://github.com/cma2015/mlPEA/blob/main/User_Manual/03_m6A_Calling_Module.md)
- [The Functional Exploration Module](https://github.com/cma2015/mlPEA/blob/main/User_Manual/04_Functional_Exploration_Module.md)

## How to cite mlPEA
- Yang, J., Song, M., Bu, Y., Zhao, H., Liu, C., Zhang, T., Zhang, C., Xu, S., & Ma, C. (2024). mlPEA:  m<sup>6</sup>A-Seq data analysis without a reference genome. (Submitted)

## How to access help
* Comments/suggestions/bugs/issues are welcome to be reported [here](https://github.com/cma2015/mlPEA/issues) or contact: Jing Yang sxxayj@nwafu.edu.cn, Minggui Song smg@nwafu.edu.cn or Chuang Ma chuangma2006@gmail.com

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


- **Step 4**: Construct assembled transcriptome and perform ML-based transcript screening

- **Step 5**: Call m<sup>6</sup>A peak and perform ML-based peak screening

- **Step 6**: Functional Exploration


## Change log
- 2024.12 Add reference genome-free m<sup>6</sup>A-Seq data analysis pipeline in tutorial
- 2024.11 Release mlPEA v1.0
- 2024.10 mlPEA web server online
- 2022.09 We launched mlPEA project
