---
title: Quick Start
---

# <span style="color:#734f96">About this Guide</span>

This guide lets you jump right into the action with minimal preparation. If you already have the latest `gfortran` or `lfortran` compiler and `fpm`, you can continue [here](./quickstart.html#get-the-code)


[TOC]


# <span style="color:#734f96">Get the Tools</span>

## Environment manager (recommended)
The use of [conda](https://docs.conda.io) is recommended for managing your environment. While strictly speaking not essential, using conda is convenient, clean, and a common way to manage coding environments and packages. You can get the latest versions of all required packages with conda from the conda-forge channel. Once conda is installed, add the channel, then create a new environment and activate it:

```
conda config --add channels conda-forge
conda create -n fsml_session
conda activate fsml_session
```

## Compilers
FSML supports the [GFortran](https://gcc.gnu.org/fortran/) and [LFortran](https://lfortran.org/) open-source compilers, but it will most likely also compile with Intel and other compilers. (No compiler specific extensions are used). GFortran is a very mature and well-established compiler, while LFortran is a new, interactive compiler that is still under heavy development. Consequently, if you're just looking for a smooth experience, the former is recommended. On the other hand, if you're looking to use FSML in a jupyter notebook and you're comfortable with some tinkering here and there, LFortran is your choice.

You can install GFortran with conda as follows:

```
conda install -c conda-forge gfortran
```

Similarly, you can get LFortran with as follows:

```
conda install -c conda-forge lfortran
```


## FPM
The [Fortran Package Manager (fpm)](https://fpm.fortran-lang.org/) lets you setup and manage your fortran projects easily. It is modelled after Rust's Cargo, and many modern Fortran projects, like FSML, are offered as FPM packages. The following will install fpm with conda:

```
conda install -c conda-forge fpm
```


## FORD (optional - primarily for developers)
[FORD (FORtran Documenter)](https://forddocs.readthedocs.io/en/stable/) is a tool that lets you document your Fortran code easily. FSML's documentation is generated with FORD. If you wish to change and document the code, add documentation, or regenerate documentation, you can get FORD with conda:
```
conda install -c conda-forge ford
```
To generate FSML's documentation, navigate to the folder containing `doc.md` and simply run:
```
ford doc.md
```


# <span style="color:#734f96">Get the Code</span>

## Download and build with FPM
If you use FPM, you will not have to download the code manually. Instead, it can be listed as a dependency and will be downloaded and compiled automatically (see next section).

## Download manually
If you prefer to include, compile or use FSML any other way, you can download the source code directly here.



# <span style="color:#734f96">Set Up Your Project</span>

## Create a new FSML project with FPM
- Create new project.
- Add fsml to your fpm toml file. When compiling/running your project, it will be downloaded and compiled without hassle.
- Build/run the project: FPM uses gfortran by default. If you want to use FSML with lfortran, you will have to specify this through the compiler flag.

## Download manually
If you prefer to include, compile or use FSML any other way, you can download the source code directly here.



# <span style="color:#734f96">Read Data</span>


# <span style="color:#734f96">Analyse Data</span>


# <span style="color:#734f96">Model Data</span>
