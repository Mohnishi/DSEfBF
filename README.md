## Readme

This code is meant to reproduce the results found in **Dynamic Structure Estimation from Bandit Feedback** by Motoya Ohnishi, Isao Ishikawa, Yuko Kuroki, and Masahiro Ikeda.

## Setup & Install

This code has been tested on Ubuntu 20.04 LTS, but should also work on different platforms (MacOS, Windows, FreeBSD) if the instructions are adapted.

The process to bring up this repo is as follows:
1. Download and install [Julia](https://julialang.org/)
2. Navigate to project and instantiate 
3. Run

The following is an example of installing Julia for Ubuntu 20.04.
```bash
cd ~/Downloads
wget https://julialang-s3.julialang.org/bin/linux/x64/1.6/julia-1.6.3-linux-x86_64.tar.gz
tar xvf julia-1.6.3-linux-x86_64.tar.gz

# the following exports can be added to your bashrc.
export JULIA_BINDIR=~/Downloads/julia-1.6.3/bin
export PATH=$JULIA_BINDIR:$PATH
export JULIA_NUM_THREADS=12

cd $directory_you_extracted_code
julia
```

Once you start Julia, regardless of platform, the following instructions may proceed:
```julia
julia> ]
(@v1.6) pkg> registry add https://github.com/JuliaRegistries/General
(@v1.6) pkg> registry add https://github.com/Lyceum/LyceumRegistry     # add Lyceum registry
(@v1.6) pkg> activate .    # activates this project
(DynamicStructureEstimation) pkg> instantiate   # the built in package manager downloads, installs dependences
(DynamicStructureEstimation) pkg> ctrl-c

julia> include("main.jl")    # this will run all of the experiments and display the results at the end; it will take only a few minutes.
```

## Notes

The results in the paper were generated with **Julia 1.6.3**, with **12 Julia threads**. This is critical to reproducibility, but not necessary for running the included algorithm; one should adapt these settings to their compute.

Also, **one may need to restart Julia to run experiments sequentially**.  To exit julia, do
```julia
julia> exit()
```
Next time you start Julia, you do not need to do instantiate but only activate.

## Code Structure

```bash
.
├── algs            # algorithm files
│   ├── eigenestimation.jl    
│   └── periodestimation.jl   
├── log             # data store
│   ├── eigen
│   ├── muperiod
│   └── period
├── main.jl         # the main file you execute
├── Project.toml    # Julia Project file for top level dependencies
├── README.md       # this file
├── scripts         # main files for estimations executed by main.jl
│   ├── eigen.jl    
│   ├── period.jl   
│   └── period_2.jl 
└── utils           # models for the experiments
    ├── ca.jl       
    ├── linear.jl   
    └── mup.jl      


```
## Code Maintenance

The codes are maintained by the authors of **Dynamic Structure Estimation from Bandit Feedback**.
Visit the [project page](https://sites.google.com/view/dsefbf/) for arXiv and Video links.
