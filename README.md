# Analysis and Design of Lean Biomolecular Controllers Repo
Here you will find the supporting code for my Bioengineering Final Year Project which is my last written coursework ever. The original code is in the [DSL_backend](DSL_backend) folder which contains the pure Julia files. I have formatted the background required and code in a readable way in the form of Jupyter notebooks for your viewing pleasure. I hope this code is helpful and that you might one day attain this level of happiness as well.

## Overview
This project was about modelling and analysing genetic circuits with biomolecular controllers in the context of cellular resources. The main idea for this repo is to introduce you to the ideas behind the model and to be able to play around with it in Julia. We will be using Jupyter Notebooks that you can also run in your favourite IDE (I use Visual Studio Code).

3 main resources used in gene expression are analysed with mathematical and computational methods:
| # | Resource               | Stage of Gene Expression | Analysis Method          |
|---|------------------------|--------------------------|--------------------------|
| 1 | RNA Polymerases (RNAP) | Transcription            | Sensitivity Analysis     |
| 2 | Ribosomes              | Translation              | Sensitivity Analysis     |
| 3 | Phosphate bonds (~ATP) | Both                     | Numerical Quantification |

## Dose Response with/without Resource Conservation
Summary of results of framework: 

Simplified framework with is able to reflect experimental results for the dose response of an activation cascade from the paper, [Resource Competition Shapes the Response of Genetic Circuits](https://pubs.acs.org/doi/10.1021/acssynbio.6b00361) by Qian. *et al*. 

See [4_activation_cascade_simulation](4_activation_cascade_simulation.ipynb). It is important to note that the interest is in the underlying dynamics not the actual quantity.

![Simulation with Resource Conservation](simulation_plots/resource-conservation-simulation.png)
![Simulation without Resource Conservation](simulation_plots/no-resource-conservation-simulation.png)

### Illustration of resource use in central dogma:
<pre>
          Transcription         Translation  
               (TX)                 (TL)  
+------------+       +------------+       +------------+  
|            |       |            |       |            |  
|  Gene, g   +--+--> |   mRNA, m  +--+--> | Protein, x |  
|            |  |    |            |  |    |            |  
+------------+  |    +------------+  |    +------------+  
                |                    |  
              RNAP, [p]            Ribosome, [r]  
              + [ATP]              + [ATP]  
</pre>

## Contents
1. Mathematical Model Definition
2. Node Visualisation in Matrix Form
3. Domain Specific Language (DSL) Showcase
4. Simulation of Activation Cascade ([Original Julia File](DSL_backend/original_cascade_simulation.jl))
5. Resource Quantification
6. Sensitivity Analysis
7. DSL_backend

## Quick Start Guide
1. [Install Julia](https://julialang.org/downloads/)
2. Open your command line
3. Run the Julia REPL
```bash
user> julia
```
4. Add Julia to Jupyter Notebook: Enter `]` to enter Julia's package manager
```julia
julia> ]
(v1.7 pkg)> add IJulia
```
5. Download the [Anaconda](https://www.anaconda.com/download/) Python distribution and open it to launch Jupyter Notebook from the Anaconda Navigator.
6. Clone this repo or just download it.
7. Navigate to the downloaded folder and open the files with Jupyter Notebook.

More info about Julia basics: https://sje30.github.io/catam-julia/intro/julia-manual.html

## Repo cloning instructions:
1. [Install Git](https://github.com/git-guides/install-git) 
2. Open up the folder you want to download the files in.
3. Open up command line
4. Type `cd` then drag the folder into the command line and press `enter`. This will navigate to the folder you want to download the files in easily.
5. Copy and paste the following into the command line
```bash
user> git clone https://github.com/AlexanderSimWS/Y4_FYP.git
```
Continue from point 7. in the Quick Start Guide.
