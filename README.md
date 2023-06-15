# Analysis and Design of Lean Biomolecular Controllers Repo
Here you will find the supporting code for my Bioenineering Final Year Project which is my last written coursework ever. I hope this code is helpful and that you might one day find inner peace as well.

## Overview
This project was about modelling and analysing genetic circuits with biomolecular controllers in the context of cellular resources. The main idea for this repo is to introduce you to the ideas behind the model and to be able to play around with it in Julia. We will be using Pluto.jl which is an interactive environment like JuPyTer

3 main resources used in gene expression are analysed with mathematical and computational methods:
| # | Resource               | Stage of Gene Expression | Analysis Method          |
|---|------------------------|--------------------------|--------------------------|
| 1 | RNA Polymerases (RNAP) | Transcription            | Sensitivity Analysis     |
| 2 | Ribosomes              | Translation              | Sensitivity Analysis     |
| 3 | Phosphate bonds (~ATP) | Both                     | Numerical Quantification |

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
3. Domain Specific Language (DSL)
4. Simulation of Activation Cascade
5. Resource Quantification
6. Sensitivity Analysis
7. System Classification

## Quick Start Guide
1. [Install Julia][https://julialang.org/downloads/]
