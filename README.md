# Analysis and Design of Lean Biomolecular Controllers Repo
Here you will find the supporting code for my Bioenineering Final Year project which is my last written coursework ever. I hope this code is helpful and that you might one day find inner peace as well.

## Overview
This project was about modelling and analysing genetic circuits with biomolecular controllers in the context of cellular resources. 3 main resources used in gene expression are analysed with mathematical and computational methods:
| # | Resource               | Stage of Gene Expression | Analysis Method          |
|---|------------------------|--------------------------|--------------------------|
| 1 | RNA Polymerases (RNAP) | Transcription            | Sensitivity Analysis     |
| 2 | Ribosomes              | Translation              | Sensitivity Analysis     |
| 3 | Phosphate bonds (~ATP) | Both                     | Numerical Quantification |

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
2. Node 
