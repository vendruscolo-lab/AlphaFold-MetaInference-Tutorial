# AlphaFold-Metainference tutorial
A tutorial to setup and run AlphaFold-Metainference simulations. The original work which explains the theory for AlphaFold-Metainference (AF-MI) can be found [here](https://www.biorxiv.org/content/10.1101/2024.11.09.622758v1).

The tutorial comprises the following sections:


#### [AlphaFold-Metainference theory](theory.md)
A brief introduction to the theory and advantage in using AlphaFold-Metainference.

#### [Case study: TDP-43 WtoA monomer dynamics](case_study.md)
We will present the case study for this tutorial, the DNA/RNA processing protein TDP-43 which undergoes both functional and pathogennic aggregation in amyotrophic lateral sclerosis (ALS) and other neurodegenerative conditions. To enable a rapid structural ensemble generation on a laptop/workstation, AF-MI is coupled to the CALVADOS-2 coarse-grained forcefield.

#### [Input files preparation](input.md)
We will discuss all the preliminary steps to set up AF-MI, from generating the AF2 distance-map restraints data, to preparing PLUMED inputes, and running energy minimization and production runs

#### [Structural ensemble generation and post-processing analysis](postprocessing.md)
Finally, we will discuss how to generate the unbiased structural ensemble of TDP-43, estimate free energy surfaces along collective variables (CVs) and predict ensemble averaged properties such as secondary structureand RMSF.

___NB:___ This tutorial assumes that you know metadynamics theory and practice in PLUMED (if not, please visit [PLUMED tutorial](https://www.plumed-tutorials.org/lessons/21/004/data/NAVIGATION.html)).


```mermaid
flowchart TB
  A[PLUMED syntax] -.-> D[Theory];
  B[Parallel Bias Metadynamics] -.->  D;
  C[Metainference] -.-> D;
  D --> E[Case study];
  E --> F[Jupyter/Colab]
  F --> G[Inputs]
  G --> F[Jupyter/Colab]
  F --> H[Post-processing]

  click A "syntax" "This lesson teaches you the basic features of the PLUMED input syntax"
  click B "PBMetaD" "This lesson teaches you how to perform and analyze Metadynamics simulations"
  click C "Metainference" "This lesson teaches you how to perform and analyze Metadynamics simulations"
  click D "theory.md" "AlphaFold-Metainference theory intro & overview"
  click E "case_study.md" "Case study: TDP-43 monomer dyamics"
  click F "AF-MIcolab" "AF-MI colab performing software installation, system preparation, AF-MI simulation and post-processing analysis"
  click G "input.md" "Input: Protein sequence,  AF-MI simulation data and parameters"
  click H "postprocessing.md" "Structural ensemble generation and post-processing analysis"
```
