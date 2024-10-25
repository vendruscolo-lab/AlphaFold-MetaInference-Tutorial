# Case study: TDP-43 monomer dynamics

Amyotrophic lateral sclerosis (ALS) and frontotemporal lobar degeneration (FTLD) are two of the most common debilitating neurodegenerative diseases (NDs) in people aged between 45 and 65 years. Thepathological hallmark of these diseases is the presence of proteinaceous cytoplasmic inclusions in degenerating neurons. TDP-43 is the primary component of the cytoplasmic inclusions in ~97$\%$ of ALS and ~45$\%$ of FTLD cases. Increasing evidence suggests that cellular changes leading to TDP-43 mislocalization and aggregation in the cytoplasm resulting in the gain of toxic functions that drive neurodegeneration in ALS and FTLD-TDP[ref](https://doi.org/10.1038/s41593-023-01341-4). TDP-43 is equipped The sequence of TDP-43 comprises 414 amino acids, which form different domains (Figure 3). These domains include a folded N-terminal domain (residues 1-76), a disordered region (residues 77-105), a folded RNA recognition motif (residues 106-176), a second disordered region (residues 177-190), another folded RNA recognition motif (residues 191-259), and an long disordered C-terminal domain (residues 274-414), which contains a glycine-rich region, is involved in protein-protein interactions, and harbors most of the mutations associated with ALS[ref](https://pubmed.ncbi.nlm.nih.gov/33177049/). 

Due to the interplay between well folded domains and disordered regions, TDP-43 monomer challenges experimental structure determination approaches. In particular a recent study of the WtoA TDP-43 revealed the SAXS intensity profile of TDP-43 monomer[Ref](https://www.cell.com/iscience/fulltext/S2589-0042(20)30344-8?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2589004220303448%3Fshowall%3Dtrue). Due to its societal relevance and challenging target for experimental structural biology, we have selected TDP-43 WtoA as a protein to simulate its structural ensemble by AF-MI and compare to experimental data in the original AF-MI [paper](https://www.biorxiv.org/content/10.1101/2023.01.19.524720v1.full).


<p align="center">
  <img src="https://github.com/riccardocapelli/VMetaD-tutorial/blob/main/img/lys_render.jpg?raw=true" alt="Alt text" width="25%">
  <br>
  <em>Render of the Lysozyme:benzene in cartoon representation. In dark gray we can see the N-terminal domain, in blue the C-terminal domain, where in red we can see benzene in its binding pocket.</em>
</p>

In this tutorial, we will use a structure based on the PDB ID [1L84](https://www.rcsb.org/structure/1L84). The force field used is multi-eGO, a structure-based potential generated via a bayesian approach, trained on all-atom simulations performed with the DES-Amber force field (details in [this paper](https://doi.org/10.26434/chemrxiv-2024-jcmgc)). 
The advantages of multi-eGO rely in three main points
* The possibility to use 5 fs time steps: multi-eGO implements a united-atom coarse-graining, which removes hydrogens;
* Absence of coulombian interactions for electrostatics: in multi-eGO, all the non-bonded interactions are included in a Lennard-Jones potential;
* Implicit solvation: multi-eGO does not have an explicit solvent, resulting in ~10x less particles in the simulation box.

Using this potential we are able to reach on consumer-level hardware the same performances of a HPC machine on all-atom potentials, allowing us to run on a workstation the entire tutorial in around 24 h. 

##### [Back to AlphaFold-Metainference home](NAVIGATION.md)