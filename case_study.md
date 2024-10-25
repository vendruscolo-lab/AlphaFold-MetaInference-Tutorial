# Case study: TDP-43 WtoA monomer dynamics

Amyotrophic lateral sclerosis (ALS) and frontotemporal lobar degeneration (FTLD) are two of the most common debilitating neurodegenerative diseases (NDs) in people aged between 45 and 65 years. Thepathological hallmark of these diseases is the presence of proteinaceous cytoplasmic inclusions in degenerating neurons. TDP-43 is the primary component of the cytoplasmic inclusions in ~97$\%$ of ALS and ~45$\%$ of FTLD cases. Increasing evidence suggests that cellular changes leading to TDP-43 mislocalization and aggregation in the cytoplasm resulting in the gain of toxic functions that drive neurodegeneration in ALS and FTLD-TDP[ref](https://doi.org/10.1038/s41593-023-01341-4). TDP-43 is equipped The sequence of TDP-43 comprises 414 amino acids, which form different domains (Figure 3). These domains include a folded N-terminal domain (residues 1-76), a disordered region (residues 77-105), a folded RNA recognition motif (residues 106-176), a second disordered region (residues 177-190), another folded RNA recognition motif (residues 191-259), and an long disordered C-terminal domain (residues 274-414), which contains a glycine-rich region, is involved in protein-protein interactions, and harbors most of the mutations associated with ALS[ref](https://pubmed.ncbi.nlm.nih.gov/33177049/). 

Due to the interplay between well folded domains and disordered regions, TDP-43 monomer challenges experimental structure determination approaches. In particular a recent study of the WtoA TDP-43 revealed the SAXS intensity profile of TDP-43 monomer[Ref](https://www.cell.com/iscience/fulltext/S2589-0042(20)30344-8?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2589004220303448%3Fshowall%3Dtrue). Due to its societal relevance and challenging target for experimental structural biology, we have selected TDP-43 WtoA as a protein to simulate its structural ensemble by AF-MI and compare to experimental data in the original AF-MI [paper](https://www.biorxiv.org/content/10.1101/2023.01.19.524720v1.full).


<p align="center">
  <img src="https://github.com/vendruscolo-lab/AlphaFold-MetaInference-Tutorial/blob/main/images/TDP-43WtoA.jpg?raw=true" alt="Alt text" width="45%">
  <br>
  <em>Render of the TDP-43WtoA AlphaFold structure. Plddt structure prediction accuracy score spans from blue (high accuracy) to red (low accuracy). Noted the low PLDDT score in the dirordered C-terminal region 274-414 as well as in the linkers between folded domains, signifying inaccuracy of AF structure prediction in such dynamic regions. </em>
</p>

In this tutorial, we start with the AF predicted structure of TDP-43 WtoA, employ CALVADOS-2 coarse grained forcefield while adding RMSD restraints to maintain the folded domains and use AF inter-residue distances as restraints within the Metainference framework. We then generate the unbiased structural ensemble, by properly accounting the MetaD bias, followed by ploting free energy surface along different collective variables, and calculation of the secondary structure elements and RMSF.  



##### [Back to AlphaFold-Metainference home](NAVIGATION.md)
