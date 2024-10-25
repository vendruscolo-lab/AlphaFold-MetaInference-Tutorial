# AlphaFold-Metainference theory

***AlphaFold distance map prediction.*** Average inter-residue distances of disordered proteins
can be predicted through the distogram head of AlphaFold using this colab [notebook](https://github.com/zshengyu14/ColabFold_distmats/blob/main/AlphaFold2.ipynb). These
distances are defined as those between the β carbon atom positions for
all amino acids except glycine, for which the α carbon atom positions
were instead used. The multiple sequence alignment (MSA) was conducted
by MMseqs2 (default setting) on BFD/MGnify and Uniclust30.
Model 1.1.1 of AlphaFold (default setting) was used for the
predictions, with no structural templates. AlphaFold describes the
distribution of inter-residue distances into 64 bins of equal width,
covering the range from 2.15625 to 21.84375 Å, with the last bin also
including distances longer than 21.84375 Å. For each pair of residues
($i$ and $j$), AlphaFold predicts the probability $p_{ij}^{b}$ that
their distance is within bin $b$. The predicted distance
${\widehat{d}}_{ij}$ of the
predicted distribution of the distances between residue $i$ and $j$ are
esimated as


$\hat{d_{ij}}=\sum_{b=1}^{64} d^{b} p_{ij}^{b} $


where $d^{b}$ represents the central value of bin $b$.

***Metainference*.** Metainference is a Bayesian inference method that
enables the determination of structural ensembles by combining prior
information and experimental data according to the maximum entropy
principle^23^. In this work, we implemented this method by using the
distance matrix $\mathbf{d}^{AF}$ predicted by AlphaFold as
pseudo-experimental data. By design, metainference can disentangle
structural heterogeneity from systematic errors, such as force field or
forward model inaccuracies, random errors in the data, and errors due to
the limited sample size of the ensemble^23^. The molecular simulations
are carried out according to the metainference energy function,
$E = - k_{B}T\log\left( p_{MI} \right)$, where $k_{B}$ is the Boltzmann
constant, *T* is the temperature, and $p_{MI}$ is the metainference,
maximum-entropy-compatible, posterior probability

$p_{MI}\left( \mathbf{X},\sigma^{SEM},\sigma^{B} \right|\mathbf{D}) = \ \prod_{r = 1}^{N_{R}}{p\ \left( X_{r} \right)\prod_{i = 1}^{N_{D}}{p(\mathbf{D|\ X},\sigma_{i}^{SEM},\sigma_{r,i}^{Β})p(\sigma_{r,i})\ \ \ \ \ \ \ \ \ \ }}$
(3)

In this formula, **X** denotes the vector comprising the atomic
coordinates of the structural ensemble, consisting of individual
replicas X~r~ (N~R~ in total), σ^SEM^ the error associated to the
limited number of replicas in the ensemble, σ~B~ the random and
systematic errors in the prior molecular dynamics forcefield as well as
in the forward model and the data, and $\mathbf{d}^{AF}$ the AF distance
matrix. Note that σ^SEM^ is calculated for each data point (σ~i~^SEM^),
while σ^B^ is computed for each data point *i* and replica *r* as
σ~r,i~^B^. The functional form of the likelihood
p($\mathbf{d}^{AF}$\|**X**, σ~i~^SEM^ , σ~r,i~^B^) is a Gaussian
function

$$p(\mathbf{d}^{AF}\mathbf{|\ X},\sigma_{i}^{SEM},\sigma_{r,i}^{Β}) = \frac{1}{\sqrt{2\pi}\ \sqrt{\left( \sigma_{r,i}^{Β} \right)^{2} + \left( \sigma_{i}^{SEM} \right)^{2}}}exp\left\lceil - \frac{1}{2}\frac{\left\lbrack d_{i,j}^{AF} - d_{ij}\left( \mathbf{X} \right) \right\rbrack^{2}}{\left( \sigma_{r,i}^{Β} \right)^{2} + \left( \sigma_{i}^{SEM} \right)^{2}} \right\rceil\ \ \ (4)$$

where d~i,j~(X) represents the forward model for data point *i,j,*
namely the *i,j* distance calculated in the ensemble. For multiple
replicas, the metainference energy function is

$$E_{MI}\left( \mathbf{X},\sigma \right) = E_{MD}\left( \mathbf{X} \right) + \frac{k_{B}T}{2}\sum_{r,i}^{N_{R},N_{D}}\frac{\left\lbrack d_{i} - f_{i}\left( X_{r} \right) \right\rbrack^{2}}{\left( \sigma_{r,i}^{Β} \right)^{2} + \left( \sigma_{i}^{SEM} \right)^{2}} + E_{\sigma}\ \ \ (5)$$

where E~σ~ corresponds to the energy term associated with the errors

$$E_{\sigma} = k_{B}T\sum_{r,i}^{N_{R},N_{D}}{- \log{\left( \sigma_{r,i}^{Β} \right) + \frac{1}{2}}\log\left\lbrack \left( \sigma_{r,i}^{Β} \right)^{2} + \left( \sigma_{i}^{SEM} \right)^{2} \right\rbrack}\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ (6).$$

Finally, E~MD~ corresponds to the potential energy function of the
molecular dynamics force field, which in this case is the CALVADOS-2
force field^44^. While the space of conformations X~r~ is sampled
through multi-replica simulations (in this study we used six replicas)
the error parameters for each datapoint σ~r,i~^B^ are sampled through a
Gibbs sampling scheme at each time step^23^. The range of the error
sampling was \[0.0001,10\] and the associated trial move error
perturbation of the Gibbs sampling was 0.1. The error parameter due to
the limited number of replicas used to estimate the forward model
(σ^SEM^) was calculated on the fly by window averaging every 200 steps
of molecular dynamics.

**AlphaFold predicted distance selection used as restraints in AlphaFold
metainference.** In Fig.SX1 plot the inter-residue AlphaFold predicted
distance maps corresponding to the systems presented in the main text.
As mentioned in the main text and shown in **Figure S1.panel2**, when
excluding AF distance pairs with $p_{i,j}(r \geq 21.84\ Å\ ) > 0.02$, we
are left with good correlation of AF distances with MD and C2
back-calculated distances. Moreover, AlphaFold Predicted Alignment Error
(PAE), can also provide with a metric of accuracy/fluctuations of a
distance pair. PAE estimates the expected positional error for each
residue in a predicted protein structure if it were aligned to a
corresponding residue in the true protein structure (see FigSX3). Hence
we carry on a sensitivity analysis throughout the 11 C2 ensembles and
quantify the sweet-spot of minimal PAE cuttoff that increases the R^2^
to the C2 back-calculated distances (see **Figure S2**). Moreover, we
report in **Figure S3** the agreement of each C2 ensemble to the
experimental SAXS data by means of D~KL~. We observe that for all
systems, the R^2^ saturates from a given PAE cuttoff onwards which in
turn signifies that adding more distances with more fluctuations and
higher predicted alignment error does not increase the correlation to
the C2 ensemble back-calculated distances.

Using as benchmarking comparison a subset of 6 proteins from Tesei et.
al. Nature 2024 where we have available experimental Rg data and
CALVADOS-2 has larger deviations from experimental Rg, we provide with a
systematic manner of selecting AF distances as restraints in Fig. SX.4.
This approach shows than when a protein has a Kyte-doolittle hydropathy
of less than -1.4, signifying hydrophilic character and at least a 5
residue stretch with plddt\>75, then selecting AF distances with PAE\<10
has better agreement with experimental Rgs (see NLS, Msh6-NRT).
Otherwise a PAE\<5 is used (see AN16,CTCF-RD,Eralpha-NTD,CortactinCRH).
CALVADOS-2 is trained based on residue hydrophobicity and alphafold
suffers in predicting structures of low plddt regions that are equipped
with flexible residues. The motivation for the combined hydropathy and
plddt criterion stems from the fact that if any corrections need to be
applied either to calvados or to alphafold so that they generate better
ensembles, these will be where a hydrophobicity-based model alone faces
challenges in predicting interactions. In Fig. SX4 we show that
CALVADOS-2 can be assisted by AF-MI in generating ensembles that better
agree with Rg data.

Therefore, we use as restraints a subset of the full
$\mathbf{d}^{AF}$matrix by considering distances between amino acids at
least 3 residues apart, and with a Predicted Alignment Error (PAE) under
5 Å or 10 Å according to the abovementioned distance selection rule and
excluding $p_{i,j}(r \geq 21.84\ Å\ ) > 0.02$. Also, based on the fact
that C2 coarse grained model used in the AF-MI simulations does maintain
secondary structure, we employed the predicted local distance difference
test (pLDDT) score in AlphaFold to select the inter-residue distances
predicted with higher confidence and hence predicted to correspond to
structured regions. Sequence regions of at least two residues with a
pLDDT score \>0.75 were hence considered structured regions and were
restrained to the AlphaFold-predicted structure by using an upper root
mean square distance (RMSD) wall. The residue distances corresponding to
these structured regions where excluded from the distance restraints
set. This action may not be necessary if one uses all-atom force fields
that are able to maintain the protein secondary structure. The PLUMED
input files for AF-MI can be found in
(<https://zenodo.org/record/7756138>). The code to prepare an AF-MI
simulation can be found in (Github
<https://github.com/vendruscolo-lab/AlphaFold-IDP>)

***Generation of structural ensembles using metainference.*** To sample
efficiently the conformational space of disordered proteins using
average inter-residue distances predicted by AlphaFold as structural
