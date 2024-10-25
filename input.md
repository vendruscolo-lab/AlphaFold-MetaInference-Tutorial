# Input files preparation

In the zip file available in the [GitHub repository](https://github.com/compsb-unimi/VMetaD-tutorial) you will find all the initial files needed to run this tutorial.

## Preliminary molecular dynamics run
As a first step for the VMetaD run, we have to choose the atoms which will constitute the reference frame on which we will calculate the position of the ligand with respect to the protein. To be sure in choosing residues that remain stable (i.e. do not belong to a moving loop), we run a 100 ns-long plain molecular dynamics (MD) simulation at 300 K (that is the temperature at which we will run and analyze all the simulations in this tutorial). This run took 2 h on 4 CPUs of a M1 Max MacBook Pro).

After the simulation, we calculate the per-residue root mean square fluctuations (RMSF), which highlights the more dynamic residues.

<p align="center">
  <img src="https://github.com/riccardocapelli/VMetaD-tutorial/blob/main/img/rmsf.jpg?raw=true" alt="Alt text" width="50%">
  <br>
  <em>Per-residue root mean square fluctuation (RMSF) computed from a plain MD simulation of lysozyme-benzene complex. The portion of the plot with light blue background represents the C-terminal domain, where the benzene binds.</em>
</p>

We can see that almost all the C-terminal domain (residues 71-162) does not show large fluctuations, and only few residues are above the (arbitrary) RMSF threshold of 0.15 nm. We thus consider to define the reference frame considering all the residues 71-162 but the five ones with RMSF > 0.15 nm (residues 135, 136, 137, 140, and 141).

## PDB file for reference frame refitting
To be sure to keep the structure stable during the simulation, we will define the reference frame of VMetaD on the residues selected for the fit. To have a reference structure to be used in PLUMED (see the next section for the input file), we convert the initial structure to the PDB format using `gmx editconf` (to ensure consistency with the numbering) and keeping only the C-alphas:
```
gmx editconf -f starting.gro -o ref_ca.pdb
grep "CA" ref_ca.pdb > tmpfile && mv tmpfile ref_ca.pdb
```
Remember to add the correct values to the b-factor columns (and/or removing unwanted atoms for the fitting procedure) telling PLUMED which are the atoms we want to use for the fitting procedure (see the `ref_ca.pdb` file in the folder for reference).

## Choice of the restraining potential size

We now need to choose the size of the restraint potential. In the [original paper](https://doi.org/10.1021/acs.jpclett.9b01183) we showed that the reliability of the estimates is not affected by the size of the potential. However, we have to keep in mind that we need a part of the box where the ligand can stay far away from the protein in order to represent the unbound state in a satisfactory way. An important point here is that the sphere constraint __must__ be inside the box, otherwise the entropic correction will not accurately account for the loss of configurational space. To visually inspect how large the potential is, and to get a feel for the possible movements of the ligand, we can visualize both the system and the restraint with VMD (downloadable [here](https://www.ks.uiuc.edu/Research/vmd/)).

We can open VMD and load the `starting.gro` structure file in the GitHub folder. After the structure is loaded and the visualization has been set up at your taste, you can open the Tk console and write
```
pbc box
```
which draws the cubic box in which the system is inserted. 

Now we can generate the atom selection we defined after checking the RMSF:
```
set sel [atomselect top "resid 71 to 134 or resid 138 to 139 or resid 142 to 162"]
```
The console should answer
```
atomselectXX
```
Where `XX` is a number. The selection for the reference frame has been defined and named `$sel`. Now we can compute the position of the center of mass of this set of atoms:
```
measure center $sel weight mass
```
The console should answer with the position of the center of mass (in Ångstrom)
```
35.28876876831055 34.06196594238281 32.041622161865234
```
Knowing this information, we can draw the sphere with a radius of (for example) 2 nm  (20 Å) with the following command
```
draw sphere {35.28876876831055 34.06196594238281 32.041622161865234} radius 20 resolution 100
```
Receiving a number as an answer from the console. Such number is the ID of the 3D object we just draw. The drawn sphere is opaque, not allowing us to see inside it; to make it transparent, we need to specify that we want a transparent material
```
draw material Transparent
```
We can see that the sphere contains the entire domain, but it is probably too small to represent the unbound state in a precise way. Let's delete the sphere using the ID of the 3D object (let's say that it is `14`), and plot a new sphere of radius 2.8 nm
```
draw delete 14
draw sphere {35.28876876831055 34.06196594238281 32.041622161865234} radius 28 resolution 100
```

You can see the expected result below
<p align="center">
  <img src="https://github.com/riccardocapelli/VMetaD-tutorial/blob/main/img/sphere_box.jpg?raw=true" alt="Alt text" width="50%">
  <br>
  <em>Cartoon representation of the lysozyme-benzene complex, including the restraining potential applied within a 2.8 nm radius of the reference frame center of mass. The boundaries of the simulation box are also highlighted to show that the sphere is entirely contained by the box. </em>
</p>

## The PLUMED input file
_(You can read the following line-by-line description keeping an eye on the `plumed.dat` file in the GitHub folder as a reference)_

### Reference frame setup
We start with the `WHOLEMOLECULES` instruction, to be sure that lysozyme (`ENTITY0`) will not be broken by the periodic boundary condition, as well as the benzene molecule (`ENTITY1`):
```plumed
WHOLEMOLECULES ENTITY0=1-1284 ENTITY1=1285-1290
```
Now that we are sure of the integrity of the structures in PLUMED, we perform the rototranslational fit of the system to make sure that the protein will be in the fixed reference frame position:
```plumed
#HIDDEN
WHOLEMOLECULES ENTITY0=1-1284 ENTITY1=1285-1290
#ENDHIDDEN
FIT_TO_TEMPLATE REFERENCE=ref_ca.pdb TYPE=OPTIMAL
```

We then start with the groups definition. We previously prepared a GROMACS index file (`index.ndx`) with all the groups named as intended. As an alternative, you can also define such groups with atom ids.
```plumed
#HIDDEN
WHOLEMOLECULES ENTITY0=1-1284 ENTITY1=1285-1290

FIT_TO_TEMPLATE REFERENCE=ref_ca.pdb TYPE=OPTIMAL
#ENDHIDDEN
prot_noh: GROUP NDX_FILE=index.ndx NDX_GROUP=Protein-H
sph: GROUP NDX_FILE=index.ndx NDX_GROUP=sphere
lig: GROUP NDX_FILE=index.ndx NDX_GROUP=ligand
```
We have three groups: 
* `prot_noh`, which contains all the non-hydrogen atoms of the protein (for our multi-eGO potential is equivalent to all the protein, but in all-atom representation it makes a difference);
* `sph`, which contains the atoms that define the reference frame (the C-alpha of the selected residues -see above-);
* `lig`, which contains the atoms of the ligand.

After the definition of the groups, to avoid that the passage in a periodic boundary conditions causes a "jump" of the ligand with respect to the protein, we add a `WRAPAROUND` instruction:
```plumed
#HIDDEN
WHOLEMOLECULES ENTITY0=1-1284 ENTITY1=1285-1290

FIT_TO_TEMPLATE REFERENCE=ref_ca.pdb TYPE=OPTIMAL

prot_noh: GROUP NDX_FILE=index.ndx NDX_GROUP=Protein-H
sph: GROUP NDX_FILE=index.ndx NDX_GROUP=sphere
lig: GROUP NDX_FILE=index.ndx NDX_GROUP=ligand
#ENDHIDDEN
WRAPAROUND ATOMS=lig AROUND=sph
```
Ending the fitting part of the PLUMED input.

### Spherical coordinates definition
We now compute the position of the center of mass of the atoms defining the reference frame ($(x,y,z)=(0,0,0)$ in our CV space), and the center of mass of the ligand:
```plumed
#HIDDEN
WHOLEMOLECULES ENTITY0=1-1284 ENTITY1=1285-1290

FIT_TO_TEMPLATE REFERENCE=ref_ca.pdb TYPE=OPTIMAL

prot_noh: GROUP NDX_FILE=index.ndx NDX_GROUP=Protein-H
sph: GROUP NDX_FILE=index.ndx NDX_GROUP=sphere
lig: GROUP NDX_FILE=index.ndx NDX_GROUP=ligand

WRAPAROUND ATOMS=lig AROUND=sph
#ENDHIDDEN
sph_center: COM ATOMS=sph
lig_center: COM ATOMS=lig

sph_coord: POSITION ATOM=sph_center NOPBC
lig_coord: POSITION ATOM=lig_center NOPBC
```
From the position, we can obtain the cartesian coordinates of the ligand in this reference frame
```plumed
#HIDDEN
WHOLEMOLECULES ENTITY0=1-1284 ENTITY1=1285-1290

FIT_TO_TEMPLATE REFERENCE=ref_ca.pdb TYPE=OPTIMAL

prot_noh: GROUP NDX_FILE=index.ndx NDX_GROUP=Protein-H
sph: GROUP NDX_FILE=index.ndx NDX_GROUP=sphere
lig: GROUP NDX_FILE=index.ndx NDX_GROUP=ligand

WRAPAROUND ATOMS=lig AROUND=sph

sph_center: COM ATOMS=sph
lig_center: COM ATOMS=lig

sph_coord: POSITION ATOM=sph_center NOPBC
lig_coord: POSITION ATOM=lig_center NOPBC
#ENDHIDDEN
abs_x: MATHEVAL ARG=lig_coord.x,sph_coord.x FUNC=x-y PERIODIC=NO
abs_y: MATHEVAL ARG=lig_coord.y,sph_coord.y FUNC=x-y PERIODIC=NO
abs_z: MATHEVAL ARG=lig_coord.z,sph_coord.z FUNC=x-y PERIODIC=NO
```
and via the usual transformation, obtain the final spherical coordinates $(\rho,\theta,\varphi)$
```plumed
#HIDDEN
WHOLEMOLECULES ENTITY0=1-1284 ENTITY1=1285-1290

FIT_TO_TEMPLATE REFERENCE=ref_ca.pdb TYPE=OPTIMAL

prot_noh: GROUP NDX_FILE=index.ndx NDX_GROUP=Protein-H
sph: GROUP NDX_FILE=index.ndx NDX_GROUP=sphere
lig: GROUP NDX_FILE=index.ndx NDX_GROUP=ligand

WRAPAROUND ATOMS=lig AROUND=sph

sph_center: COM ATOMS=sph
lig_center: COM ATOMS=lig

sph_coord: POSITION ATOM=sph_center NOPBC
lig_coord: POSITION ATOM=lig_center NOPBC

abs_x: MATHEVAL ARG=lig_coord.x,sph_coord.x FUNC=x-y PERIODIC=NO
abs_y: MATHEVAL ARG=lig_coord.y,sph_coord.y FUNC=x-y PERIODIC=NO
abs_z: MATHEVAL ARG=lig_coord.z,sph_coord.z FUNC=x-y PERIODIC=NO
#ENDHIDDEN
rho: MATHEVAL ARG=abs_x,abs_y,abs_z FUNC=sqrt(x*x+y*y+z*z) PERIODIC=NO
theta: MATHEVAL ARG=abs_z,rho FUNC=acos(x/y) PERIODIC=0.,pi
phi: MATHEVAL ARG=abs_x,abs_y FUNC=atan2(y,x) PERIODIC=-pi,pi
```
which will be our CVs.

### Restraining
We now have to impose the spherical restraint. We put a `UPPER_WALLS` which impedes the ligand to go farther than 2.8 nm:
```plumed
#HIDDEN
WHOLEMOLECULES ENTITY0=1-1284 ENTITY1=1285-1290

FIT_TO_TEMPLATE REFERENCE=ref_ca.pdb TYPE=OPTIMAL

prot_noh: GROUP NDX_FILE=index.ndx NDX_GROUP=Protein-H
sph: GROUP NDX_FILE=index.ndx NDX_GROUP=sphere
lig: GROUP NDX_FILE=index.ndx NDX_GROUP=ligand

WRAPAROUND ATOMS=lig AROUND=sph

sph_center: COM ATOMS=sph
lig_center: COM ATOMS=lig

sph_coord: POSITION ATOM=sph_center NOPBC
lig_coord: POSITION ATOM=lig_center NOPBC

abs_x: MATHEVAL ARG=lig_coord.x,sph_coord.x FUNC=x-y PERIODIC=NO
abs_y: MATHEVAL ARG=lig_coord.y,sph_coord.y FUNC=x-y PERIODIC=NO
abs_z: MATHEVAL ARG=lig_coord.z,sph_coord.z FUNC=x-y PERIODIC=NO

rho: MATHEVAL ARG=abs_x,abs_y,abs_z FUNC=sqrt(x*x+y*y+z*z) PERIODIC=NO
theta: MATHEVAL ARG=abs_z,rho FUNC=acos(x/y) PERIODIC=0.,pi
phi: MATHEVAL ARG=abs_x,abs_y FUNC=atan2(y,x) PERIODIC=-pi,pi
#ENDHIDDEN
restr: UPPER_WALLS ARG=rho AT=2.8 KAPPA=10000
```
the $k$ value is 10,000 kJ/mol/nm^2, which means that if the ligand is out by 0.1 nm he will feel a potential of 100 kJ/mol.

One effect that we should take into account is the possibility that the ligand, in advanced phases of the simulation, will try to unfold the protein, being the place occupied by it the volume portion with less history-dependent potential deposited. To limit this phenomenon we will put in place a RMSD restraining that will be removed during the reweighting procedure. To compute the RMSD we will use the same atoms included in the `ref_ca.pdb` file instruction
```plumed
#HIDDEN
WHOLEMOLECULES ENTITY0=1-1284 ENTITY1=1285-1290

FIT_TO_TEMPLATE REFERENCE=ref_ca.pdb TYPE=OPTIMAL

prot_noh: GROUP NDX_FILE=index.ndx NDX_GROUP=Protein-H
sph: GROUP NDX_FILE=index.ndx NDX_GROUP=sphere
lig: GROUP NDX_FILE=index.ndx NDX_GROUP=ligand

WRAPAROUND ATOMS=lig AROUND=sph

sph_center: COM ATOMS=sph
lig_center: COM ATOMS=lig

sph_coord: POSITION ATOM=sph_center NOPBC
lig_coord: POSITION ATOM=lig_center NOPBC

abs_x: MATHEVAL ARG=lig_coord.x,sph_coord.x FUNC=x-y PERIODIC=NO
abs_y: MATHEVAL ARG=lig_coord.y,sph_coord.y FUNC=x-y PERIODIC=NO
abs_z: MATHEVAL ARG=lig_coord.z,sph_coord.z FUNC=x-y PERIODIC=NO

rho: MATHEVAL ARG=abs_x,abs_y,abs_z FUNC=sqrt(x*x+y*y+z*z) PERIODIC=NO
theta: MATHEVAL ARG=abs_z,rho FUNC=acos(x/y) PERIODIC=0.,pi
phi: MATHEVAL ARG=abs_x,abs_y FUNC=atan2(y,x) PERIODIC=-pi,pi

restr: UPPER_WALLS ARG=rho AT=2.8 KAPPA=10000
#ENDHIDDEN
rmsd: RMSD REFERENCE=ref_ca.pdb TYPE=OPTIMAL
restr_rmsd: RESTRAINT ARG=rmsd AT=0. KAPPA=250.0
```

### Reweighting CVs
As anticipated in the theory part, to compute binding free energy differences we will need to reweight our free energy landscape on new apt CVs which will allow us in discriminating efficiently the bound and unbound states. Like in the original work, here we will use again the distance from the origin of the reference frame $\rho$ and the number of contacts between the protein and the ligand. This will guarantee that if the guest can be considered not in contact with the host (and thus defining the unbound state), even if $\rho$ is large.

The total number of contact $c$ is defined with a switching function

$$
c = \sum_{i,j} \frac{ 1 - \left(\frac{r_{ij}}{r_0}\right)^{6} } 
{1 - \left(\frac{r_{ij}}{r_0}\right)^{12} }
$$

which runs on all the (heavy) atoms of the protein and all the atoms of the ligand. This is implemented with `COOORDINATION`:
```plumed
#HIDDEN
WHOLEMOLECULES ENTITY0=1-1284 ENTITY1=1285-1290

FIT_TO_TEMPLATE REFERENCE=ref_ca.pdb TYPE=OPTIMAL

prot_noh: GROUP NDX_FILE=index.ndx NDX_GROUP=Protein-H
sph: GROUP NDX_FILE=index.ndx NDX_GROUP=sphere
lig: GROUP NDX_FILE=index.ndx NDX_GROUP=ligand

WRAPAROUND ATOMS=lig AROUND=sph

sph_center: COM ATOMS=sph
lig_center: COM ATOMS=lig

sph_coord: POSITION ATOM=sph_center NOPBC
lig_coord: POSITION ATOM=lig_center NOPBC

abs_x: MATHEVAL ARG=lig_coord.x,sph_coord.x FUNC=x-y PERIODIC=NO
abs_y: MATHEVAL ARG=lig_coord.y,sph_coord.y FUNC=x-y PERIODIC=NO
abs_z: MATHEVAL ARG=lig_coord.z,sph_coord.z FUNC=x-y PERIODIC=NO

rho: MATHEVAL ARG=abs_x,abs_y,abs_z FUNC=sqrt(x*x+y*y+z*z) PERIODIC=NO
theta: MATHEVAL ARG=abs_z,rho FUNC=acos(x/y) PERIODIC=0.,pi
phi: MATHEVAL ARG=abs_x,abs_y FUNC=atan2(y,x) PERIODIC=-pi,pi

restr: UPPER_WALLS ARG=rho AT=2.8 KAPPA=10000

rmsd: RMSD REFERENCE=ref_ca.pdb TYPE=OPTIMAL
restr_rmsd: RESTRAINT ARG=rmsd AT=0. KAPPA=250.0
#ENDHIDDEN
c: COORDINATION GROUPA=lig GROUPB=prot_noh R_0=0.45
```
where we set a $r_0$ parameter at 0.45 nm.

### Volume-based Metadynamics
We now set up the VMetaD:
```plumed
#HIDDEN
WHOLEMOLECULES ENTITY0=1-1284 ENTITY1=1285-1290

FIT_TO_TEMPLATE REFERENCE=ref_ca.pdb TYPE=OPTIMAL

prot_noh: GROUP NDX_FILE=index.ndx NDX_GROUP=Protein-H
sph: GROUP NDX_FILE=index.ndx NDX_GROUP=sphere
lig: GROUP NDX_FILE=index.ndx NDX_GROUP=ligand

WRAPAROUND ATOMS=lig AROUND=sph

sph_center: COM ATOMS=sph
lig_center: COM ATOMS=lig

sph_coord: POSITION ATOM=sph_center NOPBC
lig_coord: POSITION ATOM=lig_center NOPBC

abs_x: MATHEVAL ARG=lig_coord.x,sph_coord.x FUNC=x-y PERIODIC=NO
abs_y: MATHEVAL ARG=lig_coord.y,sph_coord.y FUNC=x-y PERIODIC=NO
abs_z: MATHEVAL ARG=lig_coord.z,sph_coord.z FUNC=x-y PERIODIC=NO

rho: MATHEVAL ARG=abs_x,abs_y,abs_z FUNC=sqrt(x*x+y*y+z*z) PERIODIC=NO
theta: MATHEVAL ARG=abs_z,rho FUNC=acos(x/y) PERIODIC=0.,pi
phi: MATHEVAL ARG=abs_x,abs_y FUNC=atan2(y,x) PERIODIC=-pi,pi

restr: UPPER_WALLS ARG=rho AT=2.8 KAPPA=10000

rmsd: RMSD REFERENCE=ref_ca.pdb TYPE=OPTIMAL
restr_rmsd: RESTRAINT ARG=rmsd AT=0. KAPPA=250.0

c: COORDINATION GROUPA=lig GROUPB=prot_noh R_0=0.45
#ENDHIDDEN
METAD ...
  ARG=rho,theta,phi
  GRID_MIN=0,0.,-pi
  GRID_MAX=3.5,pi,pi
  SIGMA=0.1,pi/16.,pi/8
  HEIGHT=0.5
  PACE=2000
  BIASFACTOR=10
  TEMP=300
  LABEL=metad
  CALC_RCT
... METAD
```
The most exotic option used is `CALC_RCT`, which allows the calculation on the fly of the [Tiwary-Parrinello estimator](https://doi.org/10.1021/jp504920s) that we will use for reweighting.

### Printing
We finally print all the relevant files that we will use for post-processing and analysis. 
```plumed
#HIDDEN
WHOLEMOLECULES ENTITY0=1-1284 ENTITY1=1285-1290

FIT_TO_TEMPLATE REFERENCE=ref_ca.pdb TYPE=OPTIMAL

prot_noh: GROUP NDX_FILE=index.ndx NDX_GROUP=Protein-H
sph: GROUP NDX_FILE=index.ndx NDX_GROUP=sphere
lig: GROUP NDX_FILE=index.ndx NDX_GROUP=ligand

WRAPAROUND ATOMS=lig AROUND=sph

sph_center: COM ATOMS=sph
lig_center: COM ATOMS=lig

sph_coord: POSITION ATOM=sph_center NOPBC
lig_coord: POSITION ATOM=lig_center NOPBC

abs_x: MATHEVAL ARG=lig_coord.x,sph_coord.x FUNC=x-y PERIODIC=NO
abs_y: MATHEVAL ARG=lig_coord.y,sph_coord.y FUNC=x-y PERIODIC=NO
abs_z: MATHEVAL ARG=lig_coord.z,sph_coord.z FUNC=x-y PERIODIC=NO

rho: MATHEVAL ARG=abs_x,abs_y,abs_z FUNC=sqrt(x*x+y*y+z*z) PERIODIC=NO
theta: MATHEVAL ARG=abs_z,rho FUNC=acos(x/y) PERIODIC=0.,pi
phi: MATHEVAL ARG=abs_x,abs_y FUNC=atan2(y,x) PERIODIC=-pi,pi

restr: UPPER_WALLS ARG=rho AT=2.8 KAPPA=10000

rmsd: RMSD REFERENCE=ref_ca.pdb TYPE=OPTIMAL
restr_rmsd: RESTRAINT ARG=rmsd AT=0. KAPPA=250.0

c: COORDINATION GROUPA=lig GROUPB=prot_noh R_0=0.45

METAD ...
  ARG=rho,theta,phi
  GRID_MIN=0,0.,-pi
  GRID_MAX=3.5,pi,pi
  SIGMA=0.1,pi/16.,pi/8
  HEIGHT=0.5
  PACE=2000
  BIASFACTOR=10
  TEMP=300
  LABEL=metad
  CALC_RCT
... METAD
#ENDHIDDEN
PRINT ARG=metad.* FILE=metad_data.dat STRIDE=200
PRINT ARG=rmsd,restr_rmsd.bias FILE=rmsd_restraint.dat STRIDE=200
PRINT ARG=restr.bias FILE=sphere_restraint.dat STRIDE=200
PRINT ARG=abs_x,abs_y,abs_z FILE=xyz_coord.dat STRIDE=200
PRINT ARG=rho,theta,phi FILE=rtp_coord.dat STRIDE=200
PRINT ARG=c,rho FILE=coord_rho.dat STRIDE=200

FLUSH STRIDE=200
```
We will have all the VMetaD quantities in `metad_data.dat`, the restraints data in `{rmsd,sphere}_restaint.dat`, the $(x,y,z)$ and $(\rho,\theta,\varphi)$ coordinates in `xyz_coord.dat` and `rtp_coord.dat`, respectively, and the reweighting CVs in `coord_rho.dat`.

___PLEASE NOTE___: all the printing frequencies are synchronized with the VMetaD `PACE` (every 10 print we deposit 1 gaussian), and it should be synchronized also with trajectory printing (I personally suggest the same frequency used for gaussian deposition). This allows us to restart safely in case of issues and re-run or re-analyze with new CVs the run, if needed.

### Last advices before launching the simulation
One issue that can be observed when launching when running VMetaD is the sudden interruption of the simulation with a cryptic error. Please check the position of the ligand: in most cases, the system reached $(\rho,\theta,\varphi)=(0,0,0)$, where the derivatives of the potential cannot be defined, and thus PLUMED sends an error. Despite being annoying, this is perfectly normal, and does not invalidate the run. Please restart it from the previous checkpoint (save a checkpoint often!).

Another consideration regarding the singularity of the reference frame is its position with respect to the actual binding site. In close proximity of it (let's say less than 2 sigmas in $\rho$) even an extremely small movement in any direction makes the ligand move of several sigmas in $\theta$ and $\varphi$, with the possibility of underestimation of the bias. This could imply the need of longer simulation times to fill up the basin. To limit this effect, we suggest to verify if the origin of the reference frame is farther that 2-3 sigmas in $\rho$ from the binding site (the perfect situation would be setting the origin in a point occupied by the host, at 4-5 Å from the binding site).

### Launch!
You can now run the VMetaD tutorial. We advice you to run it for (at least) 500 ns. In this example, we run it for 1 µs.

##### [Back to VMetad home](NAVIGATION.md)
