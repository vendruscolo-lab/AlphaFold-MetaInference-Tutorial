__Note:__ The following code can be run in google colab or jupyter notebook. [This](https://colab.research.google.com/github/vendruscolo-lab/AlphaFold-MetaInference-Tutorial/blob/main/AF-IDP_colab.ipynb) google colab provides with the full AF-MI protocol, i.e software installation, CALVADOS, PLUMED input preparation, running AF-MI, and analysis. Lets break it down and explain each part. 

## Software installation



The code below needs to run in colab. We first install condacolab which is essential in running conda google colab. This step is not necessary if one runs just in a jupyter-notebook locally. 


```python
#Install collab
!pip install -q condacolab
import condacolab
condacolab.install()
```
The code below installs various necessary packages such as PLUMED, mpi, OPENMM, pandas, gromacs, biopython, pulchra. Moreover it builds the OPENMM-PLUMED-MPI patch which is necessary for running metainference in OPENMM.  

```python
import os
home=os.getcwd()
os.chdir(home)

#Download OpenMM-Plumed-MPI patch
!conda install -y -c conda-forge mpich mpi4py openmm=8.0 plumed=2.8.2=mpi_mpich_h7ded119_0 py-plumed cmake swig pandas mdtraj biopython matplotlib gromacs
!conda install -y -c anaconda ipykernel
!conda install -y -c bioconda pulchra

!git clone https://github.com/vendruscolo-lab/OpenMM-Plumed-MPI

os.chdir(home+'/OpenMM-Plumed-MPI')

#Build openmm-plumed-mpi
!mkdir build install openmm -p plumed/include -p plumed/lib
!unzip openmm.zip -d openmm
!unzip plumed_lib.zip -d plumed/lib
!unzip plumed_include.zip -d plumed/include
os.chdir(os.getcwd()+'/build')
!cmake ..
!make
!make install
!make PythonInstall
os.chdir(os.getcwd()+'/../install/lib')
!cp -r * /usr/local/lib/
```

## Run Alpha-Fold distance map prediction


__For this example (TDP-43 WtoA):__

The AF distance map has already been calculated and it is loaded below. We here show the AF prediction for the a) Inter-residue distances, b)PLDDT score and c)PAE per residue pair. 

<p align="center">
  <img src="https://github.com/vendruscolo-lab/AlphaFold-MetaInference-Tutorial/blob/main/images/alphafold2_ptm_model_3_seed_000_distmat.png?raw=true" alt="Alt text" width="50%">
    <img src="https://github.com/vendruscolo-lab/AlphaFold-MetaInference-Tutorial/blob/main/images/tdp43_WtoA_bf4cc_plddt.png?raw=true" alt="Alt text" width="50%">
    <img src="https://github.com/vendruscolo-lab/AlphaFold-MetaInference-Tutorial/blob/main/images/tdp43_WtoA_bf4cc_pae.png?raw=true" alt="Alt text" width="50%">
  <br> a) Inter-residue distances, b)PLDDT score and c)PAE per residue pair. 
  <em> </em>
</p>



__For arbitrary protein sequences:__

*   Open [this](https://github.com/zshengyu14/ColabFold_distmats/blob/main/AlphaFold2.ipynb) link and chose colab.
*   Input the protein sequence  as query sequence.
*   The rest of the options remain default and cells are run until the end.
*   Download the link with the AF data and upload it as AF_data in AlphaFold-IDP folder

## Setup protein system in CALVADOS2 and OPENMM

```python
os.chdir(home)
!git clone https://github.com/vendruscolo-lab/AlphaFold-IDP
os.chdir(home+'/AlphaFold-IDP/prep_run')
```
In the following step we need to define ```python fasta_sequence, pH, temp, ionic, PAE_cut, NR, ordered_domains, disordered_domains, pdb_af,json_af,npy_af, mean_af```. These variables respectively stand for the protein sequence to be simulated in AF-MI, the simulation pH, the simulation temperature (in K), the ionic strength of the solution, the highest AF predicted alighment error the considered inter-residue distances will have (residue distances with higher PAE are not considered in restraints), the number of replicas, the regions of the ordered domains (regions of more than 3 residues with PLDDT>0.75) where RMSD walls wil be used, the disordered regions (usually regions of more than 3 residues with PLDDT<0.75), AF predicted pdb file, AF predicted json containing the PAE per residue pair file, per residue pair probability distribution AF prediction file, mean AF inter-residue distance map prediction file  

```python
import shutil
import csv
dir=os.getcwd()
###################### The entries below need to be adapted in each simulation ######################
#TDP-43 sequence
fasta_sequence="MSEYIRVTEDENDEPIEIPSEDDGTVLLSTVTAQFPGACGLRYRNPVSQCMRGVRLVEGILHAPDAGAGNLVYVVNYPKDNKRKMDETDASSAVKVKRAVQKTSDLIVLGLPAKTTEQDLKEYFSTFGEVLMVQVKKDLKTGHSKGFGFVRFTEYETQVKVMSQRHMIDGRACDCKLPNSKQSQDEPLRSRKVFVGRCTEDMTEDELREFFSQYGDVMDVFIPKPFRAFAFVTFADDQIAQSLCGEDLIIKGISVHISNAEPKHNSNRQLERSGRFGGNPGGFGNQGGFGNSRGGGAGLGNNQGSNMGGGMNFGAFSINPAMMAAAQAALQSSAGMMGMLASQQNQSGPSGNNQNQGNMQREPNQAFGSGNNSYSGSNSGAAIGAGSASNAGSGSGFNGGFGSSMDSKSSGAGM"
#Conditions
pH=7.4
temp=298
ionic=0.2
PAE_cut=5
NR=2
#Decide the plddt based ordered (od) and disordered (dd) regions
ordered_domains = {'od1': [3, 79], 'od2': [104,178],'od3':[191,260]}
disordered_domains={'dd1': [1,2],'dd2':[80,103],'dd3':[179,190],'dd4':[261,414]}

#AF input created in AF-distance map prediction and used in AF-MI
pdb_af='tdp43_WtoA_bf4cc_unrelaxed_rank_001_alphafold2_ptm_model_3_seed_000.pdb'
json_af='tdp43_WtoA_bf4cc_predicted_aligned_error_v1.json'
npy_af='alphafold2_ptm_model_3_seed_000_prob_distributions.npy'
mean_af='alphafold2_ptm_model_3_seed_000_mean.csv'
#Copy AF data
shutil.copy2(dir+"/../AF_DATA/"+pdb_af, dir+"/pdb_af.pdb")
shutil.copy2(dir+"/../AF_DATA/"+json_af, dir+"/pae.json")
shutil.copy2(dir+"/../AF_DATA/tdp43_WtoA_bf4cc_distmat/"+mean_af, dir+"/mean_af.csv")
shutil.copy2(dir+"/../AF_DATA/tdp43_WtoA_bf4cc_distmat/"+npy_af, dir+"/prob.npy")
####################################################################################################
```

In the following step we create the OPENMM CALVADOS2 .xml files

```python 
f = open("sequence.dat", "w")
f.write(fasta_sequence)
f.close()
#Write the csv files
with open('ordered_domains.csv', 'w') as f:  # You will need 'wb' mode in Python 2.x
    w = csv.DictWriter(f, ordered_domains.keys())
    w.writeheader()
    w.writerow(ordered_domains)
with open('disordered_domains.csv', 'w') as f:  # You will need 'wb' mode in Python 2.x
    w = csv.DictWriter(f, disordered_domains.keys())
    w.writeheader()
    w.writerow(disordered_domains)

#Copy OPENMM calvados forcefield files
shutil.copy2(dir+"/../scripts_prep/gen_xml_and_constraints.py", dir)
shutil.copy2(dir+"/../scripts_prep/residues.csv", dir)

path_gen_xml = dir+'/gen_xml_and_constraints.py sequence.dat '+str(pH)+' '+str(temp)+' '+str(ionic)
print(path_gen_xml)
os.system(f'python {path_gen_xml}')
```
Make the AF-MI PLUMED ``` plumed.dat``` input file 
```python
#Make plumed files.
#Copy and run the prep script that makes the plumed file.
#The Collective variables (CVs) in these case are chosen to be the torsion angles between structured domains.
import subprocess
shutil.copy2(dir+"/../scripts_prep/make_plumed_distmat.py", dir)
subprocess.run(['python', str(dir)+'/make_plumed_distmat.py', 'sequence.dat',str(PAE_cut), '0.2'], capture_output=True, text=True)
```
The plumed.dat file can be seen below. The ```distance_rest_domains``` are the af-distances do be restrained. For TDP-43 WtoA, the three ordered domains are RMSD restrained with RMSD1,RMSD2,RMSD3. 
Note the __FILL__ entries the user should specify for the specific system at hand. This includes the CV definition, biasfactor, sigma, number of bins in the grid for saving the FES along the PB MetaD simulation. More info on PB-MetaD can be found [here](https://www.plumed.org/doc-v2.9/user-doc/html/_p_b_m_e_t_a_d.html). A usual rule of thumb for the ```biasfactor``` value is ```10*sqrt(number of biased CVs)```. 
For TDP-43 WtoA, the ```plumed.dat``` file can be found [here](https://github.com/vendruscolo-lab/AlphaFold-IDP/blob/main/scripts_prep/plumed_TDP-43.dat)




```plumed
#SOLUTIONFILE=plumed_TDP-43.dat
MOLINFO MOLTYPE=protein STRUCTURE=input_af.pdb
WHOLEMOLECULES ENTITY0=1-414
    
distance_rest_domains:  CONTACTMAP ...
ATOMS1=1,4
ATOMS2=2,5
ATOMS3=2,20
ATOMS4=28,79
ATOMS5=31,79
ATOMS6=35,79
ATOMS7=37,79
ATOMS8=38,79
ATOMS9=39,79
ATOMS10=40,79
ATOMS11=41,79
ATOMS12=42,79
ATOMS13=74,79
ATOMS14=75,79
ATOMS15=76,79
ATOMS16=102,105
ATOMS17=102,106
ATOMS18=102,150
ATOMS19=102,151
ATOMS20=102,152
ATOMS21=102,155
ATOMS22=102,158
ATOMS23=102,159
ATOMS24=102,161
ATOMS25=103,106
ATOMS26=103,107
ATOMS27=103,108
ATOMS28=103,120
ATOMS29=103,121
ATOMS30=103,124
ATOMS31=103,130
ATOMS32=103,131
ATOMS33=103,132
ATOMS34=103,133
ATOMS35=103,134
ATOMS36=103,135
ATOMS37=103,147
ATOMS38=103,148
ATOMS39=103,149
ATOMS40=103,150
ATOMS41=103,151
ATOMS42=103,152
ATOMS43=103,153
ATOMS44=103,154
ATOMS45=103,155
ATOMS46=103,156
ATOMS47=103,157
ATOMS48=103,158
ATOMS49=103,159
ATOMS50=103,160
ATOMS51=103,161
ATOMS52=103,162
ATOMS53=103,175
ATOMS54=103,176
ATOMS55=103,177
ATOMS56=103,178
ATOMS57=104,178
ATOMS58=105,178
ATOMS59=105,179
ATOMS60=106,178
ATOMS61=106,179
ATOMS62=107,178
ATOMS63=107,179
ATOMS64=108,178
ATOMS65=108,179
ATOMS66=109,178
ATOMS67=109,179
ATOMS68=110,178
ATOMS69=111,178
ATOMS70=120,178
ATOMS71=124,178
ATOMS72=131,178
ATOMS73=132,178
ATOMS74=133,178
ATOMS75=134,178
ATOMS76=135,178
ATOMS77=136,178
ATOMS78=137,178
ATOMS79=145,178
ATOMS80=146,178
ATOMS81=147,178
ATOMS82=147,179
ATOMS83=148,178
ATOMS84=148,179
ATOMS85=149,178
ATOMS86=149,179
ATOMS87=150,178
ATOMS88=150,179
ATOMS89=151,178
ATOMS90=151,179
ATOMS91=152,178
ATOMS92=155,178
ATOMS93=158,178
ATOMS94=159,178
ATOMS95=159,179
ATOMS96=161,178
ATOMS97=161,179
ATOMS98=162,178
ATOMS99=162,179
ATOMS100=163,178
ATOMS101=164,178
ATOMS102=165,178
ATOMS103=166,178
ATOMS104=173,178
ATOMS105=174,178
ATOMS106=174,179
ATOMS107=175,178
ATOMS108=175,179
ATOMS109=176,179
ATOMS110=189,192
ATOMS111=189,233
ATOMS112=189,240
ATOMS113=190,193
ATOMS114=190,194
ATOMS115=190,195
ATOMS116=190,209
ATOMS117=190,210
ATOMS118=190,211
ATOMS119=190,212
ATOMS120=190,213
ATOMS121=190,214
ATOMS122=190,216
ATOMS123=190,217
ATOMS124=190,218
ATOMS125=190,219
ATOMS126=190,220
ATOMS127=190,230
ATOMS128=190,231
ATOMS129=190,232
ATOMS130=190,233
ATOMS131=190,234
ATOMS132=190,235
ATOMS133=190,236
ATOMS134=190,237
ATOMS135=190,238
ATOMS136=190,239
ATOMS137=190,240
ATOMS138=190,241
ATOMS139=190,242
ATOMS140=190,243
ATOMS141=190,244
ATOMS142=190,257
ATOMS143=190,258
ATOMS144=191,260
ATOMS145=192,260
ATOMS146=192,261
ATOMS147=193,260
ATOMS148=193,261
ATOMS149=193,262
ATOMS150=194,260
ATOMS151=194,261
ATOMS152=194,262
ATOMS153=195,260
ATOMS154=195,261
ATOMS155=196,260
ATOMS156=197,260
ATOMS157=211,260
ATOMS158=219,260
ATOMS159=221,260
ATOMS160=222,260
ATOMS161=223,260
ATOMS162=229,260
ATOMS163=229,261
ATOMS164=230,260
ATOMS165=230,261
ATOMS166=231,260
ATOMS167=231,261
ATOMS168=232,260
ATOMS169=232,261
ATOMS170=233,260
ATOMS171=234,260
ATOMS172=234,261
ATOMS173=236,260
ATOMS174=237,260
ATOMS175=238,260
ATOMS176=239,260
ATOMS177=240,260
ATOMS178=241,260
ATOMS179=241,261
ATOMS180=242,260
ATOMS181=243,260
ATOMS182=243,261
ATOMS183=244,260
ATOMS184=245,260
ATOMS185=246,260
ATOMS186=247,260
ATOMS187=248,260
ATOMS188=255,260
ATOMS189=256,260
ATOMS190=256,261
ATOMS191=257,260
ATOMS192=257,261
ATOMS193=258,261
ATOMS194=258,262
ATOMS195=259,262
ATOMS196=298,301
ATOMS197=300,303
ATOMS198=304,307
ATOMS199=318,321
ATOMS200=321,324
ATOMS201=321,325
ATOMS202=321,326
ATOMS203=322,325
ATOMS204=322,326
ATOMS205=323,326
ATOMS206=323,327
ATOMS207=324,327
ATOMS208=325,328
ATOMS209=326,329
ATOMS210=348,351
ATOMS211=350,353
ATOMS212=351,354
ATOMS213=352,355
ATOMS214=357,360
ATOMS215=360,363
ATOMS216=371,374
ATOMS217=376,379
ATOMS218=386,389
ATOMS219=392,395
ATOMS220=394,397
ATOMS221=396,399
ATOMS222=399,402
ATOMS223=400,403
ATOMS224=402,405
ATOMS225=402,406
ATOMS226=407,411
ATOMS227=408,411
SWITCH={CUSTOM FUNC=x R_0=1}
...

af_dist_rest: CONSTANT VALUES=1.0032847926020623,1.0311901761218905,0.7455537494271994,1.7088723132386805,1.713635583035648,1.5095473634079102,1.1722138326615095,1.2549746362492442,0.9444370301440359,1.2123430594801905,1.5226123476400972,1.3369803665205837,1.4350552991032601,1.271596952714026,0.8981292195618154,1.084396699629724,1.2897209025919436,1.5121250515803697,1.0559752460569143,1.0836912801489234,0.707867799885571,0.9799594385549426,1.2259459158405663,1.5500808618962765,0.9052655544131994,1.3762398231774569,1.5152623694390062,1.68192987870425,1.6158358810469509,1.3273277390748264,1.1055688032880424,0.7953455133363603,1.0522321155294776,1.3843332368880512,1.5040866650640965,1.9200282409787175,1.9061048422008757,1.6269824001938105,1.2354281768202782,1.0676583984866739,0.6573523612692953,0.6149551127105952,0.7450687747448683,0.7618590366095304,0.5613296514376999,1.0004315715283156,0.8841317959129811,0.5908093221485615,0.9980626434087752,1.263882938399911,1.1281702172011137,1.1415027465671301,1.3193011650815607,1.4435925820842384,1.1070262966677549,1.31037020906806,0.8523610839620233,0.6740215603262186,1.1396584818139672,0.7467208079993726,1.1898928672075273,0.4552792670205237,0.8395836766809226,1.017876618169248,1.3974865483120085,1.0051763331517576,1.1795675404369834,1.4049736192449929,1.459107712842524,1.5650913815945389,1.59471937045455,1.5574198851361871,1.3527641544118525,1.4486837401986123,1.1706370314583183,1.4376086220145226,1.1653784492984414,1.7256552413105966,1.4055226838216186,1.427379228360951,0.9732663879171014,1.2228183863684539,0.9791120953857901,1.3816555082798005,0.7472525445744396,1.2260105818510056,1.0942220810800791,1.5762166528031232,1.1099413389340045,1.5934589028358461,1.434581000171602,1.630307266302407,1.2566113555803897,1.5485478518530726,1.81100521478802,1.3977122131735087,1.7663343697786331,1.162875016592443,1.430508084408939,1.6525550663471222,1.7831156572327018,1.785999870300293,1.6451806398108602,1.4654783833771945,1.2249537132680417,1.3939366523176433,1.0378110969439147,1.3681153113022448,0.9096965923905374,1.0000388860702516,1.2344506287947297,0.8380141992121936,0.8398237504065037,1.0381346080452205,1.4053251046687365,2.167353528179228,1.9828473469242454,1.528205112926662,1.8880876734852792,2.016520819067955,1.5591710267588497,1.5430319080129267,1.4479076098650694,1.236793963611126,1.3141041703522205,1.5537533987313508,1.5692297449335457,1.121472422592342,1.1581335106864572,0.91726502366364,0.9532682375982404,1.2505696462467313,1.1205617936328052,0.6651569686830044,1.0967526400461791,1.140981236472726,0.6660797173157335,0.8100334141403437,1.2571399945765735,1.1800085892900825,0.9563456442207099,1.2749242424964906,1.1340990519151093,0.9276029605418444,0.5047679111361504,0.9407211143523457,0.6309511506929995,1.074688763730228,1.2769522689282895,0.5026340553537011,0.7528656773269176,1.0695354964584112,0.9779949204996229,1.300148520246148,1.1204699540510774,1.5327852481976152,1.4124552657827738,1.2289950992912055,1.0769147641956807,1.4084996918216348,1.2656738823279738,1.0495817463845016,1.1779037110507489,1.095422700792551,1.3921284958720208,0.6597456347197295,1.0181769583374263,0.9364564327523113,1.3884215405210854,0.9055566122755409,1.1202643141150477,1.596531630679965,1.518569000810385,1.201885962113738,1.5418776268139482,1.3921274995431305,0.8939561096951366,1.1407407995313406,1.4240611214190722,1.4680646168068052,1.164111346192658,1.5423512058332562,0.9587342431768777,1.234714117459953,1.4588959828019143,1.6841880340129138,1.6086150078102948,1.4471713358536364,1.2795427026227117,1.444847053103149,1.002940315194428,1.3102364122867585,0.9022639036178589,1.2806951073929669,0.989768156595528,0.9441032137721775,0.9244450947269796,0.9194099459797145,0.9374097302556038,0.7403792202472688,0.8557088484987617,1.1484702169895173,0.7273479776456953,0.8502402873709798,0.6981262018904091,0.8579837949946523,0.7507739521563054,0.7508660649880767,0.7793322708457708,0.9225482998415828,0.9773600473999976,0.9583356784656645,0.9616499662399293,0.9474486894905567,1.0015693807974457,0.9874830450862646,0.9377102639526129,0.9753739142790437,0.9764441156759859,0.973571441695094,0.950638056918979,0.9416109262034296,0.9724794756621122,0.966675561480224,1.2231912130489944,1.235738294199109,0.9654571769759062 NODERIV
  
af_dist2: CONSTANT VALUES=0 NODERIV
Rg: GYRATION TYPE=RADIUS ATOMS=1-414
RMSD1: RMSD REFERENCE=struct1.pdb TYPE=OPTIMAL
RMSD2: RMSD REFERENCE=struct2.pdb TYPE=OPTIMAL
RMSD3: RMSD REFERENCE=struct3.pdb TYPE=OPTIMAL
uwall: UPPER_WALLS ARG=RMSD1,RMSD2,RMSD3 AT=0.02,0.02,0.02 KAPPA=100000,100000,100000 EXP=2,2,2 EPS=1,1,1 OFFSET=0,0,0
PRINT FILE=DISTANCE_MAP_REST ARG=distance_rest_domains.* STRIDE=200
PRINT FILE=COLVAR ARG=__FILL__ STRIDE=200
PBMETAD ...
    LABEL=pb
    ARG=__FILL__
    SIGMA=__FILL__ 
    SIGMA_MIN=__FILL__ 
    SIGMA_MAX=__FILL__ 
    ADAPTIVE=DIFF
    HEIGHT=__FILL__ 
    PACE=200
    BIASFACTOR=__FILL__ 
    GRID_MIN=__FILL__ 
    GRID_MAX=__FILL__ 
    GRID_WSTRIDE=__FILL__ 
    WALKERS_MPI
    TEMP=__FILL__ 
... PBMETAD
METAINFERENCE ...
    ARG=(distance_rest_domains.*),pb.bias REWEIGHT
    PARARG=(af_dist_rest.*)
    SIGMA_MEAN0=1
    NOISETYPE=MGAUSS  OPTSIGMAMEAN=SEM AVERAGING=200
    SIGMA0=10.0 SIGMA_MIN=0.0001 SIGMA_MAX=10.0 DSIGMA=0.1
    MC_STEPS=10
    MC_CHUNKSIZE=23
    WRITE_STRIDE=10000
    TEMP=__FILL__ 
    LABEL=af_mi_rest_domains
... METAINFERENCE
FLUSH STRIDE=200
PRINT FILE=ENERGY ARG=pb.bias STRIDE=200
PRINT ARG=af_mi_rest_domains.*   STRIDE=200 FILE=BAYES_rest_domains
ENDPLUMED
```

Make the PLUMED analysis file ```python plumed_analysis.dat```

```python
#Make the plumed_analysis.dat file
shutil.copy2(dir+"/../scripts_prep/make_plumed_analysis.py", dir)
path_gen_analysis = dir+'/make_plumed_analysis.py sequence.dat'
os.system(f'python {path_gen_analysis}')

```
Similarly the plumed_analysis file is used to calculate the weights using the Torrie valeau weights as done [here](https://link.springer.com/protocol/10.1007/978-1-4939-9608-7_13). For TDP-43 WtoA, the plumed_analysis.dat can be found [here](https://github.com/vendruscolo-lab/AlphaFold-IDP/blob/main/scripts_prep/plumed_analysis_TDP-43.dat)

__Note:__ This tutorial assumes that you know [Parallel Bias Metadynamics](https://www.plumed.org/doc-v2.9/user-doc/html/_p_b_m_e_t_a_d.html) and Metainference [theory](https://link.springer.com/content/pdf/10.1007/978-1-4939-9608-7_13.pdf) and [practice](https://www.plumed.org/doc-v2.9/user-doc/html/_m_e_t_a_i_n_f_e_r_e_n_c_e.html).

### In case you are running the TDP-43 WtoA example:

You can directly use the PBMetaD and Metainference parameters just as copied below by executing the next shell. Otherwise you need to define your system specific parameters.


```python
shutil.copy2(dir+"/../scripts_prep/plumed_TDP-43.dat", dir+'/plumed.dat')
shutil.copy2(dir+"/../scripts_prep/plumed_analysis_TDP-43.dat", dir+'/plumed_analysis.dat')
```

## Short energy minimization using OPENMM and CALVADOS2
```python
shutil.copy2(dir+"/../scripts_prep/simulate_em.py", dir)
simulate_em = dir+'/simulate_em.py '+str(pH)+' '+str(temp)
os.system(f'python {simulate_em}')
```

## Run AF-MI

```python
shutil.copy2(dir+"/../scripts_prep/simulate.py", dir)
!mpirun -np {NR} python simulate.py {pH} {temp}
```
##### [Back to AlphaFold-Metainference home](NAVIGATION.md)
