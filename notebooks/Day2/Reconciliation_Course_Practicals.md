# Gene Tree - Species Tree Reconciliation
Created by: Paul Simion, Rimjhim Choudhury
Date: 2023-07-15

## Table of Contents
- [Input Dataset](#Input-Dataset)
- [Practical 1: building gene trees](#Practical-1:-building-gene-trees)
- [Practical 2: Comparing trees](#Practical-2:-Comparing-trees)
- [Practical 3: Reconciling trees](#Practical-3:-Reconciling-trees)
- [Useful Information](#Useful-Information)


################
## Input Dataset

**Summary of the input Dataset production protocol:**

A set of 28 species have been selected and their proteomes have been downloaded from then *ENSEMBL plants* ftp repository. Then, *Orthofinder2* was used in order to build sets of homologous sequences (i.e. both orthologous and paralogous). These sets of homologous sequences have then been filtered based on number of species and number of sequences thresholds (in order to facilitate the subsequent analyses, gaining time and using less electricity). For example, if a folder is named `20-28seq_20sp`, this means that the clusters of homologous sequences contain between 20 to 28 sequences and at least 20 different species.

input files location: `/data/courses/pploidy/reconciliation/`
REC_env environment location: `/home/psimion/.conda/envs/REC_env`


###################################
## Practical 1: building gene trees

**Software tools needed:**
- *mafft* (https://mafft.cbrc.jp/alignment/software/)
- *iqtree* (http://www.iqtree.org/doc/)

**Web-based tree visualization tool:**
- *iTOL* (https://itol.embl.de/upload.cgi)

**Protocol**
Using input fasta files of unaligned sequences, we will use *mafft* to align the sequences and then use *iqtree* to infer the gene tree based on the aligned sequences. input fasta files are located here: `/data/courses/pploidy/reconciliation/Practical1/orthofinder_output_seqs_20-28seq_20sp`. 

1. Randomly select one of these fasta files and proceed with the preparation of your script that you will use to launch jobs on the cluster.
```
#!/bin/bash
#SBATCH --job-name=single_gene_inference
#SBATCH --output=single_gene_inference_%j.out
#SBATCH --error=single_gene_inference_%j.err
#SBATCH --time=00:30:00			# max runtime allowed
#SBATCH --mem=1G			     	# max allocated memory 
#SBATCH --cpus-per-task=4 		# number of threads
#SBATCH --partition=pploidy 	# which cluster partition to use

module load Conda/miniconda/latest
conda activate REC_env

# Multiple Sequence Alignment (with mafft)
mafft --thread 4 --auto --anysymbol [InputFile] > [AlignedFile]
# Phylogenetic Inference (with iqtree)
iqtree -T 4 -s [AlignedFile] -b 10 -m [MODEL] -pre [OUPUT]

conda deactivate
   ```

2. Then, submit the job using the following command: `sbatch [NameOfYourScript]`. Among the output files, the finalized gene tree infered has a `.treefile` extension. You can print its content using `cat [OUTPUT.treefile]`, and then copy-paste it on the web-server tool *iTOL* in order to get a nice visualization of the gene tree: https://itol.embl.de/upload.cgi. 

3. Now, infer and visualize the gene tree for one randomly chosen gene from the `/data/courses/pploidy/reconciliation/Practical1/orthofinder_output_seqs_20-80seq_20sp` folder. For this, you will have to re-do steps 1 and 2 above and changing the `InputFile`.

4. Bonus: Many different sequence evolution models exists from simplistic to quite complex ones (http://www.iqtree.org/doc/Substitution-Models#protein-models, http://www.iqtree.org/doc/Complex-Models). Feel free to test different models to see if (and how) they impact the resulting gene trees. The simplest and least realistic model is the `JC69` model. The evolution model `GTR20+C60+R8+F*H4` is an example of a very complex one. Frequently used models include `LG+G4+F` and `GTR20+R4+F`.

**Discussion**
Take some time to look at the gene(s) tree(s) you obtained. Are they identical with the known species tree of plants ? See the reference species tree here: http://www.mobot.org/MOBOT/research/APweb/. How come we only have 50 genes (among the entire proteome of plants) in the `orthofinder_output_seqs_20-28seq_20sp` folder ? Interestingly, orthofinder outputs zero single-copy orthologous sequences, why ?



###############################
## Practical 2: Comparing trees

**Software tools needed:**
- *tqdist* (https://birc.au.dk/~cstorm/software/tqdist/)

**Protocol**

1. First, infer the phylogenetic tree of both alignments `gene_A.fas` and `gene_B.fas`, as in Practical 1 (except that the sequences have already been aligned). They are localized in the folder `/data/courses/pploidy/reconciliation/Practical2/`. After observing their gene trees using *iTOL*, find the incongruences between them, and discuss what could have caused these incongruences.

2. Comparing the topological distance between trees: we will compute the normalized quartet distance between the gene tree `gene_C.tre` (already infered from OG0010236) and the two species tree given `species_tree_1.tre` and `species_tree_2.tre` . For this, we will use *tqdist* and its `quartet_dist` function. Don't forget, as in Practical 1, to load the conda module and to activate the `REC_env` conda environment before running tdqist.
```
quartet_dist -v [treefile1] [treefile2]
   ```
Which species tree is better *supported* by the gene tree ?

3. Comparing the likelihood of the two trees given the gene sequence data (i.e. the gene alignment). We will use *iqtree* to infer two *constrained* trees. Lastly, we will also need to put the two trees under investigation in the same file. Again, don't forget to load the conda module and to activate the `REC_env` conda environment :)
```
cat [tree1] [tree2] > [trees_file] 
iqtree -T 4 -s [AlignedFile] -m [MODEL] -g [tree1] -pre [OUPUT1]
iqtree -T 4 -s [AlignedFile] -m [MODEL] -g [tree2] -pre [OUPUT2]
iqtree -T 4 -s [AlignedFile] -z [trees_file] -au
   ```
Again, which species tree is better *supported* by the gene tree ?

4. We now will check whether the previous test was sensible (or not!). For topological tests to make sense, constrained tree must not be rejected significantly when compared to the *unconstrained* tree. Luckily, you already computed the unconstrained tree in step 1 above, so you simply have to append it at the end of your `[trees_file]`
```
iqtree -T 4 -s [AlignedFile] -z [trees_file_new] -au
   ```
Was the comparison between the 2 constrained topology OK/sensible/legitimate ? 




#################################
## Practical 3: Reconciling trees

**Software tools needed:**
- *GeneRax* (https://github.com/BenoitMorel/GeneRax/wiki/GeneRax)
- *Treerecs* (https://project.inria.fr/treerecs/treerecs-options/)

**Protocol**

1. Inference of a species tree based using a set of gene trees. We pre-computed phylogenetic trees for some of the (many) sets of sequences outputed by *orthofinder2*. These 2,034 gene trees are located in the `Orthogroups_trees_60-80seq_25sp` folder. Note that you are now able to infer all these trees by yourself :)
In order to us *GeneRax*, and particularly its *SpeciesRax* function, we need to build the necessary input files that the software needs. First, we will produce the `family` file (here, using only input gene trees). For this, you will have to download the 'build_family_file.py' from *GeneRax* github page (https://github.com/BenoitMorel/GeneRax).

```
tree_folder=Orthogroups_trees_60-80seq_25sp
output_file=orthogroups_family_2034.file
python /path/to/GeneRax/scripts/build_family_file.py \
	NONE \
	$tree_folder \
	NONE \
	NONE \
	$output_file
   ```
Then, we can run the *SpeciesRax* function with the following command:
```
family_file=orthogroups_family_2034.file
output=orthogroups_family_2034_speciesrax
/path/to/GeneRax/build/bin/generax --families $family_file \
	--rec-model UndatedDTL --prefix $output \
	--strategy SKIP --si-strategy HYBRID --species-tree MiniNJ \
	--si-estimate-bl --si-quartet-support --prune-species-tree --per-family-rates
   ```
Here, the `--strategy SKIP --si-strategy HYBRID` parameterization indicate to GeneRax that the species tree is not given as an input and that reconciliation modeling will be used to infer the species tree. A first approximation of the species tree will be produce (before being refined during the run) by *MiniNJ*.

The output species tree is then localised here: `orthogroups_family_2034_speciesrax/species_trees/species_tree_quartet_support.newick`.


2. Reconciling gene trees with the newly obtained species tree using *Treerecs* (which allow for an easy visualization of the results). This software takes as input a species tree (`[InputSpeciesTree]`) and a file containing all gene trees to be reconciled (`[InputGeneTrees]`). Its use is as follows:
```
# preparing input gene trees (and adding a newline after every tree)
cat orthofinder_output_trees_20-28seq_20sp/* | sed 's/;/;\n/g' > [InputGeneTrees]

# using treerecs
treerecs -g [InputGeneTrees] -s [InputSpeciesTree] -o [OutputName] -r -f -C --fevent -O recphyloxml:svg:nhx:newick
   ```
Treerecs ouputs directly reconciliated trees as vectorial graphs (.svg extension). the ouput file ending with the extension `.tre_recs.nhx` can be opened using *iTOL* and allow to display branch metadata to visualise the location of duplication in the gene tree (data source=D).

Finally, *seaview* allows to use *Treerecs* directly on a given gene and also outputs a visualization of the results. But this can only be done one gene at a time. 



#####################
## Useful Information

# Data access and Resource allocation
Please note the following information regarding data access and resource allocation on the cluster:

1. Data Access:
   - The data for the Summer school practical sessions are located in the directory: /data/courses/pploidy/
   - Students are in the Unix group `pploidystudents` have read-only permissions for these data.

2. User Quotas:
   - Each student and user has a 4TB quota to write into the directory /data/users/loginname.
   - Each student and user has a 20GB quota to write into the directory /home/loginname.
   - Use the `lsquota` command to view the individual quota usage.

3. Cluster Resources:
   - The Slurm partition for the cluster is named `pploidy`.
   - The `pploidy` partition has 128 CPUs reserved.
   - The `pploidy` partition has 512 GB of RAM.

# Submitting and Monitoring jobs on the IBU cluster
You can create a job script to submit jobs on the cluster. The job script contains the commands that you would like to run on the cluster. You can submit the job script using the `sbatch` command. The job will be added to the queue and will be executed when the resources are available. You can check the status of the job using the `squeue` or `sacct` commands. You can cancel the job using the `scancel` command.

1. Example on how to create a job submission script (e.g. `job_to_submit.sh`), which is a text file with content such as follows:
   ```
   #!/bin/bash

   # set several important sbatch variable
   #SBATCH --job-name=polyploid_genotyping
   #SBATCH --output=polyploid_genotyping_%j.out
   #SBATCH --error=polyploid_genotyping_%j.err
   #SBATCH --time=00:30:00
   #SBATCH --mem=10G
   #SBATCH --cpus-per-task=1
   #SBATCH --partition=pploidy

   # Load conda module
   module load Conda/miniconda/latest

   # Activate a Conda environment
   conda activate [EnvironmentName]

   # Run the command (archetypical example)
   ToolName -i [InputFile] -t [NumberOfThreads] -o [Outputfile]
   ```
2. Submit the job using the `sbatch` command:
   ```
   sbatch submit.sh
   ```
3. Check the status of the job using the `squeue` command:
   ```
   squeue -u <your-username>
   ```
4. Check the status of the job using the `sacct` command:
   ```
   sacct -j <job-id>
   ```
5. Cancel the job using the `scancel` command:
   ```
   scancel <job-id>
   ```
