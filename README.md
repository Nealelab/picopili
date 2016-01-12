# picopili

### Purpose

* Expand [ricopili](https://github.com/Nealelab/ricopili) to support GWAS of family data, from sib-pairs to extended pedigrees

### New Primary Functions

* `qc_rel.py`: runs QC for family GWAS data. Designed to be consistent with ricopili trio QC, with minor adjustments to mendelian error handling and streamlined output.
* `pca_rel.py`: main driver script to run PCA on family GWAS data. Includes strict QC (with LD pruning, etc), extraction of a set of unrelated individuals to compute PCs, and projection of PCA to full dataset. Jobs are submitted via UGER, and ricopili-like success/failure emails are sent on completion. See tasks `strict_qc.py` and `imus_pca.py`.
* `admix_rel.py`: estimates relatedness for admixed samples. Relatedness is estimated by starting with a subset of individuals who are unrelated, running an Admixture analysis, selecting "exemplar" individuals for each ancestry component, using those individuals to anchor a supervised Admixture analysis fo the full data, and using those Admixture results as the basis for ancestry-adjusted relatedness estimation using REAP.
* `filter_ped.py`: uses genetic relatedness information to flag (a) cryptic relatedness across FIDs, (b) unrelated individuals within FIDs, (c) possible parent/offspring pair not identified in .fam file, and (d) parent/offspring pairs indicated by pedigree that aren't supported by genetic relatedness. Provides a suggested sample exclusion list based on weighted preferences for phenotype, pedigree, and genotyping rate.
* `gwas_rel.py`: Genome-wide association analysis (in parallelized chunks) with either `plink --dfam` or a GEE model.
* `impute_rel.py` (Coming soon): Imputation pipeline for related samples using SHAPEIT's `--duohmm` and IMPUTE2. Includes build check (with liftOver if needed), alignement to reference, phasing, imputation, and best-guess genotype calls with MAF/info score/mendelian error filtering.


### Additional Changes

* Ricopili bugfix: avoid (rarely) losing SNPs when merging chunk results (see [dameta_cat](https://github.com/Nealelab/picopili/blob/dev/rp_bin/dameta_cat))
* New task `strict_qc.py`: runs the stricter QC used for PCA/relatedness. Compared to QC included in `pcaer_20`, includes additional options for how long LD regions, indels, and strand ambiguous SNPs are handled.
* New task `imus_pca.py`: extracts a set of unrelated individuals (without relying on pedigrees and partially controlling for ancestry-based confounding), runs PCA, and projects the computed PCs back to the remaining (related) individuals. Assumes appropriate QC has already been run.
* New task `ped_confirm.py`: confirms that reported pedigrees in .fam file are consistent with genetic relatedness. Works as a wrapper to `find_expected_pedigree.pl` from PRIMUS.
* Additional component tasks for imputation/association/etc.


### Install alongside Ricopili

##### 1. The following assumes:

* You're working on the Broad UGER cluster (though can probably be adapted)

* You've already installed [ricopili](https://github.com/Nealelab/ricopili).

* You have git configured with a github account.

Note additional dependencies may exist for some modules (e.g. Admixture, REAP, PRIMUS).

##### 2. Download picopili from github

```
mkdir ~/github
cd ~/github
git clone https://github.com/Nealelab/picopili.git
git checkout dev
```

Note: the 'dev' branch is currently intended to be the stable picopili branch, with the 'master' branch reserved for pulling updates from new releases of ricopili.

##### 3. Create dotkits for ricopili and picopili

This will allow switching between ricopili and picopili as desired.

First, setup a folder for personal dotkits.

```
mkdir ~/.kits
cd ~/.kits
```

Next, create a dotkit for picopili.

```
echo "
#c personal
#d fork of ricopili for family data, loaded from github folder

dk_alter PATH $HOME/github/picopili/rp_bin
dk_alter PATH $HOME/github/picopili/rp_bin/pdfjam
" > picopili.dk
```

Also create a dotkit for ricopili (change the install location as needed)

```
echo "
#c personal
#d gwas pipeline

dk_alter PATH $HOME/rp_bin
dk_alter PATH $HOME/rp_bin/pdfjam
" > ricopili.dk
```

You should now have entries for ricopili and picopili in you list of available dotkits (`use -la`)

##### 4. Switch how ricopili loads

In `~/.my.bashrc` (or equivalent for your shell) remove the lines

```
PATH=/home/unix/<username>/rp_bin:$PATH
PATH=/home/unix/<username>/rp_bin/pdfjam:$PATH
```

(or similar, depending on where you've installed ricopili), and replace them with


and replace them with

```
use ricopili
```

##### 5. Logout and log back in.

Applies the changes in `~/.my.bashrc`.

### Using ricopili/picopili

Once the above install is done, analyses will still be done with ricopili by default. 

*  **To use picopili:** change `use ricopili` to `use picopili` in your `~/.my.bashrc` (or equivalent) and run

```
unuse ricopili
use picopili
```

* **To use ricopili:** change `use picopili` to `use ricopili` in your `~/.my.bashrc` (or equivalent) and run

```
unuse picopili
use ricopili
```

