# CRISPR/Cas AddTag Readme #

Program for identifying exclusive endogenous gRNA sites and creating unique synthetic gRNA sites.

[![Linux](https://img.shields.io/badge/Linux-✓-lightgrey.svg?logo=linux)](#)
[![Windows](https://img.shields.io/badge/Windows-✓-lightgrey.svg?logo=windows)](#)
[![macOS](https://img.shields.io/badge/macOS-✓-lightgrey.svg?logo=apple)](#)
[![Build Status](https://dev.azure.com/tdseher/addtag-project/_apis/build/status/tdseher.addtag-project?branchName=master)](https://dev.azure.com/tdseher/addtag-project/_build/latest?definitionId=1&branchName=master)
[![Code coverage](https://img.shields.io/azure-devops/coverage/tdseher/addtag-project/1)](#)

[![Python](https://img.shields.io/badge/Python-≥3.5.1-1f425f.svg?logo=python)](https://www.python.org/downloads/release/python-360/)
[![downloads](https://img.shields.io/github/downloads/tdseher/addtag-project/total.svg)](https://github.com/tdseher/addtag-project/releases)
[![PRs](https://img.shields.io/badge/PRs-welcome-brightgreen.svg)](http://makeapullrequest.com)
[![](https://img.shields.io/badge/doi-...-blue.svg)](#)

[Features](#-features) • [Requirements](#-requirements) • [Installing](#-installing-addtag) • [Usage](#-program-usage) • [Aligners](#-implemented-sequence-aligners) • [Algorithms](#-implemented-scoring-algorithms) • [Thermodynamics](#-implemented-thermodynamics-calculators) • [Citing](#-citing-addtag) • [Contributing](#-contributing)

## ☑ Features ##
Basic Features:
 * [x] Find ![Target][Target]s with 5'-adjacent PAM (such as Cas12a) or 3'-adjacent PAM (Such as Cas9) sequences
 * [x] Find ![Target][Target]s with arbitrary length ![Spacer][Spacer]s.
 * [x] Find the optimal ![Target][Target] within a ![Feature][Feature] (locus) of interest (multi-allelic, allele-specific, and allele-agnostic).
 * [x] Calculate **on-target** and **off-target** scores (see [Algorithms](#implemented-scoring-algorithms)).
 * [x] Generate unique ![Target][Target]s that don't resemble any genomic DNA (gDNA), thus maximizing **on-target** and **off-target** scores.
 * [x] Find RNA-guided nuclease (![RGN][RGN]) cut sites with arbitrary ![Spacer][Spacer]s using your favorite pairwise alignment program (see [Aligners](#implemented-sequence-aligners)).
 * [x] Perform *in silico* recombination between gDNA and exogenous, donor DNA (![dDNA][dDNA]).
 * [x] Find thermodynamic properties of arbitrary sets of ![Primer][Primer] pairs.
 * [x] See all known ![RGN][RGN] SPACER≷PAM motifs.
 * [x] Analyze gDNA with ambiguous characters or polymorphisms (`RYMKWSBDHVN`).
 * [x] Find ![Spacer][Spacer]s and ![Primer][Primer]s while respecting case-masked gDNA.
 
Advanced Features:
 * [x] Design homology-aware ![dDNA][dDNA]s (multi-allelic, allele-specific, and allele-agnostic).
 * [x] Create ![dDNA][dDNA]s with exceptional ![Target][Target]s while introducing minimal amounts of extrinsic DNA.
 * [x] Engineer a single set of conservative PCR (cPCR) ![Primer][Primer]s that work for all genotypes (wild type, knock-out, and add-back) to validate if a ![Feature][Feature] was engineered correctly.
 
## 📋 Requirements ##

Base operation of AddTag requires the following:

 * Python ≥ 3.5.1 ([source](https://www.python.org/downloads/), [binaries](https://www.python.org/downloads/), [documentation](https://docs.python.org/3/))

 * regex Python module ([source](https://bitbucket.org/mrabarnett/mrab-regex), [whls](https://pypi.org/project/regex/), [documentation](https://pypi.org/project/regex/))
   
     note: The easy way to install this is through `pip`.
     
     ```
     pip3 install regex
     ```

 * BLAST+ ≥ 2.6.0 ([source](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST), [binaries](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST), [documentation](https://www.ncbi.nlm.nih.gov/books/NBK279690/))

 * Bowtie 2 ≥ 2.3.4.1 ([source](https://github.com/BenLangmead/bowtie2), [binaries](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/), [documentation](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml))

Certain optional AddTag functionality depends on the following:

 * Git ≥ 1.7.1 ([source](https://github.com/git/git), [binaries](https://git-scm.com/downloads), [documentation](https://git-scm.com/doc))

For oligo design, AddTag requires one of the following third-party thermodynamics solutions to be installed:

 * UNAFold ≥ 3.8 ([source](http://rnaspace.sourceforge.net/software/unafold-3.8.tar.gz), [documentation](http://unafold.rna.albany.edu/)) with [patch440](http://unafold.rna.albany.edu/?q=node/440)

 * primer3-py Python module ([source](https://github.com/libnano/primer3-py), [whls](https://pypi.org/project/primer3-py/), [documentation](https://libnano.github.io/primer3-py/))
   
     note: The easy way to install this is through `pip`.
     
     ```
     pip3 install primer3-py
     ```

 * ViennaRNA Python module ([source](https://github.com/ViennaRNA/ViennaRNA), [official binaries](https://www.tbi.univie.ac.at/RNA/), [bioconda binaries](https://anaconda.org/bioconda/viennarna/files), [documentation](https://www.tbi.univie.ac.at/RNA/documentation.html))

Certain Target scoring algorithms have additional requirements:

 * Keras Python module ([source](https://github.com/keras-team/keras), [whls](https://pypi.org/project/Keras/), [documentation](https://keras.io/))
     
     note: Keras is only required if you want to calculate CINDEL scores. The easiest way to install is through `pip` also.
     
     ```
     pip3 install Keras
     ```
 
 * Theano Python module ([source](https://github.com/Theano/Theano), [whls](https://pypi.org/project/Theano/), [documentation](http://deeplearning.net/software/theano/))
     
     note: Theano is only required if you want to calculate CINDEL scores. You may install it with `pip` as well.
     
     ```
     pip3 install Theano
     ```
 
 * Azimuth 3 Python module ([source](https://github.com/milescsmith/Azimuth), [documentation](https://www.microsoft.com/en-us/research/project/crispr/))
   
     note: Either Azimuth2 or Azimuth 3 is only required if you want to calculate Azimuth scores.
     
     note: requires specific versions of numpy, scikit-learn, and pandas.
     Other dependencies include click, biopython, scipy, GPy, hyperopt, paramz, theanets, glmnet_py, dill, matplotlib, pytz, python-dateutil, six, tqdm, future, networkx, pymongo, decorator, downhill, theano, nose-parameterized, joblib, kiwisolver, cycler, pyparsing, setuptools, glmnet-py.
     
     You can install it using `pip`, assuming `git` is available in the `PATH` environmental variable.
     
     ```
     pip3 install git+https://github.com/milescsmith/Azimuth.git
     ```

 * 3.0.0 > Python ≥ 2.7.10 ([source](https://www.python.org/downloads/), [binaries](https://www.python.org/downloads/), [documentation](https://docs.python.org/2/))
   with Azimuth 2 Python module ([source](https://github.com/MicrosoftResearch/Azimuth), [documentation](https://www.microsoft.com/en-us/research/project/crispr/))
      
     note: Either Azimuth2 or Azimuth 3 is only required if you want to calculate Azimuth scores.
     
     note: requires python-tk to be installed. Also requires specific versions of scipy, numpy, matplotlib, nose, scikit-learn, pandas, biopython, pyparsing, cycler, six, pytz, python-dateutil, functools32, subprocess32.
     
     The easiest way to install it and take care of all dependencies is to use `pip`, assuming `git` is available in the `PATH` environmental variable.
     
     ```
     pip2.7 install git+https://github.com/MicrosoftResearch/Azimuth.git
     ```

<!--
##### The following third-party integrations are currently incomplete #####

For speed, we recommend at least one third-party pairwise nucleotide sequence alignment program:

 * BLAT ([source](https://genome.ucsc.edu/goldenPath/help/blatSpec.html), [binaries](http://hgdownload.cse.ucsc.edu/admin/exe/), [documentation](https://genome.ucsc.edu/goldenPath/help/blatSpec.html))

 * Bowtie ([source](https://github.com/BenLangmead/bowtie), [binaries](https://sourceforge.net/projects/bowtie-bio/files/bowtie/), [documentation](http://bowtie-bio.sourceforge.net/manual.shtml))

 * BWA ([source](https://github.com/lh3/bwa), [documentation](http://bio-bwa.sourceforge.net/bwa.shtml))

 * Cas-OFFinder ([source](https://github.com/snugel/cas-offinder), [binaries](https://sourceforge.net/projects/cas-offinder/files/Binaries/), [documentation](http://www.rgenome.net/cas-offinder/portable))
-->

## ⤵ Installing AddTag ##
You can download the latest version of AddTag over HTTPS using `git` with the following command.
```sh
git clone https://github.com/tdseher/addtag-project.git
```

This will download AddTag into a folder called `addtag-project/` in your current working directory. Go ahead and change the working directory into the AddTag folder.
```sh
cd addtag-project/
```

`git` should automatically make the `addtag` program executable. If it does not, you can use the following command to do it. 
```sh
chmod +x addtag
```

## 🔁 Updating AddTag ##
The commands in this section assume the working directory is the AddTag folder.
```sh
cd addtag-project/
```

If you would like to update your local copy to the newest version available, use the following command from within the `addtag-project/` directory.
```sh
./addtag update
```

If you want the newest version, but you made changes to the source code, then you can first discard your changes, and then update. Use the following command from inside the `addtag-project/` folder.
```sh
./addtag update --discard_local_changes
```

Alternatively, if you want to keep the local modifications, you can use the `--keep_local_changes` option to stash, pull, then reapply them afterwards.
```sh
./addtag update --keep_local_changes
```

Each one of these methods uses `git`, and may require your Atlassian login credentials.

## 💻 Program usage ##
Because AddTag is being updated regularly, the most current feature set and usage can be viewed by running AddTag with the `--help` command line option.

The following commands assume the current working directory is the AddTag folder `addtag-project/`. This will print out command line parameter descriptions and examples.
```sh
./addtag --help
```

Additionally, you may view the included man page, which is probably not up-to-date.
```sh
man ./addtag.1
```

## 📈 Implemented scoring Algorithms ##

The following scoring algorithms are subclasses of `SingleSequenceAlgorithm`.
     
 * Azimuth ([Doench, Fusi, et al (2016)](http://dx.doi.org/10.1038/nbt.3437))
 * Doench-2014 ([Doench, et al (2014)](http://dx.doi.org/10.1038/nbt.3026))
 * GC
 * PolyT
 * Housden ([Housden, et al (2015)](http://dx.doi.org/10.1126/scisignal.aab3729))
 * Moreno-Mateos ([Moreno-Mateos, et al (2015)](http://dx.doi.org/10.1038/nmeth.3543))
 * DeepCpf1/CINDEL ([Kim, Song, et al (2016)](http://dx.doi.org/10.1038/nmeth.4104))
 * PAM Identity

The following scoring algorithms are subclasses of `PairedSequenceAlgorithm`.

 * CFD ([Doench, Fusi, et al (2016)](http://dx.doi.org/10.1038/nbt.3437))
 * Substitutions, Insertions, Deletions, Errors ([Needleman, Wunsch (1970)](https://dx.doi.org/10.1016/0022-2836%2870%2990057-4))
 * Hsu-Zhang ([Hsu, et al (2013)](http://dx.doi.org/10.1038/nbt.2647))
 * Linear
 * CRISPRater ([Labuhn, et al. (2018)](http://dx.doi.org/10.1093/nar/gkx1268))

## 📐 Implemented sequence Aligners ##

 * Bowtie2
 * BLAST+

## 🌡 Implemented thermodynamics calculators ##

 * Primer3
 * UNAFold
 * ViennaRNA

## 📝 Citing AddTag ##
If you use AddTag for your research, please cite us.

## ✍ Authors ##
Who do I talk to?
 * Aaron Hernday (🔬 PI leading the project)
 * Thaddeus D. Seher (💻 programmer) (💬[@tdseher][tdseher])

See also the list of [contributors](https://github.com/tdseher/addtag-project/graphs/contributors) who participated in this project.

## 👥 Contributing ##
<details>
  <summary>🐞 How do I submit a bug report?</summary>

  First, check to see if the problem you are having has already been added to the [issue tracker](https://github.com/tdseher/addtag-project/issues).
  If not, then please submit a new issue.
</details>

<details>
  <summary>⚠ How do I make a feature request?</summary>

  Send a message to [@tdseher][tdseher].
</details>

<details>
  <summary>⤴ How do I add my code to the AddTag software?</summary>

  Please submit a [pull request](https://github.com/tdseher/addtag-project/pulls).
</details>

### Adding scoring Algorithms ###
Scoring Algorithms have been broken down into two general types.

 * `SingleSequenceAlgorithm` objects calculate scores by comparing a potential ![Spacer][Spacer] to a model trained on empirical data.
 * `PairedSequenceAlgorithm` instances generate scores that compare a potential ![Spacer][Spacer] to a target using a model.
 
To add a new scoring algorithm, you must subclass one of the the above types, and add it to a `*.py` file in the `source/algorithms/` subdirectory. AddTag will automatically calculate the score on every generated ![Spacer][Spacer].

We welcome any `git pull` requests to widen the repertoire of scoring algorithms available to AddTag. The easiest way to get started is to copy and modify one of the provided subclasses.

### Adding sequence Aligners ###
AddTag comes with wrappers for several alignment programs. Depending on your experimental design and computing system, you may decide to use an aligner with no included wrapper. To implement your own, create a subclass of `Aligner`, and put it in a `*.py` file in the `source/aligners/` subdirectory. AddTag will automatically make that aligner available for you.

Share your code with us so we can make it available to all AddTag users.

### Adding Thermodynamics calculators ###
Several wrappers to popular oligonucleotide conformation, free energy, and melting temperature calculation programs are included. You can add your own by subclassing the `Oligo` class, and then adding its `*.py` file to the `source/thermodynamics/` subdirectory.

If you create your own wrapper, please submit a `git pull` request so we can add it to the next version of the software.

## 📖 License ##
Please see the [LICENSE.md](LICENSE.md) file.

## Notes ##
Here is some miscellaneous information:
 * The ![RGN][RGN] protein you use should be engineered specifically for your organism. It should be codon-optomized, and if using eukarya, contain an appropriate nuclear localization sequence.
 * Sequences in FASTA files should have unique names. In other words, the primary sequence header--everything following the '`>`' character and preceding the first whitespace/tab '` `' character--should exist only once across all input `*.fasta` files.
 * By default, AddTag will avoid designing homology regions and Targets against polymorphisms whenever possible.

[tdseher]:https://twitter.com/tdseher
[Spacer]:docs/spacer.svg
[Target]:docs/target.svg
[dDNA]:docs/ddna.svg
[RGN]:docs/rgn.svg
[Feature]:docs/feature.svg
[Primer]:docs/primer.svg