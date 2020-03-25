# CRISPR/Cas AddTag Readme #

Program for identifying exclusive endogenous gRNA sites and creating unique synthetic gRNA sites.

[![Linux](https://img.shields.io/badge/Linux-✓-darkgreen.svg?logo=linux)](#)
[![Windows](https://img.shields.io/badge/Windows-✓-darkgreen.svg?logo=windows)](#)
[![macOS](https://img.shields.io/badge/macOS-✓-darkgreen.svg?logo=apple)](#)
[![Build Status](https://dev.azure.com/tdseher/addtag-project/_apis/build/status/tdseher.addtag-project?branchName=master)](https://dev.azure.com/tdseher/addtag-project/_build/latest?definitionId=1&branchName=master)
[![Code coverage](https://img.shields.io/azure-devops/coverage/tdseher/addtag-project/1)](#)

[![Python](https://img.shields.io/badge/Python-≥3.5.1-1f425f.svg?logo=python)](https://www.python.org/downloads/release/python-360/)
[![downloads](https://img.shields.io/github/downloads/tdseher/addtag-project/total.svg)](https://github.com/tdseher/addtag-project/releases)
[![PRs](https://img.shields.io/badge/PRs-welcome-brightgreen.svg)](http://makeapullrequest.com)
[![](https://img.shields.io/badge/doi-...-blue.svg)](#)

[Features](#-features) • [Requirements](#-requirements) • [Installing](#-installing-addtag) • [Usage](#-program-usage) • [Aligners](#-supported-sequence-aligners) • [Thermodynamics](#-supported-thermodynamics-calculators) • [Algorithms](#-supported-scoring-algorithms) • [Citing](#-citing-addtag) • [Contributing](#-contributing)

## ☑ Features ##
Basic Features:
 * [x] Analyzes any arbitrary genomic DNA (gDNA).
   * [x] Fully supports ambiguous characters or polymorphisms (`RYMKWSBDHVN`).
   * [x] Respects case-masked gDNA for ![Target][Target] and ![Primer][Primer] identification.
 * [x] Locates RNA-guided nuclease (![RGN][RGN]) cut sites (![Target][Target]s) within a ![Feature][Feature] (locus of interest) for optimal gRNA ![Spacer][Spacer]s.
   * [x] Fully supports ambiguous bases (`RYMKWSBDHVN`) in ![Spacer][Spacer] or ![PAM][PAM].
   * [x] Accepts 3'-adjacent ![PAM][PAM] sequences, such as Cas9 (`>NGG`).
   * [x] Accepts 5'-adjacent ![PAM][PAM] sequences, such as Cas12a (`TTTN<`).
   * [x] Supports arbitrary ![Spacer][Spacer] length and composition constraints (`G{,2}N{19,20}`).
   * [x] Supports arbitrary ![PAM][PAM] sequences (MAD7: `YTTN<`, Cas12d: `TA<`, BlCas9: `>NGGNCNDD`, etc).
   * [x] Uses stranded forward (`/`), reverse (`\`) and unstranded (`|`) cut sites.
   * [x] Calculates any number of **on-target** and **off-target** scores (see [Algorithms](#supported-scoring-algorithms)).
   * [x] Finds homology-aware ![Target][Target]s (**multi-allelic**, **allele-specific**, and **allele-agnostic**).
   * [x] Searches for ![Target][Target]s using selectable pairwise alignment program (see [Aligners](#supported-sequence-aligners)).
 * [x] Generates exogenous, donor DNAs (![dDNA][dDNA]s) to modify the same locus successively.
   * [x] Assembles unique ![Target][Target] sites (on ![dDNA][dDNA]s), thus maximizing **on-target** and **off-target** scores (because they don't resemble any input gDNA).
   * [x] Adds unique ![Target][Target]s to ![dDNA][dDNA]s without inserting sequences (or while introducing minimal amounts of extrinsic DNA).
   * [x] Engineers a single set of conservative PCR (cPCR) ![Primer][Primer]s that work for all genotypes (wild type, knock-out, and add-back) to validate if a ![Feature][Feature] was engineered correctly.
   * [x] Produces homology-aware ![dDNA][dDNA]s (**multi-allelic**, **allele-specific**, and **allele-agnostic**).
 * [x] Performs *in silico* recombination between gDNA and ![dDNA][dDNA]s.
 * [x] Determines thermodynamic properties of sets of ![Primer][Primer] pairs (Tm, minimum ΔG, amplicon size, etc).
 * [x] Displays all known ![RGN][RGN] ![Spacer][Spacer] and ![PAM][PAM] combinations.

## 📋 Requirements ##

Below are lists AddTag requirements. Each entry is marked with a ☑ or ☐, indicating whether or not an additional download/setup is required:

 * [x] All requirements included in AddTag
 * [ ] Additional download/setup required

### Basic requirements ###

Base operation of AddTag requires the following:

 * [ ] Python ≥ 3.5.1 ([source](https://www.python.org/downloads/), [binaries](https://www.python.org/downloads/), [documentation](https://docs.python.org/3/))

 * [ ] regex Python module ([source](https://bitbucket.org/mrabarnett/mrab-regex), [whls](https://pypi.org/project/regex/), [documentation](https://pypi.org/project/regex/))

Certain optional AddTag functionality (version information, and software updates) depends on the following:

 * [ ] Git ≥ 1.7.1 ([source](https://github.com/git/git), [binaries](https://git-scm.com/downloads), [documentation](https://git-scm.com/doc))

### 📐 Supported sequence Aligners ###

One pairwise sequence aligner is required:

 * [ ] BLAST+ ≥ 2.6.0 ([source](https://bit.ly/2Ouoqkx), [binaries](https://bit.ly/2Ouoqkx), [documentation](https://www.ncbi.nlm.nih.gov/books/NBK279690/))

 * [ ] Bowtie 2 ≥ 2.3.4.1 ([source](https://github.com/BenLangmead/bowtie2), [binaries](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/), [documentation](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml))
 
 * [ ] BWA ≥ 0.7.12 ([source](https://github.com/lh3/bwa), [ugene binaries](http://ugene.net/download-all.html#en_data_analysis_tools), [bioconda binaries](https://anaconda.org/bioconda/bwa/files), [documentation](http://bio-bwa.sourceforge.net/bwa.shtml))
 
 * [ ] Cas-OFFinder ≥ 2.4 ([source](https://github.com/snugel/cas-offinder), [binaries](https://sourceforge.net/projects/cas-offinder/files/Binaries/), [documentation](http://www.rgenome.net/cas-offinder/portable))

<!--
##### The following third-party integrations are currently incomplete #####

For speed, we recommend at least one third-party pairwise nucleotide sequence alignment program:

 * [ ] BLAT ([source](https://genome.ucsc.edu/goldenPath/help/blatSpec.html), [binaries](http://hgdownload.cse.ucsc.edu/admin/exe/), [documentation](https://genome.ucsc.edu/goldenPath/help/blatSpec.html))

 * [ ] Bowtie ([source](https://github.com/BenLangmead/bowtie), [binaries](https://sourceforge.net/projects/bowtie-bio/files/bowtie/), [documentation](http://bowtie-bio.sourceforge.net/manual.shtml))

 * [ ] Usearch
-->

For polymorphism-aware ![Feature][Feature] expansion, one multiple sequence aligner is required:

 * [ ] MAFFT ([source](https://mafft.cbrc.jp/alignment/software/source.html), [binaries](https://mafft.cbrc.jp/alignment/software/), [documentation](https://mafft.cbrc.jp/alignment/software/manual/manual.html))

### 🌡 Supported thermodynamics calculators ###

For oligo design, AddTag requires one of the following third-party thermodynamics solutions to be installed:

 * [ ] UNAFold ≥ 3.8 ([source](http://rnaspace.sourceforge.net/software/unafold-3.8.tar.gz), [documentation](http://unafold.rna.albany.edu/)) with [patch440](http://unafold.rna.albany.edu/?q=node/440)

 * [ ] primer3-py Python module ([source](https://github.com/libnano/primer3-py), [whls](https://pypi.org/project/primer3-py/), [documentation](https://libnano.github.io/primer3-py/))

 * [ ] ViennaRNA Python module ([source](https://github.com/ViennaRNA/ViennaRNA), [official binaries](https://www.tbi.univie.ac.at/RNA/index.html#download), [bioconda binaries](https://anaconda.org/bioconda/viennarna/files), [documentation](https://www.tbi.univie.ac.at/RNA/documentation.html))

### 📈 Supported scoring Algorithms ###

The following scoring algorithms are subclasses of `SingleSequenceAlgorithm`.

 * [ ] Azimuth ([Doench, Fusi, et al (2016)](http://dx.doi.org/10.1038/nbt.3437))

   note: Either Azimuth 2 or Azimuth 3 can be used to calculate Azimuth scores. There is no need to have both installed.
   
    * Azimuth 3 Python module ([source](https://github.com/milescsmith/Azimuth), [documentation](https://www.microsoft.com/en-us/research/project/crispr/))
      
      note: requires specific versions of numpy, scikit-learn, and pandas.
      Other dependencies include click, biopython, scipy, GPy, hyperopt, paramz, theanets, glmnet_py, dill, matplotlib, pytz, python-dateutil, six, tqdm, future, networkx, pymongo, decorator, downhill, theano, nose-parameterized, joblib, kiwisolver, cycler, pyparsing, setuptools, glmnet-py.

    * Azimuth 2 Python module ([source](https://github.com/MicrosoftResearch/Azimuth), [documentation](https://www.microsoft.com/en-us/research/project/crispr/))
      on 2.7.10 ≤ Python < 3.0.0 ([source](https://www.python.org/downloads/), [binaries](https://www.python.org/downloads/), [documentation](https://docs.python.org/2/)) 
      
      note: requires python-tk to be installed. Also requires specific versions of scipy, numpy, matplotlib, nose, scikit-learn, pandas, biopython, pyparsing, cycler, six, pytz, python-dateutil, functools32, subprocess32.

 * [ ] DeepCpf1/CINDEL ([Kim, Song, et al (2016)](http://dx.doi.org/10.1038/nmeth.4104))
       
   note: Requires both Keras and Theano Python modules.
   
    * Keras Python module ([source](https://github.com/keras-team/keras), [whls](https://pypi.org/project/Keras/), [documentation](https://keras.io/))
    * Theano Python module ([source](https://github.com/Theano/Theano), [whls](https://pypi.org/project/Theano/), [documentation](http://deeplearning.net/software/theano/))
 
 * [x] Doench-2014 ([Doench, et al (2014)](http://dx.doi.org/10.1038/nbt.3026))
 
 * [x] Housden ([Housden, et al (2015)](http://dx.doi.org/10.1126/scisignal.aab3729))
 
 * [x] Moreno-Mateos ([Moreno-Mateos, et al (2015)](http://dx.doi.org/10.1038/nmeth.3543))
 
 * [x] GC
 
 * [x] PolyT
  
 * [x] PAM Identity

The following scoring algorithms are subclasses of `PairedSequenceAlgorithm`.

 * [x] CFD ([Doench, Fusi, et al (2016)](http://dx.doi.org/10.1038/nbt.3437))
 
 * [x] Substitutions, Insertions, Deletions, Errors ([Needleman, Wunsch (1970)](https://dx.doi.org/10.1016/0022-2836%2870%2990057-4))
 
 * [x] Hsu-Zhang ([Hsu, et al (2013)](http://dx.doi.org/10.1038/nbt.2647))
 
 * [x] CRISPRater ([Labuhn, et al. (2018)](http://dx.doi.org/10.1093/nar/gkx1268))
 
 * [x] Linear

### Python package setup ###
There are several standard ways to make modules available to your Python installation. The easy way to install a package this is through `pip`.

For example, the following code will download and setup the `regex` package from [PYPI](https://pypi.org/) into your default Python installation.
```sh
pip install regex
```

If you want to make the module available to a specific Python installation, use a command like this:
```sh
/path/to/python -m pip install regex
```

Often, the package is not available on PYPI, or you need a development version. In these cases, you can direct `pip` to download and setup a package from a code repository. The easiest way to install it and take care of all dependencies is to use `pip`, assuming `git` is available in the `PATH` environmental variable. Here is how to install the `Azimuth` package from [GitHub](https://github.com/).
```sh
pip2.7 install git+https://github.com/MicrosoftResearch/Azimuth.git
```

Some Python packages are available through [bioconda](https://anaconda.org/bioconda). To install `viennarna` using `conda`, use this command:
```sh
conda install -c bioconda viennarna
```

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

Each one of these commands assumes `git` is available on the `PATH` environment variable.

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

## 📝 Citing AddTag ##
If you use AddTag for your research, please cite us. Because the manuscript is currently in preparation, you will need to cite the code repository instead.

 > Thaddeus D. Seher and Aaron D. Hernday. Addtag: Program for identifying exclusive endogenous gRNA sites and creating unique synthetic gRNA sites. University of California, Merced. Retrieved from \<https://github.com/tdseher/addtag-project> (2019).

## ✍ Authors ##
Who do I talk to?
 * Aaron D. Hernday (🔬 PI leading the project)
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

 * `SingleSequenceAlgorithm` objects calculate scores by comparing a potential RNA ![Spacer][Spacer] or DNA ![Target][Target] to a model trained on empirical data.
 * `PairedSequenceAlgorithm` instances generate scores that compare a potential RNA ![Spacer][Spacer] to a DNA ![Target][Target].
 
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

## Limitations ##
 * AddTag makes no effort to restrict which Target motifs the user can use according to the selected Algorithms. Therefore, the user needs to independently verify which Target motifs are compatible with the selected Algorithms.
 * Right now AddTag can only handle linear chromosomes. If you want to analyze a circular chromosome, then you will need to artificially concatenate the ends of the chromosome together and adjust any annotations before running AddTag. An additional complication the software does not address is circular chromosomes. Features and their flanking regions cannot span the junction created when the contig end is concatenated to the start (typically the starting position on a contig is labeled the ORIGIN). To address this, the user should manually shift the coordinates of the experimental Features, and wrap the contigs as appropriate.
 * AddTag assumes one Feature copy per contig. The current implementation of AddTag assumes homology regions around Features are not repeated across any one contig. This means that is will fail to generate cPCR oligos for a large proportion of genes in transposon-rich genomes such as wheat (https://dx.doi.org/10.1186/s13059-018-1479-0).
 * A single feature cannot span two or more contigs. AddTag assumes that the entire feature sequence, and any flanking regions, are not in terminal regions of the reference contig. 
 * AddTag does not address overlapping genes, such as when an intron contains an exon for another gene, or when the same DNA encodes for genes on opposite strands. Everything between the Feature bounds is removed in the first engineering step. Currently, if the selected Feature overlaps with any other feature, only the selected Feature is considered. The other Feature will be disrupted. AddTag will report a warning that these other Features may be disrupted, but it does not attempt to reconcile this in any way. However, AddTag does have the ability to limit Feature expansion to keep the deletion outside of neighboring Features.
 * AddTag was not designed to perform paired Cas design, such as FokI-dCas9 nickase Users would need to run the program and select two gRNAs designed for opposite strands within a certain distance from each other. One way to mitigate errors is to use PAM-out nickases. This requires Cas9 cutting by two targets to get double-stranded break. This significantly decreases off-target genome editing. However, this initial AddTag version does not explicitly facilitate this.
 * AddTag can identify cut sites for Cas enzymes which have the PAM site. No functionality is provided for finding sites without an adjacent PAM sequence. AddTag requires motifs to define a PAM sequence. Therefore Cas14a is not supported. The number of CRISPR/Cas genome editing technologies are rapidly growing. With the recent discovery of Cas14a, which targets single-stranded DNA (ssDNA) molecules without requiring a PAM site (https://dx.doi.org/10.1126/science.aav4294), the expanded prevalence of CRISPR/Cas methods in biological sciences is assured. However, often researchers wish to edit sites on double-stranded DNA (dsDNA) using an RGN (such as Cas9 or Cas12a) that requires binding to a PAM motif. 
 * Please note, that at this time, no special restriction sites will be taken into account when designing primers.
 * For simplicity, all calculated scores ignore terms dealing with proximity to exon/CDS/ORF sequences. In cases such as the Stemmer and Azimuth calculations, the authors attempted to include the risk of disrupting genes neighboring potential targets in their models. We don’t attempt to do this.
 * Additionally, some scoring Algorithms take chromatin structure (DNA accessibility) into account. For simplicity, AddTag treats all input gDNA as equally accessible.
 * During the course of writing this software, a paper was published that outlines how hairpins can be inserted into the pre-spacer and spacer regions of the gRNA in order to increase specificity (https://dx.doi.org/10.1038/s41587-019-0095-1). AddTag does not model pre-spacer sequences.
 * AddTag assumes the RGN template type is dsDNA. AddTag was designed specifically to enable efficient gDNA editing. It does not use predictive models for ssDNA or RNA templates.
 * Since Bartag motifs are user-specified, simple pre-computed lists of compatible 'bartag' sequences would be incomplete. Thus we implemented a greedy 'bartag' generation algorithm. When evaluating candidate 'bartag' sequences, AddTag will keep 'bartags' that satisfy all edit distance requirements with all previously-accepted 'bartags'. To limit runtime to a reasonable amount, we limited the total number of Features and 'bartags' that can be generated.
 * Of special note are things the Primer design does not explicitly consider, such as characteristics of the cPCR template molecule. AddTag does not exploit the differential nature of template sequence composition (e.g. H. sapiens compared to E. coli). Also, AddTag does not use information on the presence of known secondary modifications to the template, such as methylated residues or oxidative damage.

[tdseher]:https://twitter.com/tdseher
[Spacer]:docs/spacer.svg
[PAM]:docs/pam.svg
[Target]:docs/target.svg
[dDNA]:docs/ddna.svg
[RGN]:docs/rgn.svg
[Feature]:docs/feature.svg
[Primer]:docs/primer.svg
[BLAST+ source]:ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST
[BLAST+ binaries]:ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST
