# CRISPR/Cas AddTag Readme #

Program for identifying exclusive endogenous gRNA sites and creating unique synthetic gRNA sites.

[![Linux](https://img.shields.io/badge/Linux-✓-darkgreen.svg?logo=linux)](#)
[![Windows](https://img.shields.io/badge/Windows-✓-darkgreen.svg?logo=windows)](#)
[![macOS](https://img.shields.io/badge/macOS-✓-darkgreen.svg?logo=apple)](#)
[![Build Status](https://dev.azure.com/tdseher/addtag-project/_apis/build/status/tdseher.addtag-project?branchName=master)](https://dev.azure.com/tdseher/addtag-project/_build/latest?definitionId=1&branchName=master)
[![Code coverage](https://img.shields.io/azure-devops/coverage/tdseher/addtag-project/1)](#)

[![Python](https://img.shields.io/badge/Python-≥3.5.1-1f425f.svg?logo=python)](https://www.python.org/downloads/release/python-360/)
[![PRs](https://img.shields.io/badge/PRs-welcome-brightgreen.svg)](http://makeapullrequest.com)
[![](https://img.shields.io/badge/doi-...-blue.svg)](#)

[![Twitter](https://img.shields.io/badge/Twitter-@tdseher-blue?logo=twitter)](https://twitter.com/tdseher)
[![Open Source](https://img.shields.io/badge/Open%20Source-❤-teal)](#)
[![free software](https://img.shields.io/badge/Software-%F0%9F%97%BD%20free-teal)](https://www.gnu.org/philosophy/free-sw.html)

[Features](#-features) • [Requirements](#-requirements) • [Installing](#-installing-addtag) • [Usage](#-program-instructions) • [Aligners](#-supported-sequence-aligners) • [Thermodynamics](#-supported-thermodynamics-calculators) • [Algorithms](#-supported-scoring-algorithms) • [Citing](#-citing-addtag) • [Contributing](#-contributing)

## ☑ Features ##
Basic Features:
 * [x] Analyzes any arbitrary genomic DNA (gDNA).
   * [x] Fully supports ambiguous characters or polymorphisms (`RYMKWSBDHVN`).
   * [x] Respects case-masked gDNA for [![Target][Target]](#) and [![Primer][Primer]](#) identification.
 * [x] Locates RNA-guided nuclease ([![RGN][RGN]](#)) cut sites ([![Target][Target]](#)s) within a [![Feature][Feature]](#) (locus of interest) for optimal gRNA [![Spacer][Spacer]](#)s.
   * [x] Fully supports ambiguous bases (`RYMKWSBDHVN`) in [![Spacer][Spacer]](#) or [![PAM][PAM]](#).
   * [x] Accepts 3'-adjacent [![PAM][PAM]](#) sequences, such as Cas9 (`>NGG`).
   * [x] Accepts 5'-adjacent [![PAM][PAM]](#) sequences, such as Cas12a (`TTTN<`).
   * [x] Supports arbitrary [![Spacer][Spacer]](#) length and composition constraints, such as for plant experiments (`G{,2}N{19,20}`).
   * [x] Supports arbitrary [![PAM][PAM]](#) sequences (MAD7: `YTTN<`, Cas12d: `TA<`, BlCas9: `>NGGNCNDD`, etc).
   * [x] Uses stranded forward (`/`), reverse (`\`) and unstranded (`|`) cut sites.
   * [x] Supports [![PAM][PAM]](#) sequences defined by complex nested logic, such as xCas9 (`>(N{1,2}G,GAW,CAA)`)
   * [x] Calculates any number of **on-target** and **off-target** scores (see [Algorithms](#-supported-scoring-algorithms)).
   * [x] Finds homology-aware [![Target][Target]](#)s (**multi-allelic**, **allele-specific**, and **allele-agnostic**).
   * [x] Searches for [![Target][Target]](#)s using selectable pairwise alignment program (see [Aligners](#-supported-sequence-aligners)).
 * [x] Generates exogenous, donor DNAs ([![dDNA][dDNA]](#)s) to modify the same locus successively.
   * [x] Assembles unique [![Target][Target]](#) sites (on [![dDNA][dDNA]](#)s), thus maximizing **on-target** and **off-target** scores (because they don't resemble any input gDNA).
   * [x] Adds unique [![Target][Target]](#)s to [![dDNA][dDNA]](#)s without inserting sequence (or while introducing minimal amounts of extrinsic DNA).
   * [x] Engineers a single set of conservative PCR (cPCR) [![Primer][Primer]](#)s that work for all genotypes (wild type, knock-out, and add-back) to validate if a [![Feature][Feature]](#) was engineered correctly.
   * [x] Produces homology-aware [![dDNA][dDNA]](#)s (**multi-allelic**, **allele-specific**, and **allele-agnostic**).
 * [x] Performs *in silico* recombination between gDNA and [![dDNA][dDNA]](#)s.
 * [x] Determines thermodynamic properties of sets of [![Primer][Primer]](#) pairs (Tm, minimum ΔG, amplicon size, etc).
 * [x] Displays all known [![RGN][RGN]](#) [![Spacer][Spacer]](#) and [![PAM][PAM]](#) combinations.

## 📋 Requirements ##

### Hardware recommendations ###

Processor:
 * ≥ 4 cores, ≥ 3 GHz

Computations scale fairly linearly, so the more computational cores you can assign to the task, the faster it will go.

Memory:
 * ≥ 4 Gb (for [![Target][Target]](#) evaluation)
 * ≥ 4 Gb (for [![Primer][Primer]](#) evaluation)

See [Notes](#notes) for tips on memory optimization.

### Software requirements ###

Below are lists AddTag requirements. Each entry is marked with a ☑ or ☐, indicating whether or not an additional download/setup is required:

 * [x] All requirements included in AddTag
 * [ ] Additional download/setup required

For tips on setting up AddTag requirements, please review the commands in the `.azure-pipelines.yml` file.

#### Basic prerequisites ####

Base operation of AddTag requires the following:

 * [ ] Python ≥ 3.5.1 ([source](https://www.python.org/downloads/), [binaries](https://www.python.org/downloads/), [documentation](https://docs.python.org/3/))

 * [ ] regex Python module ([source](https://bitbucket.org/mrabarnett/mrab-regex), [whls](https://pypi.org/project/regex/), [documentation](https://pypi.org/project/regex/))

Certain optional AddTag functionality (version information, and software updates) depends on the following:

 * [ ] Git ≥ 1.7.1 ([source](https://github.com/git/git), [binaries](https://git-scm.com/downloads), [documentation](https://git-scm.com/doc))

#### 📐 Supported sequence Aligners ####

One pairwise sequence aligner is required:

 * [ ] BLAST+ ≥ 2.6.0 ([source](https://bit.ly/2Ouoqkx), [binaries](https://bit.ly/2Ouoqkx), [documentation](https://www.ncbi.nlm.nih.gov/books/NBK279690/))

 * [ ] Bowtie 2 ≥ 2.3.4.1 ([source](https://github.com/BenLangmead/bowtie2), [binaries](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/), [documentation](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml))

 * [ ] BWA ≥ 0.7.12 ([source](https://github.com/lh3/bwa), [ugene binaries](http://ugene.net/download-all.html#en_data_analysis_tools), [bioconda binaries](https://anaconda.org/bioconda/bwa/files), [documentation](http://bio-bwa.sourceforge.net/bwa.shtml))

 * [ ] Cas-OFFinder ≥ 2.4 ([source](https://github.com/snugel/cas-offinder), [binaries](https://sourceforge.net/projects/cas-offinder/files/Binaries/), [documentation](http://www.rgenome.net/cas-offinder/portable))

<!--
##### The following third-party integrations are currently incomplete #####

 * [ ] BLAT ([source](https://genome.ucsc.edu/goldenPath/help/blatSpec.html), [binaries](http://hgdownload.cse.ucsc.edu/admin/exe/), [documentation](https://genome.ucsc.edu/goldenPath/help/blatSpec.html))

 * [ ] Bowtie ([source](https://github.com/BenLangmead/bowtie), [binaries](https://sourceforge.net/projects/bowtie-bio/files/bowtie/), [documentation](http://bowtie-bio.sourceforge.net/manual.shtml))

 * [ ] Usearch ([binaries](https://www.drive5.com/usearch/download.html), [documentation](https://drive5.com/usearch/manual/cmds_all.html))
-->

For polymorphism-aware [![Feature][Feature]](#) expansion (using the `--homologs` option), one multiple sequence aligner is required:

 * [ ] MAFFT ([source](https://mafft.cbrc.jp/alignment/software/source.html), [binaries](https://mafft.cbrc.jp/alignment/software/), [documentation](https://mafft.cbrc.jp/alignment/software/manual/manual.html))

#### 🌡 Supported thermodynamics calculators ####

For oligo design, AddTag requires one of the following third-party thermodynamics solutions to be installed:

 * [ ] UNAFold ≥ 3.8 ([source](http://rnaspace.sourceforge.net/software/unafold-3.8.tar.gz), [documentation](http://unafold.rna.albany.edu/)) with [patch440](http://unafold.rna.albany.edu/?q=node/440)

 * [ ] primer3-py Python module ([source](https://github.com/libnano/primer3-py), [whls](https://pypi.org/project/primer3-py/), [documentation](https://libnano.github.io/primer3-py/))

 * [ ] ViennaRNA Python module ([source](https://github.com/ViennaRNA/ViennaRNA), [official binaries](https://www.tbi.univie.ac.at/RNA/index.html#download), [bioconda binaries](https://anaconda.org/bioconda/viennarna/files), [documentation](https://www.tbi.univie.ac.at/RNA/documentation.html))

#### 📈 Supported scoring Algorithms ####

The following scoring algorithms are subclasses of `SingleSequenceAlgorithm`.

 * [ ] Azimuth ([Doench, Fusi, et al (2016)](http://dx.doi.org/10.1038/nbt.3437))

   note: Either Azimuth 2 or Azimuth 3 can be used to calculate Azimuth scores. There is no need to have both installed.
   
    * Azimuth 3 Python module ([source](https://github.com/milescsmith/Azimuth), [documentation](https://www.microsoft.com/en-us/research/project/crispr/))
      
      note: requires specific versions of numpy, scikit-learn, and pandas.
      Other dependencies include click, biopython, scipy, GPy, hyperopt, paramz, theanets, glmnet_py, dill, matplotlib, pytz, python-dateutil, six, tqdm, future, networkx, pymongo, decorator, downhill, theano, nose-parameterized, joblib, kiwisolver, cycler, pyparsing, setuptools, glmnet-py.

    * Azimuth 2 Python module ([source](https://github.com/MicrosoftResearch/Azimuth), [documentation](https://www.microsoft.com/en-us/research/project/crispr/))
      on 2.7.10 ≤ Python < 3.0.0 ([source](https://www.python.org/downloads/), [binaries](https://www.python.org/downloads/), [documentation](https://docs.python.org/2/)) 
      
      note: requires python-tk to be installed. Also requires specific versions of scipy, numpy, matplotlib, nose, scikit-learn, pandas, biopython, pyparsing, cycler, six, pytz, python-dateutil, functools32, subprocess32.

 * [ ] CINDEL/DeepCpf1 ([Kim, Song, et al (2016)](http://dx.doi.org/10.1038/nmeth.4104), [Kim, Song, et al (2018)](https://doi.org/10.1038/nbt.4061))
       
   note: Requires both Keras and Theano Python modules.
   
    * Keras Python module ([source](https://github.com/keras-team/keras), [whls](https://pypi.org/project/Keras/), [documentation](https://keras.io/))
    * Theano Python module ([source](https://github.com/Theano/Theano), [whls](https://pypi.org/project/Theano/), [documentation](http://deeplearning.net/software/theano/))

 * [x] Doench-2014 ([Doench, et al (2014)](http://dx.doi.org/10.1038/nbt.3026))

 * [x] Housden ([Housden, et al (2015)](http://dx.doi.org/10.1126/scisignal.aab3729))

 * [x] Moreno-Mateos ([Moreno-Mateos, et al (2015)](http://dx.doi.org/10.1038/nmeth.3543))

 * [x] CRISPRater ([Labuhn, et al. (2018)](http://dx.doi.org/10.1093/nar/gkx1268))

 * [x] GC ([Wang, et al (2014)](http://dx.doi.org/10.1126/science.1246981'))

 * [x] Homopolymer ([Hough, et al. (2017)](https://doi.org/10.1186/s12859-017-1581-4))

 * [x] ProximalG

 * [x] PolyT

 * [x] PAM Identity

 * [x] Position

The following scoring algorithms are subclasses of `PairedSequenceAlgorithm`.

 * [x] CFD ([Doench, Fusi, et al (2016)](http://dx.doi.org/10.1038/nbt.3437))

 * [x] Substitutions, Insertions, Deletions, Errors ([Needleman, Wunsch (1970)](https://dx.doi.org/10.1016/0022-2836%2870%2990057-4))

 * [x] Hsu-Zhang ([Hsu, et al (2013)](http://dx.doi.org/10.1038/nbt.2647))

 * [x] Linear

#### Python package setup ####
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

To make the AddTag executable accessible from any working directory, you can add the absolute path of the current working directory to the `PATH` variable.

On Windows, run:
```sh
set PATH=%PATH%;%CD%
```

On Linux or macOS, run:
```sh
export PATH=$PATH:$PWD
```

If you run AddTag with no parameters, you should get the following output:
```
usage: addtag [-h] [-v] action ...
```

#### Special note ####
One way to obtain AddTag is by downloading and extracting the code directly from GitHub:
```sh
wget https://github.com/tdseher/addtag-project/archive/master.zip
unzip master.zip
cd addtag-project-master/
```

If you try running `addtag`, you will get a message similar to the following:
```sh
./addtag
```
> ```
> fatal: Not a git repository (or any parent up to mount point /media/sf_VirtualBox_share)
> Stopping at filesystem boundary (GIT_DISCOVERY_ACROSS_FILESYSTEM not set).
> ```
This message means that the AddTag directory isn't a valid `git` repository (it is missing the `.git` subfolder).
As a consequence, the version information will not be accessible.
```sh
./addtag --version
```
> ```
> addtag missing (revision missing)
> ```
To fix this, simply ensure `git` is installed and available in the `PATH` environment variable
(See [Software prerequisites](#basic-prerequisites)), and run the following:
```sh
./addtag update
```

Now, when you run `addtag`, you should not receive the warnings, and the version field will be populated.
```sh
./addtag --version
```
> ```
> addtag 9e8748b (revision 460)
> ```

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

## 💻 Program Instructions ##
### Displaying the usage ###
<details>
<summary>Click to expand/collapse</summary>
<table><tbody><tr><td>

Because AddTag is being updated regularly, the most current feature set and usage can be viewed by running AddTag with the `--help` command line option.

The following commands assume the current working directory is the AddTag folder `addtag-project/`. This will print out command line parameter descriptions and examples.
```sh
./addtag --help
```

Additionally, you may view the included man page, which is probably not up-to-date.
```sh
man ./addtag.1
```

</td></tr></tbody></table>
</details>

### Format of input data ###
<details>
<summary>Click to expand/collapse</summary>
<table><tbody><tr><td>

#### FASTA input ####
AddTag requires a FASTA genome of the organism you wish to manipulate. FASTA files resemble the following:
> ```fasta
> >primary_header1 attribute1=value1 attribute2=value2
> NNNNCGAAATCGGCGCATAGGCCTAAGAGCTCCTATAGAGATCGATATAAAC
> GCTAGAAGATAGAGAGAGGCTCGCGCTCGATCGCGATAAGAGAGCTCTCGGC
> CGATAGAGGAATCTCGggctcgcatatatyhcgcggcatatGGCCTAGAGGA
> CCAATAAAGATATATAGCCTAAAGGAATATATAGAGAGATATATATAGNNNN
> >primary_header2 attribute1=value1 attribute2=value2
> AGCTAGAGACWWWCTCCTCTCCTAGAGASSSAGAGGAGAGCTCTCCGAGAGA
> CGCTCGCTCGTATGCCTCTATATCGATATATAGGAGAATCCTCGATATATAG
> ```
FASTA files are plain text files that use newline (`\n` or `\r\n`) characters as delimiters. If a line begins with a greater than (`>`) symbol, it represents the start of a new sequence record. All characters between the `>` and `\n` are considered the 'header' of the record. Everything between the `>` and the first whitespace character (` ` or `\t`), if one exists, is considered the 'primary identifier' for the record. All subsequent lines until the next 'header' line contain the sequence information for that record. Therefore FASTA files can contain many sequence records. Each record in a genome assembly's FASTA file is called a 'contig'. 

Typically, the DNA sequence information in FASTA files are list a bunch of canonical nucletide abbreviations (`ACGT`). However, FASTA files can contain any number of ambiguous characters (`RYMKWSBDHVN`), which can represent allelic variation expected within the sample or sequencing uncertainty. FASTA files can also contain a mix of `UPPER` and `lower` cased characters. Typical use for `lower` case characters is to exclude these residues from [![Target][Target]](#) or [![Primer][Primer]](#) identification. 

#### GFF input ####
AddTag requires a GFF file containing annotations for the Features you wish to manipulate ([technical specifications of GFF format](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)). GFF files resemble the following:
> ```gff
> # seqid	source	feature	start	end	score	strand	frame	attribute
> C1A	DB	gene	3489	5146	.	+	.	ID=C1A_001;Name=C1_001;Gene=GENE1
> C1A	DB	mRNA	3489	5146	.	+	.	ID=C1A_001-T;Parent=C1A_001
> C1A	DB	exon	3489	5146	.	+	.	ID=C1A_001-T-E1;Parent=C1A_001-T
> C1A	DB	CDS	3489	5146	.	+	0	ID=C1A_001-P;Parent=C1A_001-T
> 
> C1B	DB	gene	3267	4924	.	+	.	ID=C1B_001
> C1B	DB	mRNA	3267	4924	.	+	.	ID=C1B_001-T;Parent=C1B_001
> C1B	DB	exon	3267	4924	.	+	.	ID=C1B_001-T-E1;Parent=C1B_001-T
> C1B	DB	CDS	3267	4924	.	+	0	ID=C1B_001-P;Parent=C1B_001-T
> ```
GFF files describe the contig locations of important genomic Features. Empty lines and lines that begin with the pound (`#`) symbol are ignored. Of note is the far-right `attribute` column, which AddTag assumes is a semicolon-delimited set of key/value pairs. AddTag assumes each Feature has a unique identifier. By default, it uses the `ID` attribute as the unique name for each Feature. If your GFF file does not have an `ID` attribute, then you can select a different one with the `--tag` command line option. 

Typical AddTag analyses require at least one GFF file. AddTag can handle GFF files in two ways.
 * For the first method, all Features matching the selected type, designated by the `--features` command line argument, will be included for analysis. By default, only lines in the GFF file containing `gene` in `feature` column will be considered. This system is useful if your GFF file contains only the Features you wish to manipulate.
 * If your GFF file contains all annotations for the entire genome (which is typical), the second approach requires you to select only the few Features you want to edit using the `--selection` command line argument.
 
Often, you will have a GFF file with annotations for the entire genome. The `attributes` column is not often structured intuitively, and can prove cumbersome to search (`grep`) or sort (`sort`) manually. To make it easy to identify the desired lines of a GFF file, AddTag includes the `find_feature` subroutine. Here is an example that tries to find all lines associated with `HSP90` by searching several attribute tags, and outputting a GFF with a commented line containing field names:
```sh
addtag find_feature --gff genome.fasta --query HSP90 --linked_tags Name Alias Parent Gene --header > features.gff
```

#### Target motif input  ####
The Target motif is written from 5' to 3'. Use a greater than (`>`) symbol if your [![RGN][RGN]](#) has a 3'-adjacent PAM, and use a less than (`<`) symbol if your [![RGN][RGN]](#) has a 5'-adjacent PAM. Ambiguous nucleotide characters are accepted. `{a,b}` are quantifiers. `(a,b,…)` are permitted alternatives. `/` is a sense strand cut, `\` is an antisense strand cut, and `|` is a double-strand cut. `.` is a base used for positional information, but not enzymatic recognition. Be sure to enclose each motif in quotes so your shell does not interpret `STDIN`/`STDOUT` redirection.

You can specify any number of Target motifs to be considered 'on-target' using the `--motifs` command line option. You can also designate any number of Target motifs to be considered 'off-target' using the `--off_target_motigs` command line option. 

To see an exhaustive list of all identified Target motifs for each known [![RGN][RGN]](#), run the following command:
```sh
addtag list_motifs
```

#### Homologs input ####
Some researchers are lucky enough to get to work on organisms with phased genomes. This means that full haplotype information is known for each chromosome. AddTag can accommodate haploid, diploid, and polyploid genomes when homologous Features are linked by the addition of the `--homologs` command line option. The 'homologs' file has the following format:
> ```homologs
> # group	hom_a	hom_b	hom_c
> GENE1	C1A_001	C1B_001
> GENE2	C1A_002	C1B_002	C1C_002
> ```
Each Feature identifier has its contig start and end position defined in the input GFF file. The 'homologs' file merely links them together. Columns in the homologs file are delimited by the `\t` character. The first column is the name of the group of Features. Every subsequent column should contain the identifier of a Feature to consider as a homolog. Homolog groups can each have any number of Features. If a Feature identifier appears on multiple lines, then all those Features are linked together as one homolog group. The identifier can be changed with the `--tag` command line option.  

</td></tr></tbody></table>
</details>

### Format of output data ###
<details>
<summary>Click to expand/collapse</summary>
<table><tbody><tr><td>

AddTag outputs most of the experimental results you need to `STDOUT`. However, for simplicity sequences are output to `FASTA` files. Please note that the output table formats are not consistent among AddTag versions--more recent releases are more thorough and useful.

#### STDOUT ####
The final data are printed to `STDOUT` as tab-delimited tables. Lines containing column headers start with a `#` character.

The `reTarget results` table contains information on optimal Targets that exist within the ★tag insert on the r1-gDNA.

The `exTarget results` table contains information on optimal Targets that exist within the extended Feature on the r0-gDNA.

The `AmpF/AmpR results` table contains information on optimal Primer Pairs for amplifying the Feature to create the r2-dDNAs.

The `reDonor results` table lists information on r2-dDNAs.

The `exDonor results` table lists information on r1-dDNAs.

The `Region definitions` table lists the genome, contig, start and end coordinates for where cPCR Primers will be selected from.

The `Primer sequences` table lists the optimal PrimerSets by weight order, with non-redundant Primers

The `PrimerPairs` table lists the PrimerPair attributes for each amplicon.

The `Amplicon diagram` succinctly relates the primer names to the regions in the genomes they bind to and amplify. 

The `In silico recombination` table lists where in the gDNAs the dDNAs were incorporated.

#### STDERR ####
If the AddTag software fails for any reason, error messages will be printed to `STDERR`. If you pipe `STDERR` into a file, and the file size is nonzero, then this indicates that an error occurred.

Often, errors happen if required AddTag arguments are missing, or input data is improperly formatted. 

#### log.txt ####
AddTag outputs intermediate calculations and computation status to the `log.txt` file. This includes the exact commands used when calling any external programs (such as [Aligners](#-supported-sequence-aligners)), alignments of Target sequences to dDNA sequences, and timestamps.

#### excision-dDNAs.fasta ####
The `excision-dDNAs.fasta` file contains the dDNA sequences for creating the intermediary genome that are referenced by the tables from `STDOUT`. These dDNA sequences contain the `mintag`, `addtag`, `unitag`, `bartag`, or `sigtag` as requested by the AddTag invocation arguments.

An example of a nominal `mintag` that targets both alleles of a diploid chromosome:
> ```fasta
> >exDonor-0 spacers=4 C1A_002:C1A:+:272323..272373::274197..274247 C1B_002:C1B:+:272338..272388::274212..274262
> ACTAAAATGAAAACCACATACAGCAGTAATAGTACTAGCCAACTCACTATTTTGATTTTGGGAACGGAGTTGAGCGGTATATGTGACAACAGTGACTATG
> ```

An example of an `addtag` experiment:
> ```fasta
> >exDonor-0 spacers=1 C1A_003:C1A:+:109972..110010:ctccgctctcgcctagactcggg:112195..112234 C1B_003:C1B:+:109967..110005:ctccgctctcgcctagactcggg:112220..112259
> GCATAGGCTAGAGATAGTCCTCAGATAATAATAGAGCTctccgctctcgcctagactcgggAATATAAGATCAGTCTCTCCCGACTAGAATCTCTAGCAA
> >exDonor-1 spacers=1 C1A_003:C1A:+:109972..110010:cccgagtctaggcgagagcggag:112195..112234 C1B_003:C1B:+:109967..110005:cccgagtctaggcgagagcggag:112220..112259
> GCATAGGCTAGAGATAGTCCTCAGATAATAATAGAGCTcccgagtctaggcgagagcggagAATATAAGATCAGTCTCTCCCGACTAGAATCTCTAGCAA
> ```
These dDNAs each are predicted to recombine with contigs `C1A` and `C1B`. Note that each dDNA incorporates the exogenous `addtag` sequence in an opposite orientation.

#### excision-targets.fasta ####
This file contains only the Target sequences that are contained within the Feature, but in `FASTA` format. For the most part, the `exTarget results` table from `STDOUT` contains more information. We intend this file to be used as input to the `find_header` subroutine.

#### reversion-dDNAs.fasta ####
This file is structured identically to the `excision-dDNAs.fasta` file.

If you direct AddTag to find Primers to amplify the wild type Feature, then their amplicon sequences will be stored in the `reversion-dDNAs.fasta` file. If you do not have AddTag find the AmpF/AmpR primers, then the entire region containing the Feature, upstream, and downstream sequences is written to the `reversion-dDNAs.fasta` file. 

This example shows that polymorphisms at the Feature and its flanking sequences mean there are two possible dDNAs:
> ```fasta
> >reDonor-0 spacers=0 C3A_005:C3A:+:1722491..1722834
> TTTTTTTTGGTTAACCACTTTGTGTCCCTTGCATACTTTTACATTGGAAACATACATACACTAACATTCACACTCAATAC
> ACTCATATTATTTACCATTTTTGTTGTGAAGATACACGTATTTATTGAGTATTCCTTCATAACATTTAATTTATATTCCA
> AGAGTTAATTGATTAAACAACTTGGTCCAAACAAACATAAACATAAACAAAAACGTTTTCTTTTTTTGCATAATATCTAT
> CTATGTATATGTATATATATGTGTGTAAGTCATTGTCTTTTCCATTTTCTTTTCCATTTTCTTTTTTTTTTAGTTTTGTT
> TTCAAGTGTGTAATAATAATAAT
> >reDonor-1 spacers=0 C3B_005:C3B:+:1723088..1723418
> TTTTTTTTGGTTAACCCCTTTGTGTCCCTTGCATACTTTTACATTGGAAACATACATACACTAACATTCACACTCAATAC
> ACTCACATTATTTACCATTTTTGTTGTGAAGATACACGTATTTATTGAGTATTCCTTCATAACATTTAATTTATATTCCA
> AGAGTTAATTGATTAAACAACTTGGTCCAAAAAACAAAAACGTTTTCTTTTTTTGCATAATATCTATCTATGTATATGTA
> TATATATGTGTGTAAGTCATTGTCTTTTCCATTTTCTTTTCCATTTTCTTTTCTTTTTAGTTTTGTTTTCAAGTGTGTAA
> TAATAATAAT
> ```

#### reversion-targets.fasta ####
This file contains only the RGN Target sequences compatible with the `exDonor` sequences (and by extension, the intermediary genome). For the most part, the `reTarget results` table from `STDOUT` contains more information.

#### genome-rN.fasta ####
In silico recombination will integrate the input dDNAs into their respective loci within the input genome. Contig names (primary identifiers) are modified with the incorporated dDNAs as well as the round.

For example, `genome-r0.fasta` may resemble the following:
> ```fasta
> >contig_001
> GCTAAGCGCATCGCGCATAGGGCGGCAAAAAAGCGCTAGAGACTCAGAGGAGCGCTAGCG
> GCTCGAATATAATAGATAGCTATAGCCTAGGAGATAGGAAACTCAGAAATAGACCATAAA
> >contig_002
> AATAAGCTCAGATAATATAGCTCGCTCTCTCGATAGCTCTAGACTCCCTAGAGCCCTAAG
> CCCGCTCGCGAATAGATCCTCTAGACTAGATGAGAGCCGGCCCTCGCGCGCGATAGAGAA
> ```

If the first round dDNA contains the following:
> ```fasta
> >dDNA1
> GCTCGAATATAATAGATAGCTATAGcccgggAGGAAACTCAGAAATAGACCATAAA
> ```

After the first round of *in silico* recombination, `genome-r1.fasta` will be:
> ```fasta
> >contig_001-r1[dDNA1]
> GCTAAGCGCATCGCGCATAGGGCGGCAAAAAAGCGCTAGAGACTCAGAGGAGCGCTAGCG
> GCTCGAATATAATAGATAGCTATAGcccgggAGGAAACTCAGAAATAGACCATAAA
> >contig_002
> AATAAGCTCAGATAATATAGCTCGCTCTCTCGATAGCTCTAGACTCCCTAGAGCCCTAAG
> CCCGCTCGCGAATAGATCCTCTAGACTAGATGAGAGCCGGCCCTCGCGCGCGATAGAGAA
> ```

</td></tr></tbody></table>
</details>

### Available subroutines ###
<details>
<summary>Click to expand/collapse</summary>
<table><tbody><tr><td>

The AddTag program contains a set of subroutines that can be run independently. There are four categories of subroutines.

 * The `evaluate_*` subroutines run only a very specific analysis on input data.
 * The `find_*` subroutines are used to search input files for specific things, so the user can easily learn the correct parameters to use for AddTag input.
 * The `generate_*` subroutines perform the deep computational analyses.
 * The `list_*` subroutines just print information the user might find useful.

</td></tr></tbody></table>
</details>

### Available RGN scoring Algorithms ###
<details>
<summary>Click to expand/collapse</summary>
<table><tbody><tr><td>

Over the past few years, several Algorithms have been proposed to describe [![RGN][RGN]](#) behavior within certain biological contexts. We implemented most of the commonly-used ones into the AddTag software. To view information about each, use the following command:
```sh
addtag list_algorithms
``` 
This will write the pertinent information for all implemented Algorithms to `STDOUT`.

If an Algorithm is used for pre-alignment filtering (`Prefilter`) or post-alignment filtering (`Postfilter`), then the score of the Target must lie between the `Min` and `Max` values to be continued on through the analysis. For instance, the 'off-target' scoring `CFD` Algorithm has a `Min` of `1.0`. This means that some positions with significant sequence similarity to the query Target (because they are identified in the Alignment step) will not contribute to the final 'off-target' score if their score is less than `1.0`.

</td></tr></tbody></table>
</details>

### Available oligonucleotide thermodynamics calculators ###
<details>
<summary>Click to expand/collapse</summary>
<table><tbody><tr><td>

To view which thermodynamics calculators are available on your system, use the following command:
```sh
addtag list_thermodynamics
```

</td></tr></tbody></table>
</details>

### Workflow for editing loci in the manuscript ###

#### ADE2_CDS ####
<details>
<summary>Click to expand/collapse</summary>
<table><tbody><tr><td>

*~ Section incomplete ~*

</td></tr></tbody></table>
</details>

#### EFG1_CDS ####
<details>
<summary>Click to expand/collapse</summary>
<table><tbody><tr><td>

*~ Section incomplete ~*

</td></tr></tbody></table>
</details>

#### BRG1_CDS ####
<details>
<summary>Click to expand/collapse</summary>
<table><tbody><tr><td>

*~ Section incomplete ~*

</td></tr></tbody></table>
</details>

#### ZAP1_US ####
<details>
<summary>Click to expand/collapse</summary>
<table><tbody><tr><td>

*~ Section incomplete ~*

</td></tr></tbody></table>
</details>

#### ZRT2_US ####
<details>
<summary>Click to expand/collapse</summary>
<table><tbody><tr><td>

*~ Section incomplete ~*

</td></tr></tbody></table>
</details>


### Typical workflows ###

#### 1-step deletion of a single Feature ####
<details>
<summary>Click to expand/collapse</summary>
<table><tbody><tr><td>

In this simplest of examples, we will choose a Feature to delete from a genome, identify the optimal Target to design the gRNA against, create the necessary dDNA, and generate the set of Primers to validate the deletion.

This process uses a 'nominal' `mintag`, which means the generated dDNA consists of homology arms concatenated together with no insert.

The first step is to obtain input data. Let's download the sequences (FASTA) and annotations (GFF) for a haploid *C. albicans* assembly into the current working directory:
```sh
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/Candida_albicans/all_assembly_versions/GCF_000182965.3_ASM18296v3/GCF_000182965.3_ASM18296v3_genomic.fna.gz
gunzip GCF_000182965.3_ASM18296v3_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/Candida_albicans/all_assembly_versions/GCF_000182965.3_ASM18296v3/GCF_000182965.3_ASM18296v3_genomic.gff.gz
gunzip GCF_000182965.3_ASM18296v3_genomic.gff.gz
```

For convenience, let's use a variable to abbreviate these paths:
```sh
GENOME=GCF_000182965.3_ASM18296v3_genomic
```

The 1-step approach is appropriate when the Feature you wish to remove contains a high quality Target within it. We will select a Feature from the GFF file using the `--selection` option.

Let's pretend we are interested in the gene `GCN20`. Let's store it into a variable.
```sh
GENE=GCN20
```
For the purposes of this walkthrough, `GCN20` is interesting because all its potential Cas9 Targets have several off-targets across the genome.
Because there is no precise Target, the Algorithm weight is especially useful for balancing the on-target and off-target scores.

If we know its gene ID, we can directly include the option
`--selection ID`. However, we don't know the ID for this gene, so we can search for it. To do this, we will use
the `addtag find_feature` subroutine to find all Features associated with `GCN20`:
```sh
addtag find_feature --linked_tags --header --query ${GENE} --gff ${GENOME}.gff
```
> ```
> # seqid	source	feature	start	end	score	strand	frame	attribute
> NC_032089.1	RefSeq	CDS	75573	77828	.	-	0	ID=cds-XP_719022.1;Parent=rna-XM_713929.2;Dbxref=CGD:CAL0000181616,GeneID:3639314,Genbank:XP_719022.1;Name=XP_719022.1;Note=YEF3-subfamily ABC family protein%2C predicted not to be a transporter;gbkey=CDS;gene=GCN20;locus_tag=CAALFM_C100480CA;orig_transcript_id=gnl|WGS:AACQ|mrna_CAALFM_C100480CA;product=putative AAA family ATPase;protein_id=XP_719022.1;transl_table=12
> NC_032089.1	RefSeq	exon	75573	77828	.	-	.	ID=exon-XM_713929.2-1;Parent=rna-XM_713929.2;Dbxref=GeneID:3639314,Genbank:XM_713929.2;end_range=77828,.;gbkey=mRNA;gene=GCN20;locus_tag=CAALFM_C100480CA;orig_protein_id=gnl|WGS:AACQ|CAALFM_C100480CA;orig_transcript_id=gnl|WGS:AACQ|mrna_CAALFM_C100480CA;partial=true;product=putative AAA family ATPase;start_range=.,75573;transcript_id=XM_713929.2
> NC_032089.1	RefSeq	gene	75573	77828	.	-	.	ID=gene-CAALFM_C100480CA;Dbxref=GeneID:3639314;Name=GCN20;end_range=77828,.;gbkey=Gene;gene=GCN20;gene_biotype=protein_coding;locus_tag=CAALFM_C100480CA;partial=true;start_range=.,75573
> NC_032089.1	RefSeq	mRNA	75573	77828	.	-	.	ID=rna-XM_713929.2;Parent=gene-CAALFM_C100480CA;Dbxref=GeneID:3639314,Genbank:XM_713929.2;Name=XM_713929.2;end_range=77828,.;gbkey=mRNA;gene=GCN20;locus_tag=CAALFM_C100480CA;orig_protein_id=gnl|WGS:AACQ|CAALFM_C100480CA;orig_transcript_id=gnl|WGS:AACQ|mrna_CAALFM_C100480CA;partial=true;product=putative AAA family ATPase;start_range=.,75573;transcript_id=XM_713929.2
> ```

We see there are 4 annotations associated with `GCN20`, each a different Feature type (`CDS`, `exon`, `gene`, `mRNA`), 
and they all point toward the same 2256 nt on chromosome 1.

Let's choose the Feature type `gene`, and its corresponding attribute ID `gene-CAALFM_C100480CA`.

We will use a Target motif, an on-target score, and an off-target score each appropriate for Cas9. We use default score weights for both `Azimuth` and `CFD`. We want to narrow the specificity by broadening the number of sequences that can be considered off-target, so we specify the `--off_target_motifs` option.

We will keep the rest of the AddTag default options. Our final command to identify the best Target sequences and generate the dDNA is the following:
```sh
addtag generate_all \
  --features gene \
  --selection gene-CAALFM_C100480CA \
  --motifs 'N{17}|N{3}>NGG' \
  --off_target_motifs 'N{17}|N{3}>NAG' \
  --ontargetfilters Azimuth \
  --offtargetfilters CFD \
  --excise_insert_lengths 0 0 \
  --ko-gRNA \
  --ko-dDNA mintag \
  --fasta ${GENOME}.fna \
  --gff ${GENOME}.gff \
  --folder ${GENE}g > ${GENE}g.out 2> ${GENE}g.err
```

This will output a single table, with the best Targets in the top of the output, and the worst toward the bottom.
```sh
head ${GENE}g.out 
```
> ```
> # exTarget results
> # gene	features	weight	exTarget name	exTarget sequence	OT:CFD	Azimuth	reDonors	None	feature:contig:strand:start..end	warnings
> gene-CAALFM_C100480CA	gene-CAALFM_C100480CA	0.876615628563555	exTarget-96	CCAACGAAACAGTTTTCAGG>GGG	71.63	64.65		None	gene-CAALFM_C100480CA:NC_032089.1:-:76719..76742	None
> gene-CAALFM_C100480CA	gene-CAALFM_C100480CA	0.6877504215295679	exTarget-110	CATTATTACGTGCCTTGTCG>AGG	57.83	57.84		None	gene-CAALFM_C100480CA:NC_032089.1:-:77090..77113	None
> gene-CAALFM_C100480CA	gene-CAALFM_C100480CA	0.650732417867368	exTarget-86	CTCTTTCTATGCAACTCGTG>AGG	49.97	59.21		None	gene-CAALFM_C100480CA:NC_032089.1:-:76456..76479	None
> gene-CAALFM_C100480CA	gene-CAALFM_C100480CA	0.6497399064714101	exTarget-84	ACAGTCTCGTATCAAGAAGT>TGG	47.83	61.06		None	gene-CAALFM_C100480CA:NC_032089.1:-:76324..76347	None
> gene-CAALFM_C100480CA	gene-CAALFM_C100480CA	0.6373121738894936	exTarget-21	GACTTTCGTATTCACGACGT>TGG	61.75	55.93		None	gene-CAALFM_C100480CA:NC_032089.1:+:76420..76443	None
> gene-CAALFM_C100480CA	gene-CAALFM_C100480CA	0.6105741443075372	exTarget-117	GAGCGAGGCGTCATTGACAT>TGG	61.2	55.21		None	gene-CAALFM_C100480CA:NC_032089.1:-:77164..77187	None
> gene-CAALFM_C100480CA	gene-CAALFM_C100480CA	0.49714752260403006	exTarget-56	GGATGAACCGTCCAATCACT>TGG	49.7	54.11		None	gene-CAALFM_C100480CA:NC_032089.1:-:75787..75810	None
> gene-CAALFM_C100480CA	gene-CAALFM_C100480CA	0.4582728976377312	exTarget-89	AGATATAATCCATCAACACT>CGG	40.04	66.96		None	gene-CAALFM_C100480CA:NC_032089.1:-:76513..76536	None
> ```

If you run this command again, but omit the `--off_target_motifs` option, you get the following:
> ```
> # exTarget results
> # gene	features	weight	exTarget name	exTarget sequence	OT:CFD	Azimuth	reDonors	None	feature:contig:strand:start..end	warnings
> gene-CAALFM_C100480CA	gene-CAALFM_C100480CA	0.8776689032281014	exTarget-96	CCAACGAAACAGTTTTCAGG>GGG	74.29	64.65		None	gene-CAALFM_C100480CA:NC_032089.1:-:76719..76742	None
> gene-CAALFM_C100480CA	gene-CAALFM_C100480CA	0.8360644970889776	exTarget-95	CAACGAAACAGTTTTCAGGG>GGG	50.0	74.38		None	gene-CAALFM_C100480CA:NC_032089.1:-:76718..76741	None
> gene-CAALFM_C100480CA	gene-CAALFM_C100480CA	0.8056514544867718	exTarget-84	ACAGTCTCGTATCAAGAAGT>TGG	91.67	61.06		None	gene-CAALFM_C100480CA:NC_032089.1:-:76324..76347	None
> gene-CAALFM_C100480CA	gene-CAALFM_C100480CA	0.7897087206837166	exTarget-30	GTTTAACTCTCTCCTCGACA>AGG	49.97	67.38		None	gene-CAALFM_C100480CA:NC_032089.1:+:77078..77101	None
> gene-CAALFM_C100480CA	gene-CAALFM_C100480CA	0.7553979473976637	exTarget-86	CTCTTTCTATGCAACTCGTG>AGG	76.79	59.21		None	gene-CAALFM_C100480CA:NC_032089.1:-:76456..76479	None
> gene-CAALFM_C100480CA	gene-CAALFM_C100480CA	0.7126907697600093	exTarget-110	CATTATTACGTGCCTTGTCG>AGG	73.09	57.84		None	gene-CAALFM_C100480CA:NC_032089.1:-:77090..77113	None
> gene-CAALFM_C100480CA	gene-CAALFM_C100480CA	0.6433552008729548	exTarget-100	GTATGGTTTGGGTTTCACAA>AGG	72.38	55.81		None	gene-CAALFM_C100480CA:NC_032089.1:-:76756..76779	None
> gene-CAALFM_C100480CA	gene-CAALFM_C100480CA	0.6373121738894936	exTarget-21	GACTTTCGTATTCACGACGT>TGG	61.75	55.93		None	gene-CAALFM_C100480CA:NC_032089.1:+:76420..76443	None
> ```

Notice that by including the additional off-target motif, we see generally lower off-target scores (the `OT:CFD` column).

Next we will identify the best cPCR primers for verifying the 'GCN20' full CDS deletion.


</td></tr></tbody></table>
</details>

#### 2-step deletion of a single Feature ####
<details>
<summary>Click to expand/collapse</summary>
<table><tbody><tr><td>

We will delete a Feature that has no Target within it.

*~ Section incomplete ~*

</td></tr></tbody></table>
</details>

#### 1-step editing of a single Feature ####
<details>
<summary>Click to expand/collapse</summary>
<table><tbody><tr><td>

We will edit a Feature

*~ Section incomplete ~*

</td></tr></tbody></table>
</details>

#### 2-step editing of a single Feature ####
<details>
<summary>Click to expand/collapse</summary>
<table><tbody><tr><td>

We will edit a Feature that has no Target within it.

*~ Section incomplete ~*

</td></tr></tbody></table>
</details>

#### 2-step deletion and add-back of a single Feature ####
<details>
<summary>Click to expand/collapse</summary>
<table><tbody><tr><td>

In this example, we will go through creating a nominal `mintag` to knock-out a single input Feature, then creating primers necessary to revert back to the wild type Feature.

The standard procedure is to first run `addtag generate_all`, and use its output as input for `addtag generate_primers`.

For simplicity, We will assume the name of the Feature you are interested in is `GENE`.

The first thing you will want to do, is compose a Target motif for the [![RGN][RGN]](#) your biological system uses. To see a list of commonly-used Target motifs, run the following:
```sh
addtag list_motifs
```
Let's pretend our biological system uses the 'AsCpf1' [![RGN][RGN]](#). So we will use the associated `TTTN<N{19}/.{4,6}\` Target motif. Thus, we will add `--motifs 'TTTN<N{19}/.{4,6}\'` to the `addtag generate_all` command.

The next step is to select one or more Algorithms to calculate the 'on-target' and 'off-target' scores for this [![RGN][RGN]](#). To see a list of all implemented Algorithms, run the following:
```sh
addtag list_algorithms
``` 
Let's choose the `DeepCpf1` Algorithm for our 'on-target' score. Let's also choose the `Linear` Algorithm for the 'off-target' score, whose implicit behavior severely penalizes insertions and deletions at 'off-target' sites, but is explicitly less biased against mismatches. Therefore we add `--ontargetfilters DeepCpf1 --offtargetfilters Linear` to the command. Because we would like the output Target sites to be ranked based on their specificity, and because the `Linear` algorithm does not have a default weight, we define a weight for it using the `--weights` command line option.

Let's use the `mintag` method for creating an RGN Target on the dDNA we generate for creating the intermediary genome. Because we don't want to add any extra bases--only remove the feature--we include `--excise_insert_lengths 0 0`. Finally, we want merely to revert back to wild type at the native locus, so we direct AddTag to generate the optimal AmpF/AmpR primers using the `--revert_amplification_primers` option.

Because our input genome is a phased diploid assembly, and we want our gRNAs to target both alleles, we use the default `--target_specificity`. Because we want a single dDNA to repair both alleles, we also use the default `--donor_specificity`. Since we want the computer to use all available compute power, we use the default number of processors (which automatically selects all available). Let's also use the default thermodynamics calculator and the default aligner.

Let's store all the output in paths that start with `GENEg`, where '`g`' is for '`generate_all`'.
 
To identify the best Target locations within our Feature of interest, and to generate dDNA for knock-out, we run the full command:
```sh
addtag generate_all \
  --motifs 'TTTN<N{19}/.{4,6}\' \
  --ontargetfilters DeepCpf1 \
  --offtargetfilters Linear \
  --weights Linear:85+1.7 \
  --excise_insert_lengths 0 0 \
  --ko-gRNA \
  --ko-dDNA mintag \
  --revert_amplification_primers \
  --fasta genome.fasta \
  --gff genome.gff \
  --folder GENEg > GENEg.out 2> GENEg.err
```

This writes 4 output tables to the `GENEg.out` file. Each of these tables refers to sequences in output FASTA files. Please note that certain sequence Aligners, such as 'Bowtie2' can have non-deterministic output. Therefore, your results may vary from what is presented here.

Now would be a good time to explain the terminology you will see in the AddTag input and output. For simplicity in text processing, we use different labels than what are presented in the manuscript, though they are equivalent.
```
OUTPUT           PAPER    DESCRIPTION
r0-gDNA          +gDNA    Wild type genome
r1-gDNA          ΔgDNA    Intermediary genome
r2-gDNA          AgDNA    Final genome
exTarget         +Target  Target site in wild type +Feature that is used to 'excise' the feature
reTarget         ΔTarget  Target site introduced with ★tag insert that is used to 'revert' the genotype
exDonor/r1-dDNA  ΔdDNA    Excision, or knock out dDNA (ko-dDNA)
reDonor/r2-dDNA  AdDNA    Reversion, add-back, or knock-in dDNA (ki-dDNA)  
```

Thus, we refer to the first round of genome engineering (r1) as the knock-out round, and the second round (r2) as the knock-in round.

From the first table, we select the highest-weighed `reTarget` ('reversion Target', abbreviated), and one of its corresponding dDNA sequences. Then we store these two sequences in their own FASTA files.

```sh
addtag find_header \
  --fasta GENEg/reversion-spacers.fasta 
  --query 'reTarget-3\b' > ki-target.fasta
```

Each `reTarget` can target one or more identified `reDonor` dDNA sequences.

In this example, we expect only a single excision dDNA with the header `exDonor-0`, so we extract that sequence, and store it in a conveniently-accessible place.

```sh
addtag find_header \
  --fasta GENEg/excision-dDNAs.fasta \
  --query 'exDonor-0\b' > ko-dDNA.fasta
```

From the second table, we select the highest-weighted `exTarget` ('excision Target' abbreviated), which is used for excising the input Feature from the input gDNA:

```sh
TARGET=$(grep '# exTarget results' -A 2 GENEg.out | tail -n +3 | cut -f 4)
addtag find_header \
  --fasta GENEg/excision-spacers.fasta 
  --query "${TARGET}\b" > ko-target.fasta
```

Finally, we identify the highest-weight dDNA for reverting back to the wild type, and put it in its own FASTA file:
```sh
addtag find_header --fasta GENEg/reversion-dDNAs.fasta --query 'reDonor-0\b' > ki-dDNA.fasta
```

Next we need to identify a single cPCR verification primer design. Let's use the default pairwise sequence aligner.

```sh
addtag generate_primers \
  --fasta genome.fasta \
  --dDNAs ko-dDNA.fasta ki-dDNA.fasta \
  --folder GENEc > GENEc.out 2> GENEc.err
``` 

</td></tr></tbody></table>
</details>

#### 2-step deletion and add-back of a single, phased Feature ####
<details>
<summary>Click to expand/collapse</summary>
<table><tbody><tr><td>

*~ Section incomplete ~*

</td></tr></tbody></table>
</details>

#### 2-step editing of several Features ####
<details>
<summary>Click to expand/collapse</summary>
<table><tbody><tr><td>

*~ Section incomplete ~*

</td></tr></tbody></table>
</details>

#### Multiplexed, 2-step editing of several Features ####
<details>
<summary>Click to expand/collapse</summary>
<table><tbody><tr><td>

All Features in input GFF file will be evaluated simultaneously.

*~ Section incomplete ~*

</td></tr></tbody></table>
</details>

## 📝 Citing AddTag ##
If you use AddTag for your research, please cite us. Because the manuscript is currently pending review, you will need to cite the code repository instead.

 > Thaddeus D. Seher, Diana Ramos, Namkha Nguyen, Priyanka Bapat, Clarissa J. Nobile, Suzanne S. Sindi, and Aaron D. Hernday. AddTag: software for automated design and validation of precision CRISPR/Cas genome edits. University of California, Merced. Retrieved from <[https://github.com/tdseher/addtag-project](https://github.com/tdseher/addtag-project)> (2019).

## ✍ Authors ##
Who do I talk to?
 * Aaron D. Hernday (🔬 PI leading the project)
 * Thaddeus D. Seher (💻 programmer) (💬[@tdseher][tdseher])

See also the list of [contributors](https://github.com/tdseher/addtag-project/graphs/contributors) who participated in this project.

## 👥 Contributing ##

### 🤔 What can I do to help improve AddTag? ###
<details>
<summary>Click to expand/collapse</summary>
<table><tbody><tr><td>

We are always looking for ways to broaden the usability of the AddTag software. Here is a list of things that would be great contributions.
 * Improvements to the documentation, such as additional example workflows.
 * More Target motifs (SPACER≷PAM combinations) from new CRISPR/Cas literature to add to the `motifs.txt` file.
 * Support for additional pairwise sequence Aligners.
 * Support for additional scoring Algorithms.
 * Support for additional thermodynamics calculators.
 * Running AddTag on different types of genomes with different parameters to test proper logic and assess compatibilities.

</td></tr></tbody></table>
</details>

### 🐞 How do I submit a bug report? ###
<details>
<summary>Click to expand/collapse</summary>
<table><tbody><tr><td>

First, check to see if the problem you are having has already been added to the [issue tracker](https://github.com/tdseher/addtag-project/issues).
If not, then please submit a new issue.

</td></tr></tbody></table>
</details>

### ⚠ How do I make a feature request? ###
<details>
<summary>Click to expand/collapse</summary>
<table><tbody><tr><td>

Send a message to [@tdseher][tdseher].

</td></tr></tbody></table>
</details>

### ⤴ How do I add my code to the AddTag software? ###
<details>
<summary>Click to expand/collapse</summary>
<table><tbody><tr><td>

Please submit a [pull request](https://github.com/tdseher/addtag-project/pulls).

</td></tr></tbody></table>
</details>

### 📈 Adding scoring Algorithms ###
<details>
<summary>Click to expand/collapse</summary>
<table><tbody><tr><td>

Scoring Algorithms have been broken down into two general types.

 * `SingleSequenceAlgorithm` objects calculate scores by comparing a potential RNA [![Spacer][Spacer]](#) or DNA [![Target][Target]](#) to a model trained on empirical data.
 * `PairedSequenceAlgorithm` instances generate scores that compare a potential RNA [![Spacer][Spacer]](#) to a DNA [![Target][Target]](#).

To add a new scoring algorithm, you must subclass one of the the above types, and add it to a `*.py` file in the `source/algorithms/` subdirectory. AddTag will automatically calculate the score on every generated [![Spacer][Spacer]](#).

We welcome any `git pull` requests to widen the repertoire of scoring algorithms available to AddTag. The easiest way to get started is to copy and modify one of the provided subclasses.

</td></tr></tbody></table>
</details>

### 📐 Adding sequence Aligners ###
<details>
<summary>Click to expand/collapse</summary>
<table><tbody><tr><td>

AddTag comes with wrappers for several alignment programs. Depending on your experimental design and computing system, you may decide to use an aligner with no included wrapper. To implement your own, create a subclass of `Aligner`, and put it in a `*.py` file in the `source/aligners/` subdirectory. AddTag will automatically make that aligner available for you.

Share your code with us so we can make it available to all AddTag users.

</td></tr></tbody></table>
</details>

### 🌡 Adding Thermodynamics calculators ###
<details>
<summary>Click to expand/collapse</summary>
<table><tbody><tr><td>

Several wrappers to popular oligonucleotide conformation, free energy, and melting temperature calculation programs are included. You can add your own by subclassing the `Oligo` class, and then adding its `*.py` file to the `source/thermodynamics/` subdirectory.

If you create your own wrapper, please submit a `git pull` request so we can add it to the next version of the software.

</td></tr></tbody></table>
</details>

## 📖 License ##
Please see the [LICENSE.md](LICENSE.md) file.

## Notes ##
Below are tips and descriptions of AddTag limitations that will help you make successful designs.
<details>
<summary>Click to expand/collapse</summary>
<table><tbody><tr><td>

 * If you are identifying cPCR primers, then it is often useful to use the `--cache` option. This lets you decrease the stringency of the PCR conditions and run the `generate_primers` subroutine again, pointing to the same `--folder`, and AddTag will use the results from the previous calculations when it can instead of doing the computations from scratch.
 * The [![RGN][RGN]](#) protein you use should be engineered specifically for your organism. It should be codon-optomized, and if using eukarya, contain an appropriate nuclear localization sequence.
 * By default, AddTag will avoid designing homology regions and Targets against polymorphisms whenever possible.
 * Sequences in FASTA files should have unique names. In other words, the primary sequence identifier--everything following the '`>`' character and preceding the first whitespace/tab '` `' character--should exist only once across all input `*.fasta` files.
 * AddTag makes no effort to restrict which Target motifs the user can use according to the selected Algorithms. Therefore, the user needs to independently verify which Target motifs are compatible with the selected Algorithms.
 * Right now AddTag can only handle linear chromosomes. If you want to analyze a circular chromosome, then you will need to artificially concatenate the ends of the chromosome together and adjust any annotations before running AddTag. An additional complication the software does not address is circular chromosomes. Features and their flanking regions cannot span the junction created when the contig end is concatenated to the start (typically the starting position on a contig is labeled the ORIGIN). To address this, the user should manually shift the coordinates of the experimental Features, and wrap the contigs as appropriate.
 * AddTag assumes one Feature copy per contig. The current implementation of AddTag assumes homology regions around Features are not repeated across any one contig. This means that is will fail to generate cPCR oligos for a large proportion of genes in transposon-rich genomes such as [wheat](https://dx.doi.org/10.1186/s13059-018-1479-0). This limitation is currently a result of both the *in silico* recombination and the primer identification routine. If there are tandem Features on a contig, then the sF and sR primers are likely duplicated across these adjacent loci. The shared primers thus can't specifically amplify one of the tandem duplications and not the other.
 * AddTag uses the *in silico* recombination phase of `generate_primers` subroutine to determine if flanking homology regions of dDNAs are too repetetive across the genome (ideally, this would be performed in the `generate_all` subroutine).
 * A single Feature cannot span two or more contigs (partially a limitation of the GFF format). AddTag assumes that the entire feature sequence, and any flanking regions, are not in terminal regions of the reference contig. 
 * AddTag does not address overlapping genes, such as when an intron contains an exon for another gene, or when the same DNA encodes for genes on opposite strands. Everything between the Feature bounds is removed in the first engineering step. Currently, if the selected Feature overlaps with any other feature, only the selected Feature is considered. The other Feature will be disrupted. AddTag will report a warning that these other Features may be disrupted, but it does not attempt to reconcile this in any way. However, AddTag does have the ability to limit Feature expansion to keep the deletion outside of neighboring Features.
 * AddTag was not designed to perform paired Cas design, such as FokI-dCas9 nickase. You would need to run the program and select two gRNAs designed for opposite strands within a certain distance from each other. Alternatively, you could probably make some really-long Target motif. One way to mitigate errors is to use PAM-out nickases. This requires Cas9 cutting by two targets to get double-stranded break. This significantly decreases off-target genome editing. However, this initial AddTag version does not explicitly facilitate this.
 * AddTag can identify cut sites for Cas enzymes which have the PAM site. No functionality is provided for finding sites without an adjacent PAM sequence. AddTag requires motifs to define a PAM sequence. Therefore Cas14a is not supported. This can be probably be circumvented by using an `N` character as the PAM sequence, but this hasn't been tested. The number of CRISPR/Cas genome editing technologies are rapidly growing. With the recent discovery of [Cas14a](https://dx.doi.org/10.1126/science.aav4294), which targets single-stranded DNA (ssDNA) molecules without requiring a PAM site, the expanded prevalence of CRISPR/Cas methods in biological sciences is assured. However, often researchers wish to edit sites on double-stranded DNA (dsDNA) using an RGN (such as Cas9 or Cas12a) that requires binding to a PAM motif. 
 * Please note, that at this time, no special restriction sites will be taken into account when designing primers.
 * For simplicity, all calculated scores ignore terms dealing with proximity to exon/CDS/ORF sequences. In cases such as the Stemmer and Azimuth calculations, the authors attempted to include the risk of disrupting genes neighboring potential targets in their models. We don’t attempt to do this.
 * Additionally, some scoring Algorithms take chromatin structure (DNA accessibility) into account. For simplicity, AddTag treats all input gDNA as equally accessible.
 * During the course of writing this software, a [paper](https://dx.doi.org/10.1038/s41587-019-0095-1) was published that outlines how hairpins can be inserted into the pre-spacer and spacer regions of the gRNA in order to increase specificity. AddTag does not model pre-spacer sequences.
 * AddTag assumes the RGN template type is dsDNA. AddTag was designed specifically to enable efficient gDNA editing. It does not use predictive models for ssDNA or RNA templates.
 * A corollary of this is that AddTag assumes all input sequences are DNA sequences. So the `--fasta` file specified will be treated as a DNA template. Thus, if there are any non-DNA residues, such as `U`, AddTag will probably fail. Also, since the Primer thermodynamics calculators are all set to estimate DNA:DNA hybridization (not DNA:RNA or RNA:RNA), any resulting calculations will be incorrect. 
 * Since Bartag motifs are user-specified, simple pre-computed lists of compatible 'bartag' sequences would be incomplete. Thus we implemented a greedy 'bartag' generation algorithm. When evaluating candidate 'bartag' sequences, AddTag will keep 'bartags' that satisfy all edit distance requirements with all previously-accepted 'bartags'. To limit runtime to a reasonable amount, we limited the total number of Features and 'bartags' that can be generated.
 * Of special note are things the Primer design does not explicitly consider, such as characteristics of the cPCR template molecule. AddTag does not exploit the differential nature of template sequence composition (e.g. H. sapiens compared to E. coli). Also, AddTag does not use information on the presence of known secondary modifications to the template, such as methylated residues or oxidative damage.
 * One of the big limitations of this version of AddTag is that the Primer attribute stringencies are held uniform across all regions. You specify this using the `--cycle_start N` and `--cycle_stop N` options. If any one of the desired Primer Pairs is not found under the selected stringency, then no simulated annealing is performed. Cycles range from `N` of 0 to 21, with 0 being the most restrictive, and 21 being the most permissive. Due to the brute-force nature of the Primer Pair calculations, increasing `N` will exponentially increase the amount of memory needed to evaluate primers. So if you increase the cycles, be sure to monitor system RAM.
 * To facilitate more straightforward programming, AddTag outputs 0-based genomic coordinates (as opposed to traditional 1-based coordinates). All input data, such as `GFF` files, are expected to use 1-based genomic coordinates.
 * If Algorithm columns in the `STDOUT` of the `generate_all` subroutine, return `0.0`, then a likely cause is that the Algorithm prerequisites are not correctly installed. For instance, if `Azimuth` scores are all `0.0` on a Linux machine, then you might be missing the `python-tk` system package. In this case, try to run `source/algorithms/addtag_wrapper.py` in isolation to troubleshoot the problem.
 * The 'forward' and 'reverse' cognomina are absolute to the input contig coordinates. For instance, the 'sF' and 'sR' Primers are not relative to the orientation of the Feature defined in the input GFF. Instead, 'sF' is earlier in the contig (lower number), and 'sR' is later in the contig (higher number).
 * In this current version, AddTag's `generate_primers` subroutine assumes that there is a double-stranded break (DSB) in the gDNA between the the locations the dDNA homology arms are similar to. Furthermore, it assumes these DSBs are repaired perfectly though homology-directed repair (HDR). This makes sense in our experimental biological system *C. albicans*. In other systems, such as *H. sapiens*, there is a higher amount of error-prone DSB repair.
 * Run time is a function of the number of potential primers that need to be analyzed. Thus, genes that are longer have more potential primers. Also, the number of potential primers actually analyzed depends on the sequence composition of each region. If a region has great complexity, then more primers will be analyzed with the full suite of filters, and the analysis will take longer. If a region has little complexity, then more potential primers will be discarded at early filters, and the analysis will take less time.
 * In rare cases, endogenous RNA may bind to the RGN to drive cutting at non-target loci. AddTag does not screen the input gDNA for this possibility because it does not analyze the scaffold section of the gRNA. 

</td></tr></tbody></table>
</details>

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
