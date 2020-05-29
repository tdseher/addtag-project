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

[Features](#-features) • [Requirements](#-requirements) • [Installing](#-installing-addtag) • [Usage](#-program-instructions) • [Aligners](#-supported-sequence-aligners) • [Thermodynamics](#-supported-thermodynamics-calculators) • [Algorithms](#-supported-scoring-algorithms) • [Citing](#-citing-addtag) • [Contributing](#-contributing)

## ☑ Features ##
Basic Features:
 * [x] Analyzes any arbitrary genomic DNA (gDNA).
   * [x] Fully supports ambiguous characters or polymorphisms (`RYMKWSBDHVN`).
   * [x] Respects case-masked gDNA for ![Target][Target] and ![Primer][Primer] identification.
 * [x] Locates RNA-guided nuclease (![RGN][RGN]) cut sites (![Target][Target]s) within a ![Feature][Feature] (locus of interest) for optimal gRNA ![Spacer][Spacer]s.
   * [x] Fully supports ambiguous bases (`RYMKWSBDHVN`) in ![Spacer][Spacer] or ![PAM][PAM].
   * [x] Accepts 3'-adjacent ![PAM][PAM] sequences, such as Cas9 (`>NGG`).
   * [x] Accepts 5'-adjacent ![PAM][PAM] sequences, such as Cas12a (`TTTN<`).
   * [x] Supports arbitrary ![Spacer][Spacer] length and composition constraints, such as for plant experiments (`G{,2}N{19,20}`).
   * [x] Supports arbitrary ![PAM][PAM] sequences (MAD7: `YTTN<`, Cas12d: `TA<`, BlCas9: `>NGGNCNDD`, etc).
   * [x] Uses stranded forward (`/`), reverse (`\`) and unstranded (`|`) cut sites.
   * [x] Supports complex nested logic ![PAM][PAM] sequences, such as xCas9 (`>(N{1,2}G,GAW,CAA)`)
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

For tips on setting up AddTag requirements, please review the commands in the `.azure-pipelines.yml` file.

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

 * [x] GC ([Wang, et al (2014)](http://dx.doi.org/10.1126/science.1246981'))

 * [x] Homopolymer ([Hough, et al. (2017)](https://doi.org/10.1186/s12859-017-1581-4))

 * [x] PolyT

 * [x] PAM Identity

 * [x] Position

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

To make the AddTag executable accessible from any working directory, you can add the absolute path of the current working directory to the `PATH` variable.

On Windows, run:
```sh
set PATH=%PATH%;%CD%
```

On Linux or macOS, run:
```sh
export PATH=$PATH:$PWD
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

## 💻 Program Instructions ##
### Displaying the usage ###
Because AddTag is being updated regularly, the most current feature set and usage can be viewed by running AddTag with the `--help` command line option.

The following commands assume the current working directory is the AddTag folder `addtag-project/`. This will print out command line parameter descriptions and examples.
```sh
./addtag --help
```

Additionally, you may view the included man page, which is probably not up-to-date.
```sh
man ./addtag.1
```

### Formatting input data ###

#### FASTA input ####
AddTag requires a FASTA genome of organism you wish to manipulate. FASTA files resemble the following:
```fasta
>primary_header1 attribute1=value1 attribute2=value2
NNNNCGAAATCGGCGCATAGGCCTAAGAGCTCCTATAGAGATCGATATAAAC
GCTAGAAGATAGAGAGAGGCTCGCGCTCGATCGCGATAAGAGAGCTCTCGGC
CGATAGAGGAATCTCGggctcgcatatatyhcgcggcatatGGCCTAGAGGA
CCAATAAAGATATATAGCCTAAAGGAATATATAGAGAGATATATATAGNNNN
>primary_header2 attribute1=value1 attribute2=value2
AGCTAGAGACWWWCTCCTCTCCTAGAGASSSAGAGGAGAGCTCTCCGAGAGA
CGCTCGCTCGTATGCCTCTATATCGATATATAGGAGAATCCTCGATATATAG
```
FASTA files are plain text files that use newline (`\n` or `\r\n`) characters as delimiters. If a line begins with a greater than (`>`) symbol, it represents the start of a new sequence record. All characters between the `>` and `\n` are considered the 'header' of the record. Everything between the `>` and the first whitespace character (` ` or `\t`), if one exists, is considered the 'primary header' for the record. All subsequent lines until the next 'header' line contain the sequence information for that record. Therefore FASTA files can contain many sequence records. Each record in the FASTA file in an assembly is called a 'contig'. 

FASTA files can contain any number of ambiguous characters (`RYMKWSBDHVN`), which can represent allelic variation expected within the sample or sequencing uncertainty. FASTA files can also contain a mix of `UPPER` and `lower` cased characters. Typical use for `lower` case characters is to exclude these residues from ![Target][Target] or ![Primer][Primer] identification. 

#### GFF input ####
AddTag requires a GFF file containing annotations for the Features you wish to manipulate. GFF files resemble the following:
```gff
# seqid	source	feature	start	end	score	strand	frame	attribute
C1A	DB	gene	3489	5146	.	+	.	ID=C1A_001;Name=C1_001;Gene=GENE1
C1A	DB	mRNA	3489	5146	.	+	.	ID=C1A_001-T;Parent=C1A_001
C1A	DB	exon	3489	5146	.	+	.	ID=C1A_001-T-E1;Parent=C1A_001-T
C1A	DB	CDS	3489	5146	.	+	0	ID=C1A_001-P;Parent=C1A_001-T

C1B	DB	gene	3267	4924	.	+	.	ID=C1B_001
C1B	DB	mRNA	3267	4924	.	+	.	ID=C1B_001-T;Parent=C1B_001
C1B	DB	exon	3267	4924	.	+	.	ID=C1B_001-T-E1;Parent=C1B_001-T
C1B	DB	CDS	3267	4924	.	+	0	ID=C1B_001-P;Parent=C1B_001-T
```
GFF files describe the contig locations of important genomic Features. Empty lines and lines that begin with the pound (`#`) symbol are ignored. Of note is the far-right `attribute` column, which AddTag assumes is a semicolon-delimited set of key/value pairs. AddTag assumes each Feature has a unique identifier. By default, it uses the `ID` attribute as the unique name for each Feature. If your GFF file does not have an `ID` attribute, then you can select a different one with the `--tag` command line option. 

Typical AddTag analyses require at least one GFF file. AddTag can handle GFF files in two ways.
 * For the first method, all Features matching the selected type, designated by the `--features` command line argument, will be included for analysis. By default, only lines in the GFF file containing `gene` in `feature` column will be considered. This system is useful if your GFF file contains only the Features you wish to manipulate.
 * If your GFF file contains all annotations for the entire genome (which is typical), the second approach requires you to select only the few Features you want to edit using the `--selection` command line argument.
 
Often, you will have a GFF file with annotations for the entire genome. The `attributes` column is not often structured intuitively, and can prove cumbersome to search (`grep`) or sort (`sort`) manually. To make it easy to identify the desired lines of a GFF file, AddTag includes the `find_feature` subroutine. Here is an example that tries to find all lines associated with `HSP90` by searching several attribute tags, and outputting a GFF with a commented line containing field names:
```sh
addtag find_feature --gff genome.fasta --query HSP90 --linked_tags Name Alias Parent Gene --header > features.gff
```

#### Target motif input  ####
The Target motif is written from 5' to 3'. Use a greater than (`>`) symbol if your RGN has a 3'-adjacent PAM, and use a less than (`<`) symbol if your RGN has a 5'-adjacent PAM. Ambiguous nucleotide characters are accepted. `{a,b}` are quantifiers. `(a,b,…)` are permitted alternatives. `/` is a sense strand cut, `\` is an antisense strand cut, and `|` is a double-strand cut. `.` is a base used for positional information, but not enzymatic recognition. Be sure to enclose each motif in quotes so your shell does not interpret `STDIN`/`STDOUT` redirection.

You can specify any number of Target motifs to be considered 'on-target' using the `--motifs` command line option. You can also designate any number of Target motifs to be considered 'off-target' using the `--off_target_motigs` command line option. 

#### Homologs input ####
Some researchers are lucky enough to get to work on organisms with phased genomes. This means that full haplotype information is known for each chromosome. AddTag can accommodate haploid, diploid, and polyploid genomes when homologous Features are linked by the addition of the `--homologs` command line option. The 'homologs' file has the following format:
```homologs
# group	homolog_a	homolog_b	homolog_c
GENE1	C1A_001 	C1B_001
GENE2	C1A_002 	C1B_002 	C1C_002
```
Each Feature identifier has its contig start and end position defined in the input GFF file. The 'homologs' file merely links them together. Columns in the homologs file are delimited by the `\t` character. The first column is the name of the group of Features. Every subsequent column should contain the identifier of a Feature to consider as a homolog. Homolog groups can each have any number of Features. If a Feature identifier appears on multiple lines, then all those Features are linked together as one homolog group. The identifier can be changed with the `--tag` command line option.  

### Available subroutines ###
The AddTag program contains a set of subroutines that can be run independently. There are four categories of subroutines.

 * The `evaluate_*` subroutines run only a very specific analysis on input data.
 * The `find_*` subroutines are used to search input files for specific things, so the user can easily learn the correct parameters to use for AddTag input.
 * The `generate_*` subroutines perform the deep computational analyses.
 * The `list_*` subroutines just print information the user might find useful.

### Available RGN scoring Algorithms ###
Over the past few years, several Algorithms have been proposed to describe RGN behavior within certain biological contexts. We implemented most of the commonly-used ones into the AddTag software. To view information about each, use the following command:
```sh
addtag list_algorithms
``` 
This will write the pertinent information for all implemented Algorithms to `STDOUT`.

If an Algorithm is used for pre-alignment filtering (`Prefilter`) or post-alignment filtering (`Postfilter`), then the score of the Target must lie between the `Min` and `Max` values to be continued on through the analysis. For instance, the 'off-target' scoring `CFD` Algorithm has a `Min` of `1.0`. This means that some positions with significant sequence similarity to the query Target (because they are identified in the Alignment step) will not contribute to the final 'off-target' score if their score is less than `1.0`.

### Available oligonucleotide thermodynamics calculators ###
To view which thermodynamics calculators are available on your system, use the following command:
```sh
addtag list_thermodynamics
```

### Typical workflow for a single Feature ###
The standard procedure is to first run `addtag generate_all`, and use its output as input for `addtag generate_primers`.

The first thing you will want to do, is compose a Target motif for the RGN your biological system uses. To see a list of commonly-used Target motifs, run the following:
```sh
addtag list_motifs
```
For simplicity, let's pretend our biological system uses the 'AsCpf1' RGN. So we will use the associated `TTTN<N{19}/.{4,6}\ ` Target motif.

The next step is to select one or more Algorithms to calculate the 'on-target' and 'off-target' scores for this RGN. To see a list of all implemented Algorithms, run the following:
```sh
addtag list_algorithms
``` 
Let's choose the `DeepCpf1` Algorithm for our 'on-target' score. Let's also choose the `Linear` Algorithm for the 'off-target' score, whose implicit behavior severely penalizes insertions and deletions at 'off-target' sites, but is explicitly less biased against mismatches. Because we would like the output Target sites to be ranked based on their specificity, and because the `Linear` algorithm does not have a default weight, we define a weight for it using the `--weights` command line option.

Let's use the `mintag` method for creating an RGN Target on the dDNA we generate for creating the intermediary genome. Finally, we want merely to revert back to wild type at the native locus, so we direct AddTag to generate the optimal AmpF/AmpR primers using the `--revert_amplification_primers` option.

Because our input genome is a phased diploid assembly, and we want our gRNAs to target both alleles, we use the default `--target_specificity`. Because we want a single dDNA to repair both alleles, we also use the default `--donor_specificity`. Since we want the computer to use all available compute power, we use the default number of processors (which automatically selects all available). Let's also use the default thermodynamics calculator and the default aligner. To identify the best Target locations within our Feature of interest, and to generate dDNA for knock-out, we run the following command:
```sh
addtag generate_all --motifs 'TTTN<N{19}/.{4,6}\ ' --ontargetfilters DeepCpf1 --offtargetfilters Linear --weights Linear:85+1.7 --ko-gRNA --ko-dDNA mintag --revert_amplification_primers --fasta genome.fasta --gff genome.gff --folder GENEg > GENEg.out 2> GENEg.err
```

This writes 4 output tables to the `GENEg.out` file. Each of these tables refers to sequences in output FASTA files. Please note that certain sequence Aligners, such as 'Bowtie 2' can have non-deterministic output. Therefore, your results may vary from what is presented here.

Now would be a good time to explain the terminology you will see in the AddTag input and output. For simplicity in text processing, we use different labels than what are presented in the manuscript, though they are equivalent.
```
 r0-gDNA = +gDNA (wild type genome)
 r1-gDNA = ΔgDNA (intermediary genome)
 r2-gDNA = AgDNA (final genome)
 rN-gDNA = Nth round of genome engineering
exTarget = +Target (Target site in wild type Feature)
reTarget = ΔTarget (Target site introduced on ★tag insert)
 exDonor = r1-dDNA = ΔdDNA (ko-dDNA)
 reDonor = r2-dDNA = AdDNA (ki-dDNA)  
```

Thus, we refer to the first round of genome engineering (r1) as the knock-out round, and the second round (r2) as the knock-in round.

From the first table, we select the highest-weighed `reTarget` ('reversion Target' abbreviated), and one of its corresponding dDNA sequences. We store the sequence in its own FASTA file with the following command:
```sh
addtag find_header --fasta GENEg/dDNAs.fasta --query 'exDonor-33\b' > ko-dDNA.fasta
```
Each `reTarget` can target one or more identified `reDonor` dDNA sequences.

From the second table, we select the highest-weighted `exTarget` ('excision Target' abbreviated), which is used for excising the input Feature from the input gDNA.

Finally, we identify the highest-weight dDNA for reverting back to the wild type, and put it in its own FASTA file:
```sh
addtag find_header --fasta GENEg/dDNAs.fasta --query 'reDonor-0\b' > ko-dDNA.fasta
```


Next we need to identify a single cPCR verification primer design. Let's use the default pairwise sequence aligner.

```sh
addtag generate_primers --fasta genome.fasta --dDNAs ko-dDNA.fasta ki-dDNA.fasta --folder GENEc > GENEc.out 2> GENEc.err
``` 

### Typical workflow for multiplexed Features ###





## 📝 Citing AddTag ##
If you use AddTag for your research, please cite us. Because the manuscript is currently in preparation, you will need to cite the code repository instead.

 > Thaddeus D. Seher and Aaron D. Hernday. AddTag: Program for identifying exclusive endogenous gRNA sites and creating unique synthetic gRNA sites. University of California, Merced. Retrieved from <[https://github.com/tdseher/addtag-project](https://github.com/tdseher/addtag-project)> (2019).

## ✍ Authors ##
Who do I talk to?
 * Aaron D. Hernday (🔬 PI leading the project)
 * Thaddeus D. Seher (💻 programmer) (💬[@tdseher][tdseher])

See also the list of [contributors](https://github.com/tdseher/addtag-project/graphs/contributors) who participated in this project.

## 👥 Contributing ##
<details>
  <summary><h3>🐞 How do I submit a bug report?</h3></summary>

  First, check to see if the problem you are having has already been added to the [issue tracker](https://github.com/tdseher/addtag-project/issues).
  If not, then please submit a new issue.

  ▲
</details>

<details>
  <summary>⚠ How do I make a feature request?</summary>

  Send a message to [@tdseher][tdseher].
</details>

<details>
  <summary>⤴ How do I add my code to the AddTag software?</summary>

  Please submit a [pull request](https://github.com/tdseher/addtag-project/pulls).
</details>

<details>
  <summary>Adding scoring Algorithms</summary>

  Scoring Algorithms have been broken down into two general types.

   * `SingleSequenceAlgorithm` objects calculate scores by comparing a potential RNA ![Spacer][Spacer] or DNA ![Target][Target] to a model trained on empirical data.
   * `PairedSequenceAlgorithm` instances generate scores that compare a potential RNA ![Spacer][Spacer] to a DNA ![Target][Target].

  To add a new scoring algorithm, you must subclass one of the the above types, and add it to a `*.py` file in the `source/algorithms/` subdirectory. AddTag will automatically calculate the score on every generated ![Spacer][Spacer].

  We welcome any `git pull` requests to widen the repertoire of scoring algorithms available to AddTag. The easiest way to get started is to copy and modify one of the provided subclasses.
</details>

<details>
  <summary>Adding sequence Aligners</summary>

  AddTag comes with wrappers for several alignment programs. Depending on your experimental design and computing system, you may decide to use an aligner with no included wrapper. To implement your own, create a subclass of `Aligner`, and put it in a `*.py` file in the `source/aligners/` subdirectory. AddTag will automatically make that aligner available for you.

  Share your code with us so we can make it available to all AddTag users.
</details>

<details>
  <summary>Adding Thermodynamics calculators</summary>

  Several wrappers to popular oligonucleotide conformation, free energy, and melting temperature calculation programs are included. You can add your own by subclassing the `Oligo` class, and then adding its `*.py` file to the `source/thermodynamics/` subdirectory.

  If you create your own wrapper, please submit a `git pull` request so we can add it to the next version of the software.
</details>

## 📖 License ##
Please see the [LICENSE.md](LICENSE.md) file.

## Notes ##
<details>
  <summary>Below are tips and descriptions of AddTag limitations that will help you make successful designs.</summary>
  
   * The ![RGN][RGN] protein you use should be engineered specifically for your organism. It should be codon-optomized, and if using eukarya, contain an appropriate nuclear localization sequence.
   * By default, AddTag will avoid designing homology regions and Targets against polymorphisms whenever possible.
   * Sequences in FASTA files should have unique names. In other words, the primary sequence header--everything following the '`>`' character and preceding the first whitespace/tab '` `' character--should exist only once across all input `*.fasta` files.
   * AddTag makes no effort to restrict which Target motifs the user can use according to the selected Algorithms. Therefore, the user needs to independently verify which Target motifs are compatible with the selected Algorithms.
   * Right now AddTag can only handle linear chromosomes. If you want to analyze a circular chromosome, then you will need to artificially concatenate the ends of the chromosome together and adjust any annotations before running AddTag. An additional complication the software does not address is circular chromosomes. Features and their flanking regions cannot span the junction created when the contig end is concatenated to the start (typically the starting position on a contig is labeled the ORIGIN). To address this, the user should manually shift the coordinates of the experimental Features, and wrap the contigs as appropriate.
   * AddTag assumes one Feature copy per contig. The current implementation of AddTag assumes homology regions around Features are not repeated across any one contig. This means that is will fail to generate cPCR oligos for a large proportion of genes in transposon-rich genomes such as wheat (https://dx.doi.org/10.1186/s13059-018-1479-0).
   * A single feature cannot span two or more contigs. AddTag assumes that the entire feature sequence, and any flanking regions, are not in terminal regions of the reference contig. 
   * AddTag does not address overlapping genes, such as when an intron contains an exon for another gene, or when the same DNA encodes for genes on opposite strands. Everything between the Feature bounds is removed in the first engineering step. Currently, if the selected Feature overlaps with any other feature, only the selected Feature is considered. The other Feature will be disrupted. AddTag will report a warning that these other Features may be disrupted, but it does not attempt to reconcile this in any way. However, AddTag does have the ability to limit Feature expansion to keep the deletion outside of neighboring Features.
   * AddTag was not designed to perform paired Cas design, such as FokI-dCas9 nickase Users would need to run the program and select two gRNAs designed for opposite strands within a certain distance from each other. One way to mitigate errors is to use PAM-out nickases. This requires Cas9 cutting by two targets to get double-stranded break. This significantly decreases off-target genome editing. However, this initial AddTag version does not explicitly facilitate this.
   * AddTag can identify cut sites for Cas enzymes which have the PAM site. No functionality is provided for finding sites without an adjacent PAM sequence. AddTag requires motifs to define a PAM sequence. Therefore Cas14a is not supported. This can be probably be circumvented by using an `N` character as the PAM sequence, but this hasn't been tested. The number of CRISPR/Cas genome editing technologies are rapidly growing. With the recent discovery of Cas14a, which targets single-stranded DNA (ssDNA) molecules without requiring a PAM site (https://dx.doi.org/10.1126/science.aav4294), the expanded prevalence of CRISPR/Cas methods in biological sciences is assured. However, often researchers wish to edit sites on double-stranded DNA (dsDNA) using an RGN (such as Cas9 or Cas12a) that requires binding to a PAM motif. 
   * Please note, that at this time, no special restriction sites will be taken into account when designing primers.
   * For simplicity, all calculated scores ignore terms dealing with proximity to exon/CDS/ORF sequences. In cases such as the Stemmer and Azimuth calculations, the authors attempted to include the risk of disrupting genes neighboring potential targets in their models. We don’t attempt to do this.
   * Additionally, some scoring Algorithms take chromatin structure (DNA accessibility) into account. For simplicity, AddTag treats all input gDNA as equally accessible.
   * During the course of writing this software, a paper was published that outlines how hairpins can be inserted into the pre-spacer and spacer regions of the gRNA in order to increase specificity (https://dx.doi.org/10.1038/s41587-019-0095-1). AddTag does not model pre-spacer sequences.
   * AddTag assumes the RGN template type is dsDNA. AddTag was designed specifically to enable efficient gDNA editing. It does not use predictive models for ssDNA or RNA templates.
   * A corollary of this is that AddTag assumes all input sequences are DNA sequences. So the `--fasta` file specified will be treated as a DNA template. Thus, if there are any non-DNA residues, such as `U`, AddTag will probably fail. Also, since the Primer thermodynamics calculators are all set to estimate DNA:DNA hybridization (not DNA:RNA or RNA:RNA), any resulting calculations will be incorrect. 
   * Since Bartag motifs are user-specified, simple pre-computed lists of compatible 'bartag' sequences would be incomplete. Thus we implemented a greedy 'bartag' generation algorithm. When evaluating candidate 'bartag' sequences, AddTag will keep 'bartags' that satisfy all edit distance requirements with all previously-accepted 'bartags'. To limit runtime to a reasonable amount, we limited the total number of Features and 'bartags' that can be generated.
   * Of special note are things the Primer design does not explicitly consider, such as characteristics of the cPCR template molecule. AddTag does not exploit the differential nature of template sequence composition (e.g. H. sapiens compared to E. coli). Also, AddTag does not use information on the presence of known secondary modifications to the template, such as methylated residues or oxidative damage.
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
