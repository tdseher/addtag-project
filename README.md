# CRISPR/Cas9 AddTag README #

Program for identifying exclusive endogenous gRNA sites and creating unique synthetic gRNA sites.

### Description ###

The CRISPR/Cas9 AddTag system can be used to do the following:

 1. Generate gRNAs that target specific "features" uniquely (specified as GFF input), with low probabilities of off-target binding across the entire genome (specified as FASTA input).
 2. Create donor DNA (dDNA) sequences that contains a genome-wide unique gRNA-binding site with upstream and downstream homologous flanking regions, for use with complete gene knock-outs (excision).
 3. Design gRNAs that target these unique dDNA sites for Cas9 cutting.
 4. Make primers for "reversion" dDNA sequence amplification, with homologous sequences flanking the original "feature" specific to each site, so the "feature" can be be re-inserted back into the genome. Please note, that at this time, no special restriction sites will be taken into account.
 5. Construct conservative PCR primers for positive amplification of "feature" knock-out (excision) and knock-in (reversion).

All generated sequences can be designed as either strand-specific or strand-agnostic. Additionally, generated sequences can target homologous chromosomes with the same gRNA or "donor" DNA if the input FASTA includes IUPAC ambiguity codes for polymorphisms. By default, AddTag will try to avoid all polymorphisms whenever possible.

### Requirements ###

 * Python >= 3.5.1 ([source](https://www.python.org/downloads/), [binaries](https://www.python.org/downloads/), [documentation](https://docs.python.org/3/))
    - regex module ([source](https://bitbucket.org/mrabarnett/mrab-regex), [whls](https://pypi.org/project/regex/), [documentation](https://pypi.org/project/regex/))
      
        note: The easy way to install this is through `pip`.
        
        ```
        $ pip3 install regex
        ```

 * 3.0.0 > Python >= 2.7.10 ([source](https://www.python.org/downloads/), [binaries](https://www.python.org/downloads/), [documentation](https://docs.python.org/2/))
    - Azimuth module ([source](https://github.com/MicrosoftResearch/Azimuth), [documentation](https://www.microsoft.com/en-us/research/project/azimuth/))
      
        note: requires python-tk to be installed.
        
        note: requires specific versions of scipy, numpy, matplotlib, nose, scikit-learn, pandas, biopython, pyparsing, cycler, six, pytz, python-dateutil, functools32, subprocess32. The easiest way to install it and take care of all dependencies is to use `pip`. This process assumes `git` is available in the `PATH` environmental variable.
        
        ```
        $ pip2.7 install git+https://github.com/MicrosoftResearch/Azimuth.git
        ```

 * For speed, we recommend at least one third-party pairwise nucleotide sequence alignment program.
    - BLAST+ >= 2.6.0 ([source](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST), [binaries](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST), [documentation](https://www.ncbi.nlm.nih.gov/books/NBK279690/))
    - BLAT ([source](https://genome.ucsc.edu/goldenPath/help/blatSpec.html), [binaries](http://hgdownload.cse.ucsc.edu/admin/exe/), [documentation](https://genome.ucsc.edu/goldenPath/help/blatSpec.html))
    - Bowtie ([source](https://github.com/BenLangmead/bowtie), [binaries](https://sourceforge.net/projects/bowtie-bio/files/bowtie/), [documentation](http://bowtie-bio.sourceforge.net/manual.shtml))
    - Bowtie 2 >= 2.3.4.1 ([source](https://github.com/BenLangmead/bowtie2), [binaries](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/), [documentation](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml))
    - BWA ([source](https://github.com/lh3/bwa), [documentation](http://bio-bwa.sourceforge.net/bwa.shtml))
    - Cas-OFFinder ([source](https://github.com/snugel/cas-offinder), [binaries](https://sourceforge.net/projects/cas-offinder/files/Binaries/), [documentation](http://www.rgenome.net/cas-offinder/portable))

 * For cPCR oligo design, we recommend at least one third-party thermodynamics calculation program.
    - UNAFold ([documentation](http://unafold.rna.albany.edu/))
    - primer3-py module ([source](https://github.com/libnano/primer3-py), [whls](https://pypi.org/project/primer3-py/), [documentation](https://libnano.github.io/primer3-py/))
      
        note: The easy way to install this is through `pip`.
        
        ```
        $ pip3 install primer3-py
        ```

 * Git >= 1.7.1 ([source, binaries](https://git-scm.com/downloads)) Certain AddTag functionality depends on Git.
 * If you are using the Windows operating system, before running AddTag, it is recommended to set the environmental variable `PYTHONIOENCODING` to `utf-8`, otherwise, STDOUT may not be written correctly to files.
   
    ```> set PYTHONIOENCODING=utf-8```

### Obtaining AddTag ###
You can download the latest version of AddTag over HTTPS using `git` with the following command (replacing `username` with your BitBucket account name).
```sh
$ git clone https://username@bitbucket.org/tdseher/addtag-project.git
```

This will ask you for your Atlassian password associated with the `username`.

Or you can download AddTag over SSH without typing your password if the repository is configured with your 'publickey'.
```sh
$ git clone git@bitbucket.org:tdseher/addtag-project.git
```

Either of these options will download the AddTag into a folder called `addtag-project/` in your current working directory. Go ahead and change the working directory into the AddTag folder.
```sh
$ cd addtag-project/
```

`git` should automatically make the `addtag` program executable. If it does not, you can use the following command to do it. 
```sh
$ chmod +x addtag
```

### Updating AddTag ###
Set the working directory as the AddTag folder.
```sh
$ cd addtag-project/
```

If you would like to update your local copy to the newest version available, use the following command from within the `addtag-project/` directory.
```sh
$ git pull
```

If you want the newest version, but you made local changes, then you can first discard your changes, and then update. Use the following two commands from inside the `addtag-project/` folder.
```sh
$ git reset --hard
$ git pull
```

Alternatively, if you want to keep the local modifications, you can use stash to hide them away before pulling, then reapply them afterwards.
```sh
$ git stash
$ git pull
$ git stash pop
```

Each one of these methods may require your Atlassian login credentials.

### Program usage ###
The following commands assume the current working directory is the AddTag folder `addtag-project/`. To view the program usage, you may run AddTag with the `--help` flag. This will print out command line parameter descriptions and examples.
```sh
$ ./addtag --help
```

Additionally, you may view the included man page.
```sh
$ man ./addtag.1
```

### Scoring algorithms ###
Scoring algorithms have been broken down into two general types.

 * `SingleSequenceAlgorithm` objects calculate scores by comparing a potential spacer to a model trained on empirical data (Azimuth, Chari, Doench-2014, GC, Housden, Moreno-Mateos, Wang)
 * `PairedSequenceAlgorithm` instances generate scores that compare a potential spacer to a target using a model (CFD, Distance, Hsu-Zhang, Linear)

Adding a new scoring algorithm is as simple as subclassing the above type, and adding it to a `*.py` file in the `source/algorithms/` subdirectory. AddTag will automatically calculate the score on every generated spacer.

We welcome any pull requests to widen the repertoire of scoring algorithms available to AddTag. The easiest way to get started is to copy and modify one of the provided subclasses.

### Alignment programs ###
AddTag comes with wrappers for several alignment programs. Depending on your experimental design and computing system, you may decide to use an aligner with no included wrapper. To implement your own, create a subclass of `Aligner`, and put it in a `*.py` file in the `source/aligners/` subdirectory. AddTag will automatically make that aligner available for you.

Share your code with us so we can make it available to all AddTag users.

### Thermodynamics calculations ###
Several wrappers to popular oligonucleotide conformation, free energy, and melting temperature calculation programs are included. You can add your own by subclassing the `Oligo` class, and then adding its `*.py` file to the `source/oligos` subdirectory.

If you create your own wrapper, please submit a pull request so we can add it to the next version of the software.

### Who do I talk to? ###

 * Aaron Hernday (PI leading the project)
 * Thaddeus Seher (programmer)

### Notes ###

 * The Cas9 protein you use should be engineered specifically for your organism. It should be codon-optomized, and if using eukarya, contain an appropriate nuclear localization sequence.

### References ###
 * [Nguyen, Quail, & Hernday (2017)](http://dx.doi.org/10.1128/mSphereDirect.00149-17) for AddTag concept.
 * [Moreno-Mateos, et al (2015)](http://dx.doi.org/10.1038/nmeth.3543) for scoring algorithm
 * [Lin, et al (2014)](http://dx.doi.org/10.1093/nar/gku402) The GC content of Cas9 target may affects binding specificity; gRNA may bind off-target if it has insertions/deletions (RNA-bulge/DNA-bulge) relative to multiple genome locations.
 * [Fu, et al (2014)](http://dx.doi.org/10.1038/nbt.2808) Using a shorter gRNAs (17-19 nt) can greatly improve specificity by reducing off-target binding
 * [Vyas, et al (2015)](http://dx.doi.org/10.1126/sciadv.1500248) Anecdotally, gRNA may target sites less efficiently if they have differences within 12 nt of the PAM.
 * [Braglia, et al (2005)](http://dx.doi.org/10.1074/jbc.M412238200) sequences containing consecutive Ts may cause polymerase termination
 * [Ronda, et al (2015)](http://dx.doi.org/10.1186/s12934-015-0288-3) 60 bp flanking homology is sufficient to drive HR in *Saccharomyces cerevisiae*.
 * [Ryan, et al (2014)](http://dx.doi.org/10.7554/eLife.03703) indicates that 50 bp flanking homology is sufficient to drive homologous recombination "donor" DNA knock-in.
 * [Doench, et al (2014)](http://dx.doi.org/10.1038/nbt.3026) for scoring algorithm
 * [Doench, et al (2016)](http://dx.doi.org/10.1038/nbt.3437) for scoring algorithm
 * [Hsu, et al (2013)](http://dx.doi.org/10.1038/nbt.2647) for scoring algorithm
 * [Haeussler, et al (2016)](http://dx.doi.org/10.1186/s13059-016-1012-2) CRISPOR paper for implementation of certain scoring algorithms
