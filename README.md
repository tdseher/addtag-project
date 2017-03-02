# CRISPR AddTag README #

Program for identifying unique endogenous gRNA sites and creating unique synthetic gRNA sites.

### Description ###

The CRISPR/Cas9 AddTag system can be used to do the following:
1. Generate gRNAs that target specific "features" uniquely (specified as GFF input), with low probabilities of off-target binding across the entire genome (specified as FASTA input).
2. Create "donor" DNA that contains a genome-wide unique gRNA-binding site with upstream and downstream homologous flanking regions, for use with complete gene knock-outs.
3. Construct conservative PCR primers for positive/negative amplification of "feature" knock-out.
4. Design gRNAs that target these unique "donor" DNA sites for Cas9 cutting.
5. Make "reversion" DNA sequences, composed of homologous sequences flanking the original "features" specific to each site, so the "feature" can be be re-inserted back into the genome. Please note, that at this time, no special restriction sites will be taken into account.
6. Compose conservative PCR primers for positive/negative amplification of "feature" reversion knock-in.

All generated sequences can be designed as either strand-specific or strand-agnostic. Additionally, generated sequences can target homologous chromosomes with the same gRNA or "donor" DNA if the input FASTA includes IUPAC ambiguity codes for polymorphisms. By default, AddTag will try to avoid all polymorphisms whenever possible.

### Requirements ###

* [Python 3.5+](https://www.python.org/downloads/)
* [regex module](https://pypi.python.org/pypi/regex)
* [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/)

### Who do I talk to? ###

* Aaron Hernday (PI leading the project)
* Thaddeus Seher (programmer)

### Notes ###

* The Cas9 protein you use should be engineered specifically for your organism. It should be codon-optomized, and contain a nuclear localization sequence for eukarya.

### References ###
* [Fu et al 2014](http://dx.doi.org/10.1038/nbt.2808) Using a shorter gRNAs (17-19 nt) can greatly improve specificity by reducing off-target binding
* [Braglia et al 2005](http://dx.doi.org/10.1074/jbc.M412238200) sequences containing consecutive Ts may cause polymerase termination
* [Ryan et al 2014](http://dx.doi.org/10.7554/eLife.03703) indicates that 50 bp flanking homology is sufficient to drive homologous recombination "donor" DNA knock-in.
* [Doench et al 2014](http://dx.doi.org/10.1038/nbt.3026) for scoring algorithm
* [Doench et al 2016](http://dx.doi.org/10.1038/nbt.3437) for scoring algorithm
* [Hsu et al 2013](http://dx.doi.org/10.1038/nbt.2647) for scoring algorithm
* [Haeussler et al 2016](http://dx.doi.org/10.1186/s13059-016-1012-2) CRISPOR
