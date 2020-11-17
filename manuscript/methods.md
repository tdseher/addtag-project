# Instructions #

This file describes the commands used to obtain the Targets, dDNAs, and Primers used in the manuscript. 

In order to proceed with the strain engineering portion of the manuscript in a reasonable time frame, we used an incomplete version of AddTag at the time.
We used AddTag r284 as the base version, which may have had uncommitted changes which could potentially alter the expected output.
Please note that the commands listed here are for an outdated version of AddTag.

We recommend you use the latest version of AddTag that has the full complement of software features.
The [README.md](../README.md) file contains these same workflows but updated for the current AddTag version.
Recent versions are less error-prone, are much more streamlined, and do not strictly require viewing of the `log.txt` file.

## Setup AddTag
Download and install AddTag and its dependencies. For instructions, please refer to the [README.md](../README.md) file located at the root of the Git repository.

## Get genome data
Download the *Candida albicans* reference genome and annotations used for this study.
```sh
wget http://www.candidagenome.org/download/sequence/C_albicans_SC5314/Assembly22/archive/C_albicans_SC5314_version_A22-s07-m01-r19_chromosomes.fasta.gz
gunzip C_albicans_SC5314_version_A22-s07-m01-r19_chromosomes.fasta.gz
wget http://www.candidagenome.org/download/gff/C_albicans_SC5314/archive/C_albicans_SC5314_version_A22-s07-m01-r19_features.gff
```

Set convenience variables for referencing these two files.
```sh
GENOME_FASTA=C_albicans_SC5314_version_A22-s07-m01-r19_chromosomes.fasta
GENOME_GFF=C_albicans_SC5314_version_A22-s07-m01-r19_features.gff
```

Create the `*.homologs` file for the *C. albicans* genome.
```sh
python3 gff2homologs.py ${GENOME_FEATURES} > C_albicans_SC5314_A22_current_homologs.txt
```

Obtain the general `*.homologs` file included with the package, or create one with only the genes you are interested in manipulating.
```sh
cat << EOF > C_albicans_SC5314_A22_current_homologs.txt
ADE2	C3_04520C_A	C3_04520C_B
BRG1	C1_05140W_A	C1_05140W_B
CSR1	C4_04850C_A	C4_04850C_B
EFG1	CR_07890W_A	CR_07890W_B
ROB1	C1_13620W_A	C1_13620W_B
ZRT2	C2_02590W_A	C2_02590W_B
EOF
```

## ADE2_CDS ##

For simplicity, we use a variable to hold the label for this computational experiment.
```sh
GENE=ADE2
```

Create and enter the directory for this experiment.
```sh
mkdir ${GENE}_CDS
cd ${GENE}_CDS
```

We create additional environmental variables for simplicity
```sh
DDNA_LEN=100
MIN_INS=0
MAX_INS=4
MIN_HOM=$(((DDNA_LEN-MAX_INS)/2)) # 48
MAX_HOM=$((DDNA_LEN/2)) # 50
```

Extract the annotations for ADE2 and put them in their own `*.gff` file.
```sh
addtag feature --gff ../${GENOME_FEATURES} --header --query ${GENE} > ${GENE}.gff
```

Identify the optimal Target sites and generate potential dDNAs.
```sh
addtag generate \
  --fasta ../${GENOME_FASTA} \
  --gff ${GENE}.gff \
  --homologs ../C_albicans_SC5314_A22_current_homologs.txt \
  --ko-gRNA \
  --ko-dDNA mintag \
  --ki-gRNA \
  --ki-dDNA \
  --features gene \
  --tag ID \
  --motifs 'N{17}|N{3}>NGG' \
  --off_target_motifs 'N{17}|N{3}>NAG' \
  --excise_upstream_homology ${MIN_HOM} ${MAX_HOM} \
  --excise_downstream_homology ${MIN_HOM} ${MAX_HOM} \
  --excise_donor_lengths ${DDNA_LEN} ${DDNA_LEN} \
  --excise_insert_lengths ${MIN_INS} ${MAX_INS} \
  --revert_amplification_primers \
  --max_homology_errors 0 \
  --revert_homology_length 100 200 \
  --primer_scan_limit 120 \
  --primer_pair_limit 120 \
  --folder ${GENE}g > ${GENE}g.out 2> ${GENE}g.err
```

We select the best +Target with no predicted off-targets (Hsu-Zhang = 100.0, and CFD = 100.0). This script can fail if no Targets meet this criteria.
```sh
TARGET=$(echo -e "import sys\nwith open(sys.argv[1], 'r') as flo:\n  for line in flo:\n    if not line.startswith('#'):\n      sline = line.split('\t')\n      if ((len(sline[1].split(',')) == 2) and sline[3].startswith('exTarget-') and (sline[5] == '100.0') and (sline[6] == '100.0')):\n        print(sline[3])\n        break\n" | python3 - ${GENE}g.out)
addtag extract \
  --fasta ${GENE}g/excision-spacers.fasta \
  --query "${TARGET}\b" > ${GENE}-ko-Target.fasta
```

Select the top-ranked ΔTarget and its associated ΔdDNA. When we ran it, `exDonor-55` and `reTarget-0` were optimal.
```sh
DONOR=exDonor-55
addtag extract --fasta ${GENE}g/excision-dDNAs.fasta --query "${DONOR}\b" > ${GENE}-ko-dDNA.fasta
TARGET=reTarget-0
addtag extract --fasta ${GENE}g/reversion-spacers.fasta --query "${TARGET}\b" > ${GENE}-ki-Target.fasta
```

Manually scan the output for the top-ranked AmpF/AmpR primers, and extract their associated AdDNA. We found `reDonor-0`.

> ```
> reDonor-0       C3_04520C_A     0.005368097146465493    PrimerPair(names=('AmpF', 'AmpR'), seqs=('CTAAGAAGGGAAAAGCACCAC', 'CTCGGTACAATCTTGTCAATGAG'), amplicon_sizes=[('ADE2', 'C3_04520C_A', 0, 'Ca22chr3A_C_albicans_SC5314', 1980), ('ADE2', 'C3_04520C_B', 0, 'Ca22chr3B_C_albicans_SC5314', 1980)], Tms=(52.77, 53.05), GCs=(0.48, 0.43), min(dG)=-3.72)
> ```

```sh
DONOR=reDonor-0
addtag extract --fasta ${GENE}g/reversion-dDNAs.fasta --query "${DONOR}\b" > ${GENE}-ki-dDNA.fasta
```

Calculate a decent Primer Design for validating each genome engineering step. We killed the process after getting through the first cycle.

```sh
addtag confirm 
  --fasta ../${GENOME_FASTA} \
  --dDNAs ${GENE}-ko-dDNA.fasta ${GENE}-ki-dDNA.fasta \
  --primer_scan_limit 600 \
  --primer_pair_limit 300 \
  --i_primers_required y n y \
  --o_primers_required y n y \
  --specificity all \
  --max_number_designs_reported 1000 \ 
  --folder ${GENE}c > ${GENE}c.out 2> ${GENE}c.err \
```

After the first cycle, a usable set of primers is identified.
Manually select the `PrimerSet` with the highest `weight`.
You can extract the primer sets with the following command.

```sh
grep -n 'cycle: 1' ${GENE}c/log.txt | cut -f 1 -d ':' | xargs -I {} head -n {} ${GENE}c/log.txt | grep optimal > optimal.txt
```
Then sort them by weight.

> ```
> optimal = PrimerSet(weight=9.515753054508352e-05, [PrimerPair(names=('', ''), seqs=('GTGGTGGATTGGTATTTCTTTCTGTG', 'AAGACCCCAAACATTTTGACTCG'), amplicon_sizes=[('exDonor-55,reDonor-0', 1, 2, 'Ca22chr3B_C_albicans_SC5314-r1[exDonor-55]-r2[reDonor-0]', 2046), ('exDonor-55,reDonor-0', 0, 0, 'Ca22chr3A_C_albicans_SC5314', 2046), ('exDonor-55,reDonor-0', 1, 1, 'Ca22chr3B_C_albicans_SC5314-r1[exDonor-55]', 341), ('exDonor-55,reDonor-0', 1, 0, 'Ca22chr3B_C_albicans_SC5314', 2046), ('exDonor-55,reDonor-0', 0, 1, 'Ca22chr3A_C_albicans_SC5314-r1[exDonor-55]', 341), ('exDonor-55,reDonor-0', 0, 2, 'Ca22chr3A_C_albicans_SC5314-r1[exDonor-55]-r2[reDonor-0]', 2046)], Tms=(55.99, 55.4), GCs=(0.42, 0.43), min(dG)=-3.34), PrimerPair(names=('', ''), seqs=('GTGGTGGATTGGTATTTCTTTCTGTG', 'CATTGCCTGTCATTGGTGTTCC'), amplicon_sizes=[('exDonor-55,reDonor-0', 1, 2, 'Ca22chr3B_C_albicans_SC5314-r1[exDonor-55]-r2[reDonor-0]', 474), ('exDonor-55,reDonor-0', 0, 0, 'Ca22chr3A_C_albicans_SC5314', 474), ('exDonor-55,reDonor-0', 1, 0, 'Ca22chr3B_C_albicans_SC5314', 474), ('exDonor-55,reDonor-0', 0, 2, 'Ca22chr3A_C_albicans_SC5314-r1[exDonor-55]-r2[reDonor-0]', 474)], Tms=(55.99, 56.26), GCs=(0.42, 0.5), min(dG)=-2.53), PrimerPair(names=('', ''), seqs=('CCCCAATGTGTAACAAGTCATCG', 'AAGACCCCAAACATTTTGACTCG'), amplicon_sizes=[('exDonor-55,reDonor-0', 1, 2, 'Ca22chr3B_C_albicans_SC5314-r1[exDonor-55]-r2[reDonor-0]', 567), ('exDonor-55,reDonor-0', 0, 0, 'Ca22chr3A_C_albicans_SC5314', 567), ('exDonor-55,reDonor-0', 1, 0, 'Ca22chr3B_C_albicans_SC5314', 567), ('exDonor-55,reDonor-0', 0, 2, 'Ca22chr3A_C_albicans_SC5314-r1[exDonor-55]-r2[reDonor-0]', 567)], Tms=(55.68, 55.4), GCs=(0.48, 0.43), min(dG)=-3.58), PrimerPair(names=('', ''), seqs=('CAGAGTTGTGAGGTCTTGGTG', 'GGCGTATGATGGTAGAGGTAAC'), amplicon_sizes=[('exDonor-55,reDonor-0', 1, 0, 'Ca22chr3B_C_albicans_SC5314', 377), ('exDonor-55,reDonor-0', 0, 2, 'Ca22chr3A_C_albicans_SC5314-r1[exDonor-55]-r2[reDonor-0]', 377), ('exDonor-55,reDonor-0', 1, 2, 'Ca22chr3B_C_albicans_SC5314-r1[exDonor-55]-r2[reDonor-0]', 377), ('exDonor-55,reDonor-0', 0, 0, 'Ca22chr3A_C_albicans_SC5314', 377)], Tms=(54.42, 54.05), GCs=(0.52, 0.5), min(dG)=-2.0), None, None, None, PrimerPair(names=('', ''), seqs=('GTGGTGGATTGGTATTTCTTTCTGTG', 'CATTGCCTGTCATTGGTGTTCC'), amplicon_sizes=[('exDonor-55,reDonor-0', 1, 2, 'Ca22chr3B_C_albicans_SC5314-r1[exDonor-55]-r2[reDonor-0]', 474), ('exDonor-55,reDonor-0', 0, 0, 'Ca22chr3A_C_albicans_SC5314', 474), ('exDonor-55,reDonor-0', 1, 0, 'Ca22chr3B_C_albicans_SC5314', 474), ('exDonor-55,reDonor-0', 0, 2, 'Ca22chr3A_C_albicans_SC5314-r1[exDonor-55]-r2[reDonor-0]', 474)], Tms=(55.99, 56.26), GCs=(0.42, 0.5), min(dG)=-2.53), PrimerPair(names=('', ''), seqs=('CCCCAATGTGTAACAAGTCATCG', 'AAGACCCCAAACATTTTGACTCG'), amplicon_sizes=[('exDonor-55,reDonor-0', 1, 2, 'Ca22chr3B_C_albicans_SC5314-r1[exDonor-55]-r2[reDonor-0]', 567), ('exDonor-55,reDonor-0', 0, 0, 'Ca22chr3A_C_albicans_SC5314', 567), ('exDonor-55,reDonor-0', 1, 0, 'Ca22chr3B_C_albicans_SC5314', 567), ('exDonor-55,reDonor-0', 0, 2, 'Ca22chr3A_C_albicans_SC5314-r1[exDonor-55]-r2[reDonor-0]', 567)], Tms=(55.68, 55.4), GCs=(0.48, 0.43), min(dG)=-3.58), PrimerPair(names=('', ''), seqs=('CAGAGTTGTGAGGTCTTGGTG', 'GGCGTATGATGGTAGAGGTAAC'), amplicon_sizes=[('exDonor-55,reDonor-0', 1, 0, 'Ca22chr3B_C_albicans_SC5314', 377), ('exDonor-55,reDonor-0', 0, 2, 'Ca22chr3A_C_albicans_SC5314-r1[exDonor-55]-r2[reDonor-0]', 377), ('exDonor-55,reDonor-0', 1, 2, 'Ca22chr3B_C_albicans_SC5314-r1[exDonor-55]-r2[reDonor-0]', 377), ('exDonor-55,reDonor-0', 0, 0, 'Ca22chr3A_C_albicans_SC5314', 377)], Tms=(54.42, 54.05), GCs=(0.52, 0.5), min(dG)=-2.0)])
> ```

Finally change back to the parent folder
```sh
cd ..
```

## EFG1 CDS ##

```sh
GENE=EFG1
cd ${GENE}_CDS

DDNA_LEN=100
MIN_INS=0
MAX_INS=4
MIN_HOM=$(((DDNA_LEN-MAX_INS)/2)) # 48
MAX_HOM=$((DDNA_LEN/2)) # 50
```

```sh
addtag feature --gff ../${GENOME_FEATURES} --header --query CR_07890W --linked_tags ID > EFG1.gff
```

```sh
addtag generate \
  --fasta ../${GENOME_FASTA} \
  --gff ${GENE}.gff \
  --homologs ../C_albicans_SC5314_A22_current_homologs.txt \
  --ko-gRNA \
  --ko-dDNA mintag \
  --ki-gRNA \
  --ki-dDNA \
  --features gene \
  --tag ID \
  --motifs 'N{17}|N{3}>NGG' \
  --off_target_motifs 'N{17}|N{3}>NAG' \
  --excise_upstream_homology ${MIN_HOM} ${MAX_HOM} \
  --excise_downstream_homology ${MIN_HOM} ${MAX_HOM} \
  --excise_donor_lengths ${DDNA_LEN} ${DDNA_LEN} \
  --excise_insert_lengths ${MIN_INS} ${MAX_INS} \
  --revert_amplification_primers \
  --excise_upstream_feature_trim 0 3 \
  --excise_downstream_feature_trim 7 10 \
  --max_homology_errors 0 \
  --revert_homology_length 100 200 \
  --primer_scan_limit 120 \
  --primer_pair_limit 120 \
  --folder ${GENE}g > ${GENE}g.out 2> ${GENE}g.err
```

```sh
addtag extract --fasta ${GENE}g/excision-spacers.fasta --query 'exTarget-97\b' > ${GENE}-ko-spacer.fasta
addtag extract --fasta ${GENE}g/reversion-spacers.fasta --query 'reTarget-118\b' > ${GENE}-ki-spacer.fasta
addtag extract --fasta ${GENE}g/excision-dDNAs.fasta --query 'exDonor-7183\b' > ${GENE}-ko-dDNA.fasta
addtag extract --fasta ${GENE}g/reversion-dDNAs.fasta --query 'reDonor-0\b' > ${GENE}-ki-dDNA.fasta
```

AmpF/AmpR
> ```
> reDonor-0       CR_07890W_A     0.002814545952972088    PrimerPair(names=('AmpF', 'AmpR'), seqs=('CACATTAGTTGCTCAGGTCAC', 'GTCAATGGATTTGGGAGAAGA'), amplicon_sizes=[('EFG1', 'CR_07890W_A', 0, 'Ca22chrRA_C_albicans_SC5314', 1972), ('EFG1', 'CR_07890W_B', 0, 'Ca22chrRB_C_albicans_SC5314', 1985)], Tms=(52.89, 52.04), GCs=(0.48, 0.43), min(dG)=-3.46)
> reDonor-349     CR_07890W_B     0.002814545952972088    PrimerPair(names=('AmpF', 'AmpR'), seqs=('CACATTAGTTGCTCAGGTCAC', 'GTCAATGGATTTGGGAGAAGA'), amplicon_sizes=[('EFG1', 'CR_07890W_A', 0, 'Ca22chrRA_C_albicans_SC5314', 1972), ('EFG1', 'CR_07890W_B', 0, 'Ca22chrRB_C_albicans_SC5314', 1985)], Tms=(52.89, 52.04), GCs=(0.48, 0.43), min(dG)=-3.46)
> ```

```sh
addtag confirm \
--fasta ../${GENOME_FASTA} \
--dDNAs ${GENE}-ko-dDNA.fasta ${GENE}-ki-dDNA.fasta \
--primer_scan_limit 120 \
--primer_pair_limit 120 \
--i_primers_required y n y \
--o_primers_required y n y \
--specificity all \
--max_number_designs_reported 1000 \
--cycle_start 10 \
--folder ${GENE}c > ${GENE}c.out 2> ${GENE}c.err
```

```sh
grep -n 'cycle: 11' ${GENE}c/log.txt | cut -f 1 -d ':' | xargs -I {} head -n {} ${GENE}c/log.txt | grep optimal > optimal.txt
```

cPCR primer set identified
> ```
> optimal = PrimerSet(weight=2.8889235120047e-11, [PrimerPair(names=('', ''), seqs=('TTAACCCCTTTGTGTCCCTT',    'CCCAAATAGTATAAATTCGTTCATGTC'), amplicon_sizes=[('exDonor-7183,reDonor-0', 1, 1, 'Ca22chrRB_C_albicans_SC5314-r1[exDonor-7183]', 575), ('exDonor-7183,reDonor-0', 1, 2, 'Ca22chrRB_C_albicans_SC5314-r1[exDonor-7183]-r2[reDonor-0]', 2236), ('exDonor-7183,reDonor-0', 0, 2, 'Ca22chrRA_C_albicans_SC5314-r1[exDonor-7183]-r2[reDonor-0]', 2236), ('exDonor-7183,reDonor-0', 0, 0, 'Ca22chrRA_C_albicans_SC5314', 2236), ('exDonor-7183,reDonor-0', 0, 1, 'Ca22chrRA_C_albicans_SC5314-r1[exDonor-7183]', 574), ('exDonor-7183,reDonor-0', 1, 0, 'Ca22chrRB_C_albicans_SC5314', 2249)], Tms=(52.59, 52.55), GCs=(0.45, 0.33), min(dG)=-3.06), PrimerPair(names=('', ''), seqs=('TTAACCCCTTTGTGTCCCTT',   'GCTGTTGTTGTTGTTGTCCT'), amplicon_sizes=[('exDonor-7183,reDonor-0', 1, 2, 'Ca22chrRB_C_albicans_SC5314-r1[exDonor-7183]-r2[reDonor-0]', 657), ('exDonor-7183,reDonor-0', 0, 2, 'Ca22chrRA_C_albicans_SC5314-r1[exDonor-7183]-r2[reDonor-0]', 657), ('exDonor-7183,reDonor-0', 0, 0, 'Ca22chrRA_C_albicans_SC5314', 657), ('exDonor-7183,reDonor-0', 1, 0, 'Ca22chrRB_C_albicans_SC5314', 658)], Tms=(52.59, 52.48), GCs=(0.45, 0.45), min(dG)=-2.3), PrimerPair(names=('', ''), seqs=('ACCAATCACCCCAAGTTCAG',     'CCCAAATAGTATAAATTCGTTCATGTC'), amplicon_sizes=[('exDonor-7183,reDonor-0', 0, 0, 'Ca22chrRA_C_albicans_SC5314', 301), ('exDonor-7183,reDonor-0', 1, 0, 'Ca22chrRB_C_albicans_SC5314', 301), ('exDonor-7183,reDonor-0', 1, 2, 'Ca22chrRB_C_albicans_SC5314-r1[exDonor-7183]-r2[reDonor-0]', 301), ('exDonor-7183,reDonor-0', 0, 2, 'Ca22chrRA_C_albicans_SC5314-r1[exDonor-7183]-r2[reDonor-0]', 301)], Tms=(53.12, 52.55), GCs=(0.5, 0.33), min(dG)=-3.06), PrimerPair(names=('', ''), seqs=('CCCCCATACCTTCCAATTCTAC',   'GACACATTACTGCCACCACTG'), amplicon_sizes=[('exDonor-7183,reDonor-0', 0, 0, 'Ca22chrRA_C_albicans_SC5314', 682), ('exDonor-7183,reDonor-0', 1, 0, 'Ca22chrRB_C_albicans_SC5314', 682), ('exDonor-7183,reDonor-0', 1, 2, 'Ca22chrRB_C_albicans_SC5314-r1[exDonor-7183]-r2[reDonor-0]', 682), ('exDonor-7183,reDonor-0', 0, 2, 'Ca22chrRA_C_albicans_SC5314-r1[exDonor-7183]-r2[reDonor-0]', 682)], Tms=(53.82, 55.19), GCs=(0.5, 0.52), min(dG)=-1.59), None, None, None, PrimerPair(names=('', ''), seqs=('TTAACCCCTTTGTGTCCCTT', 'GCTGTTGTTGTTGTTGTCCT'), amplicon_sizes=[('exDonor-7183,reDonor-0', 1, 2, 'Ca22chrRB_C_albicans_SC5314-r1[exDonor-7183]-r2[reDonor-0]', 657), ('exDonor-7183,reDonor-0', 0, 2, 'Ca22chrRA_C_albicans_SC5314-r1[exDonor-7183]-r2[reDonor-0]', 657), ('exDonor-7183,reDonor-0', 0, 0, 'Ca22chrRA_C_albicans_SC5314', 657), ('exDonor-7183,reDonor-0', 1, 0, 'Ca22chrRB_C_albicans_SC5314', 658)], Tms=(52.59, 52.48), GCs=(0.45, 0.45), min(dG)=-2.3), PrimerPair(names=('', ''), seqs=('ACCAATCACCCCAAGTTCAG', 'CCCAAATAGTATAAATTCGTTCATGTC'), amplicon_sizes=[('exDonor-7183,reDonor-0', 0, 0, 'Ca22chrRA_C_albicans_SC5314', 301), ('exDonor-7183,reDonor-0', 1, 0, 'Ca22chrRB_C_albicans_SC5314', 301), ('exDonor-7183,reDonor-0', 1, 2, 'Ca22chrRB_C_albicans_SC5314-r1[exDonor-7183]-r2[reDonor-0]', 301), ('exDonor-7183,reDonor-0', 0, 2, 'Ca22chrRA_C_albicans_SC5314-r1[exDonor-7183]-r2[reDonor-0]', 301)], Tms=(53.12, 52.55), GCs=(0.5, 0.33), min(dG)=-3.06), PrimerPair(names=('', ''), seqs=('CCCCCATACCTTCCAATTCTAC', 'GACACATTACTGCCACCACTG'), amplicon_sizes=[('exDonor-7183,reDonor-0', 0, 0, 'Ca22chrRA_C_albicans_SC5314', 682), ('exDonor-7183,reDonor-0', 1, 0, 'Ca22chrRB_C_albicans_SC5314', 682), ('exDonor-7183,reDonor-0', 1, 2, 'Ca22chrRB_C_albicans_SC5314-r1[exDonor-7183]-r2[reDonor-0]', 682), ('exDonor-7183,reDonor-0', 0, 2, 'Ca22chrRA_C_albicans_SC5314-r1[exDonor-7183]-r2[reDonor-0]', 682)], Tms=(53.82, 55.19), GCs=(0.5, 0.52), min(dG)=-1.59)])
> ```


## BRG1_CDS ##

```sh
GENE=BRG1
cd ${GENE}_CDS

DDNA_LEN=100
MIN_INS=0
MAX_INS=4
MIN_HOM=$(((DDNA_LEN-MAX_INS)/2)) # 48
MAX_HOM=$((DDNA_LEN/2)) # 50
```

Include only the BRG1 features.
```sh
addtag feature --gff ../${GENOME_FEATURES} --header --query ${GENE} > ${GENE}.gff
```

With this version of AddTag, there was a bug in the feature expansion, so it is corrected here.
```sh
awk '{$4-=42}1' FS='\t' OFS='\t' ${GENE}.gff > ${GENE}ex.gff
```

Run AddTag to find Targets, make dDNAs, and find AmpF/AmpR primers.
```sh
addtag generate \
  --fasta ../${GENOME_FASTA} \
  --gff ${GENE}ex.gff \
  --homologs ../C_albicans_SC5314_A22_current_homologs.txt \
  --ko-gRNA \
  --ko-dDNA mintag \
  --ki-gRNA \
  --ki-dDNA \
  --features gene \
  --tag ID \
  --motifs 'N{17}|N{3}>NGG' \
  --off_target_motifs 'N{17}|N{3}>NAG' \
  --excise_upstream_homology ${MIN_HOM} ${MAX_HOM} \
  --excise_downstream_homology ${MIN_HOM} ${MAX_HOM} \
  --excise_donor_lengths ${DDNA_LEN} ${DDNA_LEN} \
  --excise_insert_lengths ${MIN_INS} ${MAX_INS} \
  --revert_amplification_primers \
  --max_homology_errors 0 \
  --revert_homology_length 100 200 \
  --primer_scan_limit 120 \
  --primer_pair_limit 120 \
  --folder ${GENE}g > ${GENE}g.out 2> ${GENE}g.err
```

```sh
TARGET=$(echo -e "import sys\nwith open(sys.argv[1], 'r') as flo:\n  for line in flo:\n    if not line.startswith('#'):\n      sline = line.split('\t')\n      if ((len(sline[1].split(',')) == 2) and sline[3].startswith('exTarget-') and (sline[5] == '100.0') and (sline[6] == '100.0')):\n        print(sline[3])\n        break\n" | python3 - ${GENE}g.out)
addtag extract \
  --fasta ${GENE}g/excision-spacers.fasta \
  --query "${TARGET}\b" > ${GENE}-ko-Target.fasta
```

```sh
addtag extract --fasta ${GENE}g/excision-dDNAs.fasta --query 'exDonor-440\b' > ${GENE}-ko-dDNA.fasta
addtag extract --fasta ${GENE}g/reversion-spacers.fasta --query 'reTarget-628\b' > ${GENE}-ki-spacer.fasta
addtag extract --fasta ${GENE}g/reversion-dDNAs.fasta --query 'reDonor-80\b' > ${GENE}-ki-dDNA.fasta
```

A single optimal primer pair was predicted to create one amplicon for each homologous chromosome.
> ```
> reDonor-0       C1_05140W_B     4.45941697919701e-05    PrimerPair(names=('AmpF', 'AmpR'), seqs=('TATAAATATCAGGTCATAGATCCCTG', 'GTTCAAACAACAATACTGTAGCAG'), amplicon_sizes=[('BRG1', 'C1_05140W_B', 0, 'Ca22chr1B_C_albicans_SC5314', 1608), ('BRG1', 'C1_05140W_A', 0, 'Ca22chr1A_C_albicans_SC5314', 1605)], Tms=(52.03, 52.39), GCs=(0.35, 0.38), min(dG)=-3.44)
> reDonor-80      C1_05140W_A     4.45941697919701e-05    PrimerPair(names=('AmpF', 'AmpR'), seqs=('TATAAATATCAGGTCATAGATCCCTG', 'GTTCAAACAACAATACTGTAGCAG'), amplicon_sizes=[('BRG1', 'C1_05140W_B', 0, 'Ca22chr1B_C_albicans_SC5314', 1608), ('BRG1', 'C1_05140W_A', 0, 'Ca22chr1A_C_albicans_SC5314', 1605)], Tms=(52.03, 52.39), GCs=(0.35, 0.38), min(dG)=-3.44)
> ```

Identify the set of primers to use for cPCR validation

```sh 
addtag confirm \
  --fasta ../${GENOME_FASTA} \
  --dDNAs ${GENE}-ko-dDNA.fasta ${GENE}-ki-dDNA.fasta \
  --primer_scan_limit 3600 \
  --primer_pair_limit 600 \
  --i_primers_required y n y \
  --o_primers_required y n y \
  --specificity all \
  --max_number_designs_reported 1000 \
  --cycle_start 6 
  --folder ${GENE}c > ${GENE}c.out 2> ${GENE}c.err
```

After 6 cycles, a usable set of primers was identified.
```ah
grep -n 'cycle: 7' ${GENE}c/log.txt | cut -f 1 -d ':' | xargs -I {} head -n {} ${GENE}c/log.txt | grep -c optimal > optimal.txt
```

> ```
> optimal = PrimerSet(weight=6.792474801162299e-10, [PrimerPair(names=('', ''), seqs=('TGCAGCTTTTGTACTACATTTGG', 'CCAGCTCAGGATATAATTTACAGC'), amplicon_sizes=[('exDonor-440,reDonor-80', 0, 1, 'Ca22chr1B_C_albicans_SC5314-r1[exDonor-440]', 636), ('exDonor-440,reDonor-80', 1, 0, 'Ca22chr1A_C_albicans_SC5314', 1950), ('exDonor-440,reDonor-80', 1, 2, 'Ca22chr1A_C_albicans_SC5314-r1[exDonor-440]-r2[reDonor-80]', 1950), ('exDonor-440,reDonor-80', 0, 2, 'Ca22chr1B_C_albicans_SC5314-r1[exDonor-440]-r2[reDonor-80]', 1947), ('exDonor-440,reDonor-80', 0, 0, 'Ca22chr1B_C_albicans_SC5314', 1950), ('exDonor-440,reDonor-80', 1, 1, 'Ca22chr1A_C_albicans_SC5314-r1[exDonor-440]', 639)], Tms=(53.24, 52.99), GCs=(0.39, 0.42), min(dG)=-3.26), PrimerPair(names=('', ''), seqs=('TGCAGCTTTTGTACTACATTTGG', 'ACCTCCACTAATGGTTGATCG'), amplicon_sizes=[('exDonor-440,reDonor-80', 1, 0, 'Ca22chr1A_C_albicans_SC5314', 521), ('exDonor-440,reDonor-80', 1, 2, 'Ca22chr1A_C_albicans_SC5314-r1[exDonor-440]-r2[reDonor-80]', 521), ('exDonor-440,reDonor-80', 0, 2, 'Ca22chr1B_C_albicans_SC5314-r1[exDonor-440]-r2[reDonor-80]', 518), ('exDonor-440,reDonor-80', 0, 0, 'Ca22chr1B_C_albicans_SC5314', 521)], Tms=(53.24, 52.84), GCs=(0.39, 0.48), min(dG)=-3.3), PrimerPair(names=('', ''), seqs=('GTCATTCATCAACCACCACCA', 'CCAGCTCAGGATATAATTTACAGC'), amplicon_sizes=[('exDonor-440,reDonor-80', 1, 0, 'Ca22chr1A_C_albicans_SC5314', 337), ('exDonor-440,reDonor-80', 1, 2, 'Ca22chr1A_C_albicans_SC5314-r1[exDonor-440]-r2[reDonor-80]', 337), ('exDonor-440,reDonor-80', 0, 2, 'Ca22chr1B_C_albicans_SC5314-r1[exDonor-440]-r2[reDonor-80]', 337), ('exDonor-440,reDonor-80', 0, 0, 'Ca22chr1B_C_albicans_SC5314', 337)], Tms=(54.01, 52.99), GCs=(0.48, 0.42), min(dG)=-3.05), PrimerPair(names=('', ''), seqs=('CCACCACAACAACCACAATCAG', 'CGACCGTTCTTCCCTTTTGTC'), amplicon_sizes=[('exDonor-440,reDonor-80', 1, 2, 'Ca22chr1A_C_albicans_SC5314-r1[exDonor-440]-r2[reDonor-80]', 569), ('exDonor-440,reDonor-80', 1, 0, 'Ca22chr1A_C_albicans_SC5314', 569), ('exDonor-440,reDonor-80', 0, 2, 'Ca22chr1B_C_albicans_SC5314-r1[exDonor-440]-r2[reDonor-80]', 569), ('exDonor-440,reDonor-80', 0, 0, 'Ca22chr1B_C_albicans_SC5314', 569)], Tms=(56.07, 55.5), GCs=(0.5, 0.52), min(dG)=-1.87), None, None, None, PrimerPair(names=('', ''), seqs=('TGCAGCTTTTGTACTACATTTGG', 'ACCTCCACTAATGGTTGATCG'), amplicon_sizes=[('exDonor-440,reDonor-80', 1, 0, 'Ca22chr1A_C_albicans_SC5314', 521), ('exDonor-440,reDonor-80', 1, 2, 'Ca22chr1A_C_albicans_SC5314-r1[exDonor-440]-r2[reDonor-80]', 521), ('exDonor-440,reDonor-80', 0, 2, 'Ca22chr1B_C_albicans_SC5314-r1[exDonor-440]-r2[reDonor-80]', 518), ('exDonor-440,reDonor-80', 0, 0, 'Ca22chr1B_C_albicans_SC5314', 521)], Tms=(53.24, 52.84), GCs=(0.39, 0.48), min(dG)=-3.3), PrimerPair(names=('', ''), seqs=('GTCATTCATCAACCACCACCA', 'CCAGCTCAGGATATAATTTACAGC'), amplicon_sizes=[('exDonor-440,reDonor-80', 1, 0, 'Ca22chr1A_C_albicans_SC5314', 337), ('exDonor-440,reDonor-80', 1, 2, 'Ca22chr1A_C_albicans_SC5314-r1[exDonor-440]-r2[reDonor-80]', 337), ('exDonor-440,reDonor-80', 0, 2, 'Ca22chr1B_C_albicans_SC5314-r1[exDonor-440]-r2[reDonor-80]', 337), ('exDonor-440,reDonor-80', 0, 0, 'Ca22chr1B_C_albicans_SC5314', 337)], Tms=(54.01, 52.99), GCs=(0.48, 0.42), min(dG)=-3.05), PrimerPair(names=('', ''), seqs=('CCACCACAACAACCACAATCAG', 'CGACCGTTCTTCCCTTTTGTC'), amplicon_sizes=[('exDonor-440,reDonor-80', 1, 2, 'Ca22chr1A_C_albicans_SC5314-r1[exDonor-440]-r2[reDonor-80]', 569), ('exDonor-440,reDonor-80', 1, 0, 'Ca22chr1A_C_albicans_SC5314', 569), ('exDonor-440,reDonor-80', 0, 2, 'Ca22chr1B_C_albicans_SC5314-r1[exDonor-440]-r2[reDonor-80]', 569), ('exDonor-440,reDonor-80', 0, 0, 'Ca22chr1B_C_albicans_SC5314', 569)], Tms=(56.07, 55.5), GCs=(0.5, 0.52), min(dG)=-1.87)])
> ```


Change back to the parent folder
```sh
cd ..
```

## CSR1_US (ZAP1_US) ##
Obtain or create the `*.gff` annotation for the CSR1 self-regulator.
```sh
cat << EOF > ZRT2_US.gff
# seqid	source	feature	start	end	score	strand	frame	attribute
Ca22chr4A_C_albicans_SC5314	Seher	protein_bind	1044578	1044588	.	-	.	ID=CSR1_US_A
Ca22chr4B_C_albicans_SC5314	Seher	protein_bind	1044606	1044616	.	-	.	ID=CSR1_US_B
EOF
```

Obtain or create the `*.homologs` file for this new annotation.
```sh
cat << EOF > CSR1_US_homologs.txt
CSR1_US	CSR1_US_A	CSR1_US_B
EOF
```

Mask the Csr1 transcription factor binding site by changing the characters to lower case, and save as `C_albicans_SC5314_A22_current_chromosomes_CSR1_US.fasta`

Make a FASTA file with the modified binding site
```sh
cat << EOF > CSR1_US_MOD.fasta
>CSR1_US
CAGCTGTTCTA
EOF
```

Make a FASTA file with the addtag sequence
```sh
cat << EOF > CSR1_US_UNITAG.fasta
>CSR1_US
CGTACGCTGCAGGTCGACAGTGG
EOF
```

The generated ΔdDNA are longer to accomodate the 23 nt addtag sequence.
```sh
GENE=CSR1_US
cd ${GENE}

DDNA_LEN=130
```

```sh
addtag generate \
  --fasta ../C_albicans_SC5314_A22_current_chromosomes_CSR1_US.fasta \
  --gff ${GENE}.gff ../${GENOME_FEATURES} \
  --homologs ${GENE}_homologs.txt \
  --ko-gRNA \
  --ko-dDNA ${GENE}_UNITAG.fasta \
  --ki-gRNA \
  --ki-dDNA ${GENE}_MOD.fasta \
  --features protein_bind \
  --excluded_features gene \
  --tag ID \
  --motifs 'N{17}|N{3}>NGG' \
  --off_target_motifs 'N{17}|N{3}>NAG' \
  --excise_donor_lengths ${DDNA_LEN} ${DDNA_LEN} \
  --feature_expansion_method center_both \
  --feature_expansion_pad 200 \
  --feature_expansion_lengths 300 500 \
  --revert_amplification_primers \
  --max_homology_errors 0 \
  --revert_homology_length 100 200 \
  --primer_scan_limit 120 \
  --primer_pair_limit 120 \
  --folder ${GENE}g > ${GENE}g.out 2> ${GENE}g.err
```

> ```
> reDonor-0       CSR1_US_B_derived-0     0.03497651705045208     PrimerPair(names=('AmpF', 'AmpR'), seqs=('GTGAACCACTCATCATCATTGG', 'ACGACTAATGCTATGACTGCTC'), amplicon_sizes=[('CSR1_US', 'CSR1_US_A', 0, 'Ca22chr4A_C_albicans_SC5314', 641), ('CSR1_US', 'CSR1_US_B', 0, 'Ca22chr4B_C_albicans_SC5314', 642)], Tms=(53.25, 53.32), GCs=(0.45, 0.45), min(dG)=-3.47)
> reDonor-809     CSR1_US_A_derived-0     0.03497651705045208     PrimerPair(names=('AmpF', 'AmpR'), seqs=('GTGAACCACTCATCATCATTGG', 'ACGACTAATGCTATGACTGCTC'), amplicon_sizes=[('CSR1_US', 'CSR1_US_A', 0, 'Ca22chr4A_C_albicans_SC5314', 641), ('CSR1_US', 'CSR1_US_B', 0, 'Ca22chr4B_C_albicans_SC5314', 642)], Tms=(53.25, 53.32), GCs=(0.45, 0.45), min(dG)=-3.47)
> ```

```sh
addtag extract --fasta ${GENE}g/reversion-spacers.fasta --query 'reTarget-0\b' > ${GENE}-ki-spacer.fasta
addtag extract --fasta ${GENE}g/excision-spacers.fasta --query 'exTarget-5\b' > ${GENE}-ko-spacer.fasta
addtag extract --fasta ${GENE}g/reversion-dDNAs.fasta --query 'reDonor-809\b' > ${GENE}-ki-dDNA.fasta
addtag extract --fasta ${GENE}g/excision-dDNAs.fasta --query 'exDonor-0\b' > ${GENE}-ko-dDNA.fasta
```

```
addtag confirm \
  --fasta ../${GENOME_FASTA} \
  --dDNAs ${GENE}-ko-dDNA.fasta ${GENE}-ki-dDNA.fasta \
  --primer_scan_limit 120 \
  --primer_pair_limit 120 \
  --i_primers_required y n y \
  --o_primers_required y y y \
  --specificity all \
  --max_number_designs_reported 1000 \
  --cycle_start 13 \
  --folder ${GENE}c > ${GENE}c.out 2> ${GENE}c.err
```

The `{$GENE}c/log.txt` contains the optimal primer set.
> ```
> optimal = PrimerSet(weight=1.1089802250429206e-12, [PrimerPair(names=('', ''), seqs=('CTGTGATCGTGATTATGAATGTGGC', 'ACGTTGTTCGTCTCAAGCTGG'), amplicon_sizes=[('exDonor-0,reDonor-809', 0, 2, 'Ca22chr4B_C_albicans_SC5314-r1[exDonor-0]-r2[reDonor-809]', 978), ('exDonor-0,reDonor-809', 0, 1, 'Ca22chr4B_C_albicans_SC5314-r1[exDonor-0]', 578), ('exDonor-0,reDonor-809', 1, 2, 'Ca22chr4A_C_albicans_SC5314-r1[exDonor-0]-r2[reDonor-809]', 977), ('exDonor-0,reDonor-809', 1, 0, 'Ca22chr4A_C_albicans_SC5314', 974), ('exDonor-0,reDonor-809', 0, 0, 'Ca22chr4B_C_albicans_SC5314', 977), ('exDonor-0,reDonor-809', 1, 1, 'Ca22chr4A_C_albicans_SC5314-r1[exDonor-0]', 577)], Tms=(56.19, 56.41), GCs=(0.44, 0.52), min(dG)=-3.44), PrimerPair(names=('', ''), seqs=('CTGTGATCGTGATTATGAATGTGGC', 'TACCCAGCATCATCATCATCG'), amplicon_sizes=[('exDonor-0,reDonor-809', 0, 2, 'Ca22chr4B_C_albicans_SC5314-r1[exDonor-0]-r2[reDonor-809]', 442), ('exDonor-0,reDonor-809', 1, 2, 'Ca22chr4A_C_albicans_SC5314-r1[exDonor-0]-r2[reDonor-809]', 442), ('exDonor-0,reDonor-809', 1, 0, 'Ca22chr4A_C_albicans_SC5314', 439), ('exDonor-0,reDonor-809', 0, 0, 'Ca22chr4B_C_albicans_SC5314', 442)], Tms=(56.19, 55.05), GCs=(0.44, 0.48), min(dG)=-2.94), PrimerPair(names=('', ''), seqs=('GTTGTCGATGATGATGATGCTGG', 'ACGTTGTTCGTCTCAAGCTGG'), amplicon_sizes=[('exDonor-0,reDonor-809', 1, 0, 'Ca22chr4A_C_albicans_SC5314', 561), ('exDonor-0,reDonor-809', 0, 0, 'Ca22chr4B_C_albicans_SC5314', 561), ('exDonor-0,reDonor-809', 0, 2, 'Ca22chr4B_C_albicans_SC5314-r1[exDonor-0]-r2[reDonor-809]', 562), ('exDonor-0,reDonor-809', 1, 2, 'Ca22chr4A_C_albicans_SC5314-r1[exDonor-0]-r2[reDonor-809]', 561)], Tms=(55.92, 56.41), GCs=(0.48, 0.52), min(dG)=-3.44), PrimerPair(names=('', ''), seqs=('CGATGATGATGATGCTGGGTA', 'GTGTTACTTGGTAGCACTTTGATC'), amplicon_sizes=[('exDonor-0,reDonor-809', 0, 0, 'Ca22chr4B_C_albicans_SC5314', 319), ('exDonor-0,reDonor-809', 0, 2, 'Ca22chr4B_C_albicans_SC5314-r1[exDonor-0]-r2[reDonor-809]', 320), ('exDonor-0,reDonor-809', 1, 2, 'Ca22chr4A_C_albicans_SC5314-r1[exDonor-0]-r2[reDonor-809]', 320), ('exDonor-0,reDonor-809', 1, 0, 'Ca22chr4A_C_albicans_SC5314', 320)], Tms=(55.05, 53.68), GCs=(0.48, 0.42), min(dG)=-4.24), PrimerPair(names=('', ''), seqs=('CTGTGATCGTGATTATGAATGTGGC', 'TGGTGTTACTTGGTAGCCCAC'), amplicon_sizes=[('exDonor-0,reDonor-809', 0, 1, 'Ca22chr4B_C_albicans_SC5314-r1[exDonor-0]', 343), ('exDonor-0,reDonor-809', 1, 1, 'Ca22chr4A_C_albicans_SC5314-r1[exDonor-0]', 343)], Tms=(56.19, 55.41), GCs=(0.44, 0.52), min(dG)=-3.82), PrimerPair(names=('', ''), seqs=('AGAAAGTGGCGTTTAATAAATACTACGTAC', 'ACGTTGTTCGTCTCAAGCTGG'), amplicon_sizes=[('exDonor-0,reDonor-809', 0, 1, 'Ca22chr4B_C_albicans_SC5314-r1[exDonor-0]', 300), ('exDonor-0,reDonor-809', 1, 1, 'Ca22chr4A_C_albicans_SC5314-r1[exDonor-0]', 299)], Tms=(55.13, 56.41), GCs=(0.33, 0.52), min(dG)=-4.58), None, PrimerPair(names=('', ''), seqs=('CTGTGATCGTGATTATGAATGTGGC', 'TACCCAGCATCATCATCATCG'), amplicon_sizes=[('exDonor-0,reDonor-809', 0, 2, 'Ca22chr4B_C_albicans_SC5314-r1[exDonor-0]-r2[reDonor-809]', 442), ('exDonor-0,reDonor-809', 1, 2, 'Ca22chr4A_C_albicans_SC5314-r1[exDonor-0]-r2[reDonor-809]', 442), ('exDonor-0,reDonor-809', 1, 0, 'Ca22chr4A_C_albicans_SC5314', 439), ('exDonor-0,reDonor-809', 0, 0, 'Ca22chr4B_C_albicans_SC5314', 442)], Tms=(56.19, 55.05), GCs=(0.44, 0.48), min(dG)=-2.94), PrimerPair(names=('', ''), seqs=('GTTGTCGATGATGATGATGCTGG', 'ACGTTGTTCGTCTCAAGCTGG'), amplicon_sizes=[('exDonor-0,reDonor-809', 1, 0, 'Ca22chr4A_C_albicans_SC5314', 561), ('exDonor-0,reDonor-809', 0, 0, 'Ca22chr4B_C_albicans_SC5314', 561), ('exDonor-0,reDonor-809', 0, 2, 'Ca22chr4B_C_albicans_SC5314-r1[exDonor-0]-r2[reDonor-809]', 562), ('exDonor-0,reDonor-809', 1, 2, 'Ca22chr4A_C_albicans_SC5314-r1[exDonor-0]-r2[reDonor-809]', 561)], Tms=(55.92, 56.41), GCs=(0.48, 0.52), min(dG)=-3.44), PrimerPair(names=('', ''), seqs=('CGATGATGATGATGCTGGGTA', 'GTGTTACTTGGTAGCACTTTGATC'), amplicon_sizes=[('exDonor-0,reDonor-809', 0, 0, 'Ca22chr4B_C_albicans_SC5314', 319), ('exDonor-0,reDonor-809', 0, 2, 'Ca22chr4B_C_albicans_SC5314-r1[exDonor-0]-r2[reDonor-809]', 320), ('exDonor-0,reDonor-809', 1, 2, 'Ca22chr4A_C_albicans_SC5314-r1[exDonor-0]-r2[reDonor-809]', 320), ('exDonor-0,reDonor-809', 1, 0, 'Ca22chr4A_C_albicans_SC5314', 320)], Tms=(55.05, 53.68), GCs=(0.48, 0.42), min(dG)=-4.24)])
> ```

## ZRT2 US ##

Obtain or create the `*.gff` annotations for the region that includes both CSR1 binding sites, and everything in between.
```sh
cat << EOF > ZRT2_US.gff
# seqid	source	feature	start	end	score	strand	frame	attribute
Ca22chr2A_C_albicans_SC5314	Seher	protein_bind	523078	523744	.	+	.	ID=ZRT2_US_A
Ca22chr2B_C_albicans_SC5314	Seher	protein_bind	523078	523743	.	+	.	ID=ZRT2_US_B
```

Obtain or create the `*.homologs` file for this new annotation.
```sh
cat << EOF > ZRT2_US_homologs.txt
ZRT2_US	ZRT2_US_A	ZRT2_US_B
EOF
```

Make a FASTA file with the addtag sequence
```sh
cat << EOF > ZRT2_US_UNITAG.fasta
>ZRT2_US
CGTACGCTGCAGGTCGACAGTGG
EOF
```

```sh
GENE=ZRT2_US
cd ${GENE}

DDNA_LEN=130
```

```sh
addtag generate \
  --fasta ../${GENOME_FASTA} \
  --gff ${GENE}.gff ../${GENOME_FEATURES} \
  --homologs ${GENE}_homologs.txt \
  --ko-gRNA \
  --ko-dDNA ${GENE}_UNITAG.fasta \
  --ki-gRNA \
  --ki-dDNA ${GENE}_MOD-00-A.fasta \
  --features protein_bind \
  --excluded_features gene \
  --tag ID \
  --motifs 'N{17}|N{3}>NGG' \
  --off_target_motifs 'N{17}|N{3}>NAG' \
  --excise_donor_lengths ${DDNA_LEN} ${DDNA_LEN} \
  --revert_amplification_primers \
  --max_homology_errors 0 \
  --case upper-only \
  --revert_homology_length 100 200 \
  --primer_scan_limit 120 \
  --primer_pair_limit 120 \
  --folder ${GENE}g > ${GENE}g.out 2> ${GENE}g.err
```

> ```
> reDonor-0	ZRT2_US_B	0.0021891124190610164	PrimerPair(names=('AmpF', 'AmpR'), seqs=('GCATATTTACTTGCTTGCCTG', 'TTGACAGGAATATGGAGGGTA'), amplicon_sizes=[('ZRT2_US', 'ZRT2_US_B', 0, 'Ca22chr2B_C_albicans_SC5314', 953), ('ZRT2_US', 'ZRT2_US_A', 0, 'Ca22chr2A_C_albicans_SC5314', 955)], Tms=(52.02, 52.71), GCs=(0.43, 0.43), min(dG)=-4.99)
> reDonor-407	ZRT2_US_A	0.0021891124190610164	PrimerPair(names=('AmpF', 'AmpR'), seqs=('GCATATTTACTTGCTTGCCTG', 'TTGACAGGAATATGGAGGGTA'), amplicon_sizes=[('ZRT2_US', 'ZRT2_US_B', 0, 'Ca22chr2B_C_albicans_SC5314', 953), ('ZRT2_US', 'ZRT2_US_A', 0, 'Ca22chr2A_C_albicans_SC5314', 955)], Tms=(52.02, 52.71), GCs=(0.43, 0.43), min(dG)=-4.99)
> ```

```sh
addtag extract --fasta ${GENE}g/reversion-spacers.fasta --query 'reTarget-0\b' > ${GENE}-ki-spacer.fasta
addtag extract --fasta ${GENE}g/excision-spacers.fasta --query 'exTarget-11\b' > ${GENE}-ko-spacer.fasta
addtag extract --fasta ${GENE}g/reversion-dDNAs.fasta --query 'reDonor-407\b' > ${GENE}-ki-dDNA.fasta
addtag extract --fasta ${GENE}g/excision-dDNAs.fasta --query 'exDonor-0\b' > ${GENE}-ko-dDNA.fasta
```

```sh
addtag confirm \
  --fasta ../${GENOME_FASTA} \
  --dDNAs ${GENE}-ko-dDNA.fasta ${GENE}-ki-dDNA.fasta \
  --case upper-only \
  --primer_scan_limit 120 \
  --primer_pair_limit 120 \
  --i_primers_required y n y \
  --o_primers_required y n y \
  --specificity all \
  --max_number_designs_reported 1000 \
  --folder ${GENE}c > ${GENE}c.out 2> ${GENE}c.err
```

The `{$GENE}c/log.txt` contains the optimal primer set.
> ```
> optimal = PrimerSet(weight=9.355670319407343e-06, [PrimerPair(names=('', ''), seqs=('GAACCAATCCTTCCACATAGC', 'GCTGGGAATTGATAATGAAAGC'), amplicon_sizes=[('exDonor-0,reDonor-407', 1, 0, 'Ca22chr2B_C_albicans_SC5314', 1059), ('exDonor-0,reDonor-407', 1, 2, 'Ca22chr2B_C_albicans_SC5314-r1[exDonor-0]-r2[reDonor-407]', 1061), ('exDonor-0,reDonor-407', 0, 0, 'Ca22chr2A_C_albicans_SC5314', 1061), ('exDonor-0,reDonor-407', 0, 2, 'Ca22chr2A_C_albicans_SC5314-r1[exDonor-0]-r2[reDonor-407]', 1061), ('exDonor-0,reDonor-407', 1, 1, 'Ca22chr2B_C_albicans_SC5314-r1[exDonor-0]', 416), ('exDonor-0,reDonor-407', 0, 1, 'Ca22chr2A_C_albicans_SC5314-r1[exDonor-0]', 417)], Tms=(52.7, 52.04), GCs=(0.48, 0.41), min(dG)=-3.45), PrimerPair(names=('', ''), seqs=('GAACCAATCCTTCCACATAGC', 'TTGCGTTTCGGGTATAATCAC'), amplicon_sizes=[('exDonor-0,reDonor-407', 1, 0, 'Ca22chr2B_C_albicans_SC5314', 654), ('exDonor-0,reDonor-407', 1, 2, 'Ca22chr2B_C_albicans_SC5314-r1[exDonor-0]-r2[reDonor-407]', 655), ('exDonor-0,reDonor-407', 0, 0, 'Ca22chr2A_C_albicans_SC5314', 655), ('exDonor-0,reDonor-407', 0, 2, 'Ca22chr2A_C_albicans_SC5314-r1[exDonor-0]-r2[reDonor-407]', 655)], Tms=(52.7, 52.81), GCs=(0.48, 0.43), min(dG)=-3.13), PrimerPair(names=('', ''), seqs=('TATTGGTCGGATTGGGTTAC', 'GCTGGGAATTGATAATGAAAGC'), amplicon_sizes=[('exDonor-0,reDonor-407', 0, 0, 'Ca22chr2A_C_albicans_SC5314', 369), ('exDonor-0,reDonor-407', 0, 2, 'Ca22chr2A_C_albicans_SC5314-r1[exDonor-0]-r2[reDonor-407]', 369), ('exDonor-0,reDonor-407', 1, 0, 'Ca22chr2B_C_albicans_SC5314', 368), ('exDonor-0,reDonor-407', 1, 2, 'Ca22chr2B_C_albicans_SC5314-r1[exDonor-0]-r2[reDonor-407]', 369)], Tms=(52.36, 52.04), GCs=(0.45, 0.41), min(dG)=-1.9), PrimerPair(names=('', ''), seqs=('GAGAAGAACCATAAAGTCCAAGC', 'CACCTCAAACCACACACTAC'), amplicon_sizes=[('exDonor-0,reDonor-407', 1, 2, 'Ca22chr2B_C_albicans_SC5314-r1[exDonor-0]-r2[reDonor-407]', 337), ('exDonor-0,reDonor-407', 1, 0, 'Ca22chr2B_C_albicans_SC5314', 337), ('exDonor-0,reDonor-407', 0, 2, 'Ca22chr2A_C_albicans_SC5314-r1[exDonor-0]-r2[reDonor-407]', 337), ('exDonor-0,reDonor-407', 0, 0, 'Ca22chr2A_C_albicans_SC5314', 337)], Tms=(53.16, 52.36), GCs=(0.43, 0.5), min(dG)=-2.71), None, None, None, PrimerPair(names=('', ''), seqs=('GAACCAATCCTTCCACATAGC', 'TTGCGTTTCGGGTATAATCAC'), amplicon_sizes=[('exDonor-0,reDonor-407', 1, 0, 'Ca22chr2B_C_albicans_SC5314', 654), ('exDonor-0,reDonor-407', 1, 2, 'Ca22chr2B_C_albicans_SC5314-r1[exDonor-0]-r2[reDonor-407]', 655), ('exDonor-0,reDonor-407', 0, 0, 'Ca22chr2A_C_albicans_SC5314', 655), ('exDonor-0,reDonor-407', 0, 2, 'Ca22chr2A_C_albicans_SC5314-r1[exDonor-0]-r2[reDonor-407]', 655)], Tms=(52.7, 52.81), GCs=(0.48, 0.43), min(dG)=-3.13), PrimerPair(names=('', ''), seqs=('TATTGGTCGGATTGGGTTAC', 'GCTGGGAATTGATAATGAAAGC'), amplicon_sizes=[('exDonor-0,reDonor-407', 0, 0, 'Ca22chr2A_C_albicans_SC5314', 369), ('exDonor-0,reDonor-407', 0, 2, 'Ca22chr2A_C_albicans_SC5314-r1[exDonor-0]-r2[reDonor-407]', 369), ('exDonor-0,reDonor-407', 1, 0, 'Ca22chr2B_C_albicans_SC5314', 368), ('exDonor-0,reDonor-407', 1, 2, 'Ca22chr2B_C_albicans_SC5314-r1[exDonor-0]-r2[reDonor-407]', 369)], Tms=(52.36, 52.04), GCs=(0.45, 0.41), min(dG)=-1.9), PrimerPair(names=('', ''), seqs=('GAGAAGAACCATAAAGTCCAAGC', 'CACCTCAAACCACACACTAC'), amplicon_sizes=[('exDonor-0,reDonor-407', 1, 2, 'Ca22chr2B_C_albicans_SC5314-r1[exDonor-0]-r2[reDonor-407]', 337), ('exDonor-0,reDonor-407', 1, 0, 'Ca22chr2B_C_albicans_SC5314', 337), ('exDonor-0,reDonor-407', 0, 2, 'Ca22chr2A_C_albicans_SC5314-r1[exDonor-0]-r2[reDonor-407]', 337), ('exDonor-0,reDonor-407', 0, 0, 'Ca22chr2A_C_albicans_SC5314', 337)], Tms=(53.16, 52.36), GCs=(0.43, 0.5), min(dG)=-2.71)])
> ```


## FAQ ##
### Why do overall PCR primer stringency levels differ for each gene? ###
The AddTag software tried all primers at all higher stringencies (lower numbers), and didn’t find a viable design until it reached the listed number. To save computation time, we had the program skip up to 6, 10, and 13 for BRG1, EFG1, and CSR1.

### Why are you filtering the whole-genome GFF? ###
Early versions of AddTag needed abridged GFF as input, otherwise they would attempt to design Targets and dDNAs for all Features of the selected type. Later versions of AddTag introduced the functionality to list the specific feature IDs you wish to include.

