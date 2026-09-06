codon_harmonization_tool

A lightweight Python desktop application for codon usage analysis and score-based codon harmonization between a source and a target organism.

The tool compares codon usage patterns between two organisms, calculates a normalized codon usage score for every synonymous codon, analyzes a coding sequence codon by codon, and generates a new target coding sequence in which each codon is selected to preserve the relative codon usage preference of the source sequence as closely as possible.

The application includes a graphical user interface built with tkinter and does not require external Python packages.

Overview

Codon usage differs between organisms. Even though several codons may encode the same amino acid, these synonymous codons are not necessarily used with equal frequency.

codon_harmonization_tool uses codon frequency tables from a source organism and a target organism to calculate a normalized codon score for each codon.

For each codon in the source coding sequence, the program then selects a synonymous codon in the target organism whose normalized usage score is as close as possible to the source score.

The resulting coding sequence therefore aims to preserve the relative codon usage profile of the original sequence rather than simply replacing every codon with the most frequent codon in the target organism.

The tool can additionally detect inconsistencies between an entered amino acid sequence and the supplied source coding sequence. If such a mismatch is present, the amino acid sequence is treated as the desired protein sequence and the generated target codon is selected from the codons encoding the intended amino acid.

Main features

Graphical user interface using Python tkinter

Input of:

amino acid sequence

source coding sequence

source-organism codon frequency table

target-organism codon frequency table

Direct parsing of Kazusa/Kazuza-style codon usage tables

Automatic conversion of RNA codons (U) to DNA codons (T)

Calculation of normalized codon usage scores

Calculation of Relative Frequency Numbers (RFN)

Codon-by-codon comparison of source and target usage

Score-based codon harmonization

Detection and correction of amino acid / coding-sequence mismatches

Generation of a harmonized target coding sequence

Export of:

codon-by-codon analysis as CSV

generated target CDS

calculated source and target codon metric tables

No external Python dependencies

How the algorithm works

1. Codon frequency

The input codon usage tables contain a frequency for each codon.

For example:

UUU 26.1(170666)
UUC 18.4(120510)

The program reads:

UUU -> 26.1
UUC -> 18.4

The numbers in parentheses are ignored.

Internally, RNA codons are converted to DNA notation:

UUU -> TTT
UUC -> TTC

The absolute scale of the frequency values is not important for the score calculation as long as all synonymous codons within a table are reported using the same scale.

2. Codon score

For each amino acid, the frequencies of all synonymous codons are summed.

The score of an individual codon is then calculated as:

Codon score =
frequency of the codon
------------------------------------------
sum of frequencies of all synonymous codons

For example, if phenylalanine is encoded by:

TTT = 26.1
TTC = 18.4

then:

Total = 26.1 + 18.4 = 44.5

Score(TTT) = 26.1 / 44.5 = 0.5865
Score(TTC) = 18.4 / 44.5 = 0.4135

The score therefore represents the relative usage of a codon among all codons encoding the same amino acid.

For amino acids encoded by only one codon, such as methionine (ATG) or tryptophan (TGG), the score is normally:

1.0

provided that a non-zero frequency is present in the input table.

3. Relative Frequency Number (RFN)

The program also calculates a Relative Frequency Number (RFN).

For each amino acid, synonymous codons are sorted by increasing score:

lowest score  -> RFN 1
next score    -> RFN 2
...
highest score -> highest RFN

RFN is useful for visualizing the relative rarity of codons.

Example:

Alanine codons:

GCG  score 0.10 -> RFN 1
GCT  score 0.20 -> RFN 2
GCC  score 0.30 -> RFN 3
GCA  score 0.40 -> RFN 4

The current version of the program does not use RFN for codon harmonization.

RFN is retained as an analysis parameter and for rare-codon visualization.

4. Rare-codon threshold

The GUI contains a configurable:

Rare-Codon-Schwelle (RFN <= ...)

With the default value:

RFN <= 1

the least frequently used synonymous codon for each amino acid is classified as rare.

Increasing the threshold to 2, for example, includes the two lowest-ranked synonymous codons.

This classification is used for the analysis and graphical highlighting only.

It does not influence the harmonization algorithm.

Score-based codon harmonization

The central harmonization step is based on score similarity, not RFN matching.

For every position in the source coding sequence:

The source codon is identified.

Its source-organism codon score is determined.

The desired amino acid for that position is identified.

All target-organism codons encoding that amino acid are considered.

The program calculates the absolute score difference:

Delta score = |target score - source score|

The target codon with the smallest score difference is selected.

In other words:

best target codon =
synonymous target codon minimizing
|TargetScore - SourceScore|

This attempts to preserve the relative codon usage behavior of the source sequence in the target organism.

Example

Assume that a source codon has:

Source codon: TCT
Amino acid: Serine
Source score: 0.18

and the target organism has the following serine codons:

TCT  score 0.31
TCC  score 0.19
TCA  score 0.22
TCG  score 0.05
AGT  score 0.15
AGC  score 0.08

The program compares all score differences:

|0.31 - 0.18| = 0.13
|0.19 - 0.18| = 0.01
|0.22 - 0.18| = 0.04
|0.05 - 0.18| = 0.13
|0.15 - 0.18| = 0.03
|0.08 - 0.18| = 0.10

Therefore:

TCC

is selected because its target score (0.19) is closest to the source score (0.18).

Amino acid / CDS matchmaking

The tool can use both:

an amino acid sequence

a source coding sequence

The amino acid sequence acts as the desired protein sequence.

The source CDS provides the codon usage information that should be transferred to the target organism.

For each sequence position, the software compares:

amino acid translated from source codon
vs.
amino acid supplied in the amino acid sequence

If they match

The program performs normal score-based harmonization.

If they do not match

The entered amino acid sequence is treated as authoritative.

The program:

keeps the source codon's original source score,

identifies all target codons encoding the intended amino acid,

chooses the target codon whose target score is closest to the original source score.

Example:

Desired amino acid: S
Source codon: TAA

TAA is a stop codon and therefore does not encode serine.

The program nevertheless uses the source score of TAA as the relative usage reference and searches only among the six serine codons in the target organism:

TCT
TCC
TCA
TCG
AGT
AGC

The serine codon with the closest target score is selected for the generated target sequence.

This allows the tool to correct nucleotide sequences while retaining as much of the original codon-usage pattern as possible.

Tie-breaking

If two or more target codons have exactly the same score distance from the source score, the current implementation uses the following tie-breaking order:

smallest absolute score difference

highest target codon frequency

alphabetical codon order

This makes codon selection deterministic.

Requirements

Python 3.9 or newer

tkinter

No third-party Python packages are required.

The program uses only modules from the Python standard library:

csv
io
re
tkinter

Windows

tkinter is normally included with the standard Python installer from python.org.

Linux

Depending on the distribution, Tk may have to be installed separately.

For example, on Debian/Ubuntu:

sudo apt install python3-tk

Running the program

Clone or download the repository and run:

python codon_harmonization_tool.py

On some systems:

python3 codon_harmonization_tool.py

The graphical interface will open automatically.

Input

1. Amino acid sequence

Enter a protein sequence using the one-letter amino acid code.

Example:

MAFGQRTCDNVPELGIQVAFTEYGLNSDAFPKVYGLDNRASGINVPEY

A terminal stop may be represented as:

*

Example:

MAFGQRTCDNVPELGIQVAFTEY*

Spaces and formatting characters are removed automatically.

If an amino acid sequence is supplied, its length must correspond to the number of codons in the source coding sequence.

2. Source coding sequence

Enter the coding nucleotide sequence from the source organism.

Example:

ATGGCTTTCGGCCAGCGTCGTACTTGTGATAACGTTCCGGAATTAGGCATTCAGGTT

Spaces and line breaks are accepted:

ATG GCT TTC GGC CAG CGT CGT ACT
TGT GAT AAC GTT CCG GAA

RNA bases are also accepted and automatically converted:

AUG -> ATG

The cleaned nucleotide sequence must have a length divisible by three.

Codon frequency table formats

Kazusa/Kazuza-style table

The program can directly parse blocks copied in formats such as:

UUU 26.1(170666)  UCU 23.5(153557)  UAU 18.8(122728)  UGU 8.1(52903)
UUC 18.4(120510)  UCC 14.2(92923)   UAC 14.8(96596)   UGC 4.8(31095)
UUA 26.2(170884)  UCA 18.7(122028)  UAA 1.1(6913)     UGA 0.7(4447)
UUG 27.2(177573)  UCG 8.6(55951)    UAG 0.5(3312)     UGG 10.4(67789)

For every entry:

UUU 26.1(170666)

the application uses:

codon = UUU
frequency = 26.1

and ignores:

(170666)

Simple text format

TTT 26.1
TTC 18.4
TTA 26.2
TTG 27.2

CSV format

A CSV file can contain:

codon,frequency
TTT,26.1
TTC,18.4
TTA,26.2
TTG,27.2

Accepted codon column names include:

codon
triplet
codons

Accepted frequency column names include:

frequency
freq
per_thousand
value
count

TSV format

Tab-separated tables are also supported:

codon	frequency
TTT	26.1
TTC	18.4

Source and target codon tables

Two frequency tables are required:

Source organism

Represents the codon usage distribution of the organism from which the original coding sequence originates.

Target organism

Represents the codon usage distribution of the organism in which the harmonized sequence should be expressed.

For reliable results, complete 64-codon tables are strongly recommended.

Output

After clicking:

Analyse starten

the program generates a codon-by-codon analysis.

The result table contains:

Column

Description

Position

Codon position in the sequence

AA

Desired amino acid at this position

Source Codon

Original codon from the source CDS

Source Frequency

Frequency of the source codon in the source organism

Source Score

Normalized synonymous codon usage score in the source organism

Source RFN

Relative frequency rank in the source organism

Target Codon

Selected codon for the target organism

Target Frequency

Frequency of the selected codon in the target organism

Target Score

Normalized synonymous codon usage score in the target organism

Target RFN

Relative frequency rank in the target organism

Note

Information about matching, score distance, and possible AA/CDS mismatch

Summary output

The result page also reports:

total number of codons

number of rare codons according to the selected RFN threshold

number of amino-acid/CDS mismatches

amino acid sequence translated from the source CDS

desired amino acid sequence

generated target coding sequence

translated target amino acid sequence

Export options

Analysis as CSV

Exports the complete codon-by-codon analysis.

This is useful for downstream analysis in:

Excel

Origin

R

Python

statistical software

Target CDS

The generated harmonized coding sequence can be exported as:

plain text

FASTA

Example FASTA output:

>target_codon_sequence
ATGGCTTTC...

RFN / codon metric tables

The calculated source and target metrics can also be exported.

The table includes:

organism
codon
aa
frequency
score
rfn
synonymous_total

Interpretation of the main metrics

Frequency

The codon frequency value supplied by the codon usage table.

For example:

GCT = 21.2

may represent approximately 21.2 occurrences per 1000 codons, depending on the source database.

Score

Normalized codon usage within one synonymous amino acid group.

Score =
codon frequency /
sum of frequencies of all synonymous codons

A score close to 1 means that the codon represents a large fraction of codons used for that amino acid.

A low score means that the codon is relatively uncommon compared with synonymous alternatives.

RFN

The rank of a synonymous codon after sorting by increasing score.

RFN 1 = lowest-scoring synonymous codon

RFN is descriptive in the current harmonization workflow and is not used to choose the target codon.

Why score-based harmonization?

Simple codon optimization commonly favors highly frequent codons in the target organism.

This tool follows a different principle.

Instead of maximizing codon frequency, it attempts to preserve the relative synonymous codon usage pattern of the source sequence.

Conceptually:

source rare-ish codon
        |
        v
target codon with similar relative usage

source common codon
        |
        v
target codon with similar relative usage

This makes the program suitable for exploring codon harmonization rather than conventional maximal-frequency codon optimization.

Important limitations

1. Codon usage is not the only determinant of translation

The program uses codon-frequency relationships only.

It does not model:

tRNA abundance directly

tRNA charging

wobble pairing

ribosome pausing

codon-pair bias

dicodon effects

mRNA secondary structure

GC content constraints

RNA stability

transcriptional regulation

translation initiation

internal ribosome binding sites

cryptic splice sites

repetitive sequence motifs

restriction sites

RNA degradation motifs

organism-specific regulatory signals

The generated sequence should therefore be interpreted as a codon-usage harmonization candidate, not as a guaranteed optimal expression construct.

2. Complete codon tables are strongly recommended

If a codon is missing from an input table, the current implementation assigns:

frequency = 0

This can influence:

codon scores

RFN values

target matching

For biological analyses, use complete codon usage tables whenever possible.

3. Score similarity is a heuristic

The program minimizes:

|TargetScore - SourceScore|

This is a simple and transparent harmonization strategy.

It does not establish that two codons with equal normalized frequencies necessarily produce identical translational effects in different organisms.

4. Standard genetic code

The application currently uses the standard genetic code.

Organisms or organelles using alternative genetic codes are not supported without modifying GENETIC_CODE.

5. Stop codons

Stop codons are represented as:

*

and are treated as a synonymous group containing:

TAA
TAG
TGA

If the sequence contains a stop codon, the target stop codon is selected using the same score-similarity principle.

6. Amino acid sequence length

When an amino acid sequence is provided, it must contain exactly one amino acid symbol per source codon.

For example:

60 codons -> 60 amino acid symbols

Otherwise, the analysis is stopped to avoid positional misalignment.

Recommended workflow

Obtain a complete codon usage table for the source organism.

Obtain a complete codon usage table for the target organism.

Paste the source protein sequence into Aminosäuresequenz.

Paste the original coding sequence into Coding-Sequenz im Source-Organismus.

Paste the source codon usage table into Source-Organismus.

Paste the target codon usage table into Target-Organismus.

Select the desired RFN threshold for rare-codon visualization.

Click Analyse starten.

Review:

source scores

target scores

score differences

mismatch notes

Export the harmonized target CDS and the analysis table.

Example conceptual workflow

Source protein sequence
        +
Source coding sequence
        +
Source codon usage table
        +
Target codon usage table
        |
        v
Calculate source and target codon scores
        |
        v
Analyze each source codon
        |
        v
Determine desired amino acid
        |
        v
Find synonymous target codon with
minimum |TargetScore - SourceScore|
        |
        v
Assemble harmonized target CDS

Project structure

A minimal repository can be organized as:

codon_harmonization_tool/
│
├── codon_harmonization_tool.py
├── README.md
└── LICENSE

Optional additions:

examples/
screenshots/
test_data/

Suggested repository description

GUI-based Python tool for codon usage analysis and score-based codon harmonization between source and target organisms using codon frequency tables.

Scientific use

The tool was developed as a lightweight research utility for comparing organism-specific codon usage and generating candidate coding sequences for heterologous expression workflows.

It is particularly intended for:

exploratory codon-usage analysis

codon harmonization

comparison of source and target organisms

heterologous protein-expression workflow development

research and teaching

For experimental applications, generated sequences should be evaluated together with additional sequence-design criteria and validated experimentally.

License

No license is included automatically with the script.

If the repository should be openly reusable, add an appropriate license file, for example:

MIT License

BSD 3-Clause License

GNU GPLv3

Choose the license according to the intended use and institutional requirements.

Author / Citation

If this tool is used in a scientific project, thesis, or publication, cite the corresponding GitHub repository and version or commit used for the analysis.

A suggested software citation format is:

Author. codon_harmonization_tool. Version/commit. GitHub repository, year.

Replace the placeholders with the final repository metadata.

Notes

codon_harmonization_tool is intended as a transparent and easily inspectable implementation of score-based codon harmonization. The algorithm deliberately remains simple so that every codon replacement can be traced back to the source and target codon frequency tables.
