# {{title}}

This document is the output of the [VR-verify](https://github.com/EthanHolleman/VR-verify) workflow 
and was generated on {{current_date}}. Each section details Sanger sequencing results
produced with the intent to confirm the presence of variable region inserts
within a plasmid vector. For each variable region insert, the workflow will
suggest a status for a given sample, either confirmed or unconfirmed, based
on the content and quality of the Sanger sequencing data. This determination
is made off of three primary metrics which are described in the sections below.

## BlAST alignment of Sanger read to VR insert reference sequence

First, the workflow creates a local BLAST database which includes all
Sanger reads supplied to the workflow and then queries against it
with the reference sequences of the VR inserts. Sanger reads
are then filtered by their alignment characteristics, namely by
alignment length and percent identity of alignment. For the run
that produced this document minimum alignment length was {{min_align_length}}
bp and minim percent identity was set to {{min_pident}}. Reads that
did not mean these criteria were excluded from further analysis.

## Simulated EcoRI and SacI double digest

Insertion of VR inserts into the pFC9 relies on homology between
the insert and pFC9 sequences near it's EcoRI and SacI sites. Before
ligation, pFC9 is digested with EcoRI and SacI and the large fragment
is extracted and used in the ligation reaction. If insertion is successful
these restriction sites will be reconstituted and digestion of the
ligation product is expected to produce a 223 bp fragment. This
digestion is simulated for all Sanger reads and the difference in
length between the fragment produced and the expected 223 bp
fragment is measured.

## Insert distance from start of read



