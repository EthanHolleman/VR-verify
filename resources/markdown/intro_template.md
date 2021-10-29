# {title}

This document is the output of the [VR-verify](https://github.com/EthanHolleman/VR-verify) workflow 
and was generated on {current_date}. Each section details Sanger sequencing results
produced with the intent to confirm the presence of variable region inserts
within a plasmid vector. For each variable region insert, the workflow will
suggest a status for a given sample, either confirmed or unconfirmed, based
on the content and quality of the Sanger sequencing data. This determination
is made off of three primary metrics which are described in the sections below.

## BLAST alignment of Sanger read to VR insert reference sequence

First, the workflow creates a local BLAST database which includes all
Sanger reads supplied to the workflow and then queries against it
with the reference sequences of the VR inserts. Sanger reads
are then filtered by their alignment characteristics, namely by
alignment length and percent identity of alignment. For the run
that produced this document minimum alignment length was `{min_align_length}`
bp and minim percent identity was set to `{min_pident}`. Reads that
did not mean these criteria were excluded from further analysis.

## Simulated EcoRI and SacI double digest

Insertion of VR inserts into the pFC9 relies on homology between
the insert and pFC9 sequences near it's EcoRI and SacI sites. Before
ligation, pFC9 is digested with EcoRI and SacI and the large fragment
is extracted and used in the ligation reaction. If insertion is successful
these restriction sites will be reconstituted and digestion of the
ligation product is expected to produce a 223 bp fragment. This
digestion is simulated Sanger reads that used primers that do not directly
target the insert itself (most common). There are circumstances where
primers that directly target the insert were used and in these cases 
this simulated digest is not meaningful as the recognition site is usually
within the early region of the read making it un-resolvable. For
reads this metric is calculated for, the difference in
length between the fragment produced and the expected 223 bp
fragment is measured.

## Insert distance from start of read

Since the primer used in each sequencing reaction is known, we also know
how far the insert should be from the start of the read if the ligation
into the backbone plasmid was successful. Reads that vary significantly
from this expected value may have deletions or other abnormalities.

## Summary of Sanger sequencing results

The figure below shows summary statistics for all reads and all runs available
to the workflow at the time of execution ({current_date}).

![]({summary_plot_png})

The table below shows the status (confidence a read contains an insert) 
of the best read mapped to each insert. Inserts with no mapped reads
will have a value of `No mapped reads`.

{confidence_table}

The best confidence for each insert is summarized as a "completeness" table
below showing the number of unique inserts per their highest confidence read.

{completeness_table}

If only counting high confidence reads, then the series is currently
{percent_complete}% complete.


