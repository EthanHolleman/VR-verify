## {section_title}

Sanger sequencing results for reads aligning to {insert_name}. A total
of {number_of_reads} aligned to the {insert_name} reference sequence,
a summary of these reads is shown in the table below.

{summary_table}

In addition to BLAST alignment reads were assessed by a simulated
EcoRI and SacI digest, the distance of the insert from the start of
the read given the specific primer used for the read and
read quality determined by the mean Phred score of the read's trace.
Results are shown in the figure below. If no figure is visible then
no reads aligned successfully. 

![]({summary_plot_path})

Read read was automatically accessed for basic characteristics that
indicate the likiehood of a successful ligation reaction based
on quality of the read, simulated restriction digest,
location of the insert in the read given the specific primers used,
and most heavily weighted quality of BLAST alignment to the insert
reference sequence.

{confidence_table}

