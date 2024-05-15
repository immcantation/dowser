**formatClones** - *Generate an ordered list of airrClone objects for lineage construction*

Description
--------------------

`formatClones` takes a `data.frame` or `tibble` with AIRR or 
Change-O style columns as input and masks gap positions, masks ragged ends, 
removes duplicates sequences, and merges annotations associated with duplicate
sequences. If specified, it will un-merge duplicate sequences with different 
values specified in the `traits` option. It returns a list of `airrClone`
objects ordered by number of sequences which serve as input for lineage reconstruction.


Usage
--------------------
```
formatClones(
data,
seq = "sequence_alignment",
clone = "clone_id",
subgroup = "clone_subgroup",
id = "sequence_id",
germ = "germline_alignment_d_mask",
v_call = "v_call",
j_call = "j_call",
junc_len = "junction_length",
mask_char = "N",
max_mask = 0,
pad_end = TRUE,
text_fields = NULL,
num_fields = NULL,
seq_fields = NULL,
add_count = TRUE,
verbose = FALSE,
collapse = TRUE,
cell = "cell_id",
locus = "locus",
traits = NULL,
mod3 = TRUE,
randomize = TRUE,
use_regions = TRUE,
dup_singles = FALSE,
nproc = 1,
chain = "H",
heavy = "IGH",
filterstop = FALSE,
minseq = 2,
split_light = FALSE,
light_traits = FALSE,
majoronly = FALSE,
columns = NULL
)
```

Arguments
-------------------

data
:   data.frame containing the AIRR or Change-O data for a clone.
See [makeAirrClone](makeAirrClone.md) for required columns and their defaults

seq
:   name of the column containing observed DNA sequences. All 
sequences in this column must be multiple aligned.

clone
:   name of the column containing the identifier for the clone. All 
entries in this column should be identical.

subgroup
:   name of the column containing the identifier for the subgroup.

id
:   name of the column containing sequence identifiers.

germ
:   name of the column containing germline DNA sequences. All entries 
in this column should be identical for any given clone, and they
must be multiple aligned with the data in the `seq` column.

v_call
:   name of the column containing V-segment allele assignments. All 
entries in this column should be identical to the gene level.

j_call
:   name of the column containing J-segment allele assignments. All 
entries in this column should be identical to the gene level.

junc_len
:   name of the column containing the length of the junction as a 
numeric value. All entries in this column should be identical 
for any given clone.

mask_char
:   character to use for masking and padding.

max_mask
:   maximum number of characters to mask at the leading and trailing
sequence ends. If `NULL` then the upper masking bound will 
be automatically determined from the maximum number of observed 
leading or trailing Ns amongst all sequences. If set to `0` 
(default) then masking will not be performed.

pad_end
:   if `TRUE` pad the end of each sequence with `mask_char`
to make every sequence the same length.

text_fields
:   text annotation columns to retain and merge during duplicate removal.

num_fields
:   numeric annotation columns to retain and sum during duplicate removal.

seq_fields
:   sequence annotation columns to retain and collapse during duplicate 
removal. Note, this is distinct from the `seq` and `germ` 
arguments, which contain the primary sequence data for the clone
and should not be repeated in this argument.

add_count
:   if `TRUE` add an additional annotation column called 
`COLLAPSE_COUNT` during duplicate removal that indicates the 
number of sequences that were collapsed.

verbose
:   passed on to `collapseDuplicates`. If `TRUE`, report the 
numbers of input, discarded and output sequences; otherwise, process
sequences silently.

collapse
:   collapse identical sequences?

cell
:   name of the column containing cell assignment information

locus
:   name of the column containing locus information

traits
:   column ids to keep distinct during sequence collapse

mod3
:   pad sequences to length mutliple three?

randomize
:   randomize sequence order? Important if using PHYLIP

use_regions
:   assign CDR/FWR regions?

dup_singles
:   Duplicate sequences in singleton clones to include them as trees?

nproc
:   number of cores to parallelize formating over.

chain
:   if HL, include light chain information if available.

heavy
:   name of heavy chain locus (default = "IGH")

filterstop
:   only use sequences that do not contain an in-frame stop codon

minseq
:   minimum number of sequences per clone

split_light
:   split or lump subgroups? See `resolveLightChains`.

light_traits
:   Include the traits from the light chain when concatenating and collapsing trees?

majoronly
:   only return largest subgroup and sequences without light chains

columns
:   additional data columns to include in output




Value
-------------------

A tibble of [airrClone](airrClone-class.md) objects containing modified clones.


Details
-------------------

This function is a wrapper for [makeAirrClone](makeAirrClone.md). Also removes whitespace,
;, :, and = from ids



Examples
-------------------

```R
data(ExampleAirr)
# Select two clones, for demonstration purpose
sel <- c("3170", "3184")
clones <- formatClones(ExampleAirr[ExampleAirr$clone_id %in% sel,],traits="sample_id")
```



See also
-------------------

Executes in order [makeAirrClone](makeAirrClone.md). Returns a tibble of 
[airrClone](airrClone-class.md) objects 
which serve as input to [getTrees](getTrees.md) and [findSwitches](findSwitches.md).






