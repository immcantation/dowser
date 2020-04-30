**formatClones** - *Generate an ordered list of airrClone objects for lineage construction*

Description
--------------------

`formatClones` takes a `data.frame` or `tibble` with AIRR or 
Change-O style columns as input and masks gap positions, masks ragged ends, 
removes duplicates sequences, and merges annotations associated with duplicate
sequences. If specified, it will un-merge duplicate sequences with different 
values specified in the `trait` option. It returns a list of `airrClone`
objects ordered by number of sequences which serve as input for lineage reconstruction.


Usage
--------------------
```
formatClones(
data,
id = "sequence_id",
seq = "sequence_alignment",
germ = "germline_alignment_d_mask",
vcall = "v_call",
jcall = "j_call",
junc_len = "junction_length",
clone = "clone_id",
subclone = "subclone_id",
mask_char = "N",
max_mask = 0,
pad_end = TRUE,
text_fields = NULL,
num_fields = NULL,
seq_fields = NULL,
add_count = TRUE,
verbose = FALSE,
nproc = 1,
collapse = TRUE,
region = "H",
heavy = NULL,
cell = "cell_id",
locus = "locus",
traits = NULL,
minseq = 2,
triple = TRUE,
subclones = "lump",
majoronly = FALSE
)
```

Arguments
-------------------

data
:   data.frame containing the AIRR or Change-O data for a clone. See Details
for the list of required columns and their default values.

id
:   name of the column containing sequence identifiers.

seq
:   name of the column containing observed DNA sequences. All 
sequences in this column must be multiple aligned.

germ
:   name of the column containing germline DNA sequences. All entries 
in this column should be identical for any given clone, and they
must be multiple aligned with the data in the `seq` column.

vcall
:   name of the column containing V-segment allele assignments. All 
entries in this column should be identical to the gene level.

jcall
:   name of the column containing J-segment allele assignments. All 
entries in this column should be identical to the gene level.

junc_len
:   name of the column containing the length of the junction as a 
numeric value. All entries in this column should be identical 
for any given clone.

clone
:   name of the column containing the identifier for the clone. All 
entries in this column should be identical.

subclone
:   name of the column containing the identifier for the subclone.

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

nproc
:   Number of cores to parallelize formating over.

collapse
:   iollapse identical sequences?

region
:   if HL, include light chain information if available.

heavy
:   name of heavy chain locus (default = "IGH")

cell
:   name of the column containing cell assignment information

locus
:   name of the column containing locus information

traits
:   column ids to keep distinct during sequence collapse

minseq
:   minimum numbner of sequences per clone

triple
:   pad sequences to length multiple of three

subclones
:   split or lump subclones? See `getSubclones`.

majoronly
:   only return largest subclone and sequences without light chains




Value
-------------------

A list of [airrClone](airrClone-class.md) objects containing modified clones.


Details
-------------------

This function is largely a wrapper for alakazam::makeAirrClone.
The input data.frame (`data`) must columns for each of the required column name 
arguments: `id`, `seq`, `germ`, `vcall`, `jcall`, 
`junc_len`, and `clone`.  The default values are as follows:

+ `id       = "sequence_id"`:           unique sequence identifier.
+ `seq      = "sequence_alignment"`:         IMGT-gapped sample sequence.
+ `germ     = "germline_alignment_d_mask"`:  IMGT-gapped germline sequence.
+ `vcall    = "v_call"`:                V-segment allele call.
+ `jcall    = "j_call"`:                J-segment allele call.
+ `junc_len = "junction_length"`:       junction sequence length.
+ `clone    = "clone_id"`:                 clone identifier.

Additional annotation columns specified in the `text_fields`, `num_fields` 
or `seq_fields` arguments will be retained in the `data` slot of the return 
object, but are not required. If the input data.frame `data` already contains a 
column named `sequence`, which is not used as the `seq` argument, then that 
column will not be retained.

The default columns are IMGT-gapped sequence columns, but this is not a requirement. 
However, all sequences (both observed and germline) must be multiple aligned using
some scheme for both proper duplicate removal and lineage reconstruction. 

The value for the germline sequence, V-segment gene call, J-segment gene call, 
junction length, and clone identifier are determined from the first entry in the 
`germ`, `vcall`, `jcall`, `junc_len` and `clone` columns, 
respectively. For any given clone, each value in these columns should be identical.



Examples
-------------------

```R
### Not run:
data(ExampleDb)
# clones <- formatClones(ExampleDb,trait="sample_id")
```



See also
-------------------

Executes in order `alakazam::maskSeqGaps`, `alakazam::maskSeqEnds`, 
`alakazam::padSeqEnds`, `alakazam::collapseDuplicates`,
			 and `processClones`. Returns a list of [airrClone](airrClone-class.md) objects 
			which serve as input to [getTrees](getTrees.md) and [bootstrapTrees](bootstrapTrees.md).






