**makeAirrClone** - *Generate a airrClone object for lineage construction*

Description
--------------------

`makeAirrClone` takes a data.frame with AIRR or Change-O style columns as input and 
masks gap positions, masks ragged ends, removes duplicates sequences, and merges 
annotations associated with duplicate sequences. It returns a `airrClone` 
object which serves as input for lineage reconstruction.


Usage
--------------------
```
makeAirrClone(
data,
id = "sequence_id",
seq = "sequence_alignment",
germ = "germline_alignment_d_mask",
v_call = "v_call",
j_call = "j_call",
junc_len = "junction_length",
clone = "clone_id",
subgroup = "clone_subgroup",
mask_char = "N",
max_mask = 0,
pad_end = TRUE,
text_fields = NULL,
num_fields = NULL,
seq_fields = NULL,
add_count = TRUE,
verbose = FALSE,
collapse = TRUE,
chain = "H",
heavy = NULL,
cell = "cell_id",
locus = "locus",
traits = NULL,
mod3 = TRUE,
randomize = TRUE,
use_regions = TRUE,
dup_singles = FALSE,
light_traits = FALSE
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

clone
:   name of the column containing the identifier for the clone. All 
entries in this column should be identical.

subgroup
:   name of the column containing the identifier for the subgroup.

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

chain
:   if HL, include light chain information if available.

heavy
:   name of heavy chain locus (default = "IGH")

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

light_traits
:   Include the traits from the light chain when concatenating and collapsing trees?




Value
-------------------

A [airrClone](airrClone-class.md) object containing the modified clone.


Details
-------------------

The input data.frame (`data`) must columns for each of the required column name 
arguments: `id`, `seq`, `germ`, `v_call`, `j_call`, 
`junc_len`, and `clone`.  
Additional annotation columns specified in the `traits`, `text_fields`, 
`num_fields` or `seq_fields` arguments will be retained in the `data` 
slot of the return object, but are not required. These options differ by their behavior
among collapsed sequences. Identical sequences that differ by any values specified in the
`traits` option will be kept distinct. Identical sequences that differ only by
values in the `num_fields` option will be collapsed and the values of their 
`num_fields` columns will be added together. Similar behavior occurs with 
`text_fields` but the unique values will concatenated with a comma.

The default columns are IMGT-gapped sequence columns, but this is not a requirement. 
However, all sequences (both observed and germline) must be multiple aligned using
some scheme for both proper duplicate removal and lineage reconstruction. 

The value for the germline sequence, V-segment gene call, J-segment gene call, 
junction length, and clone identifier are determined from the first entry in the 
`germ`, `v_call`, `j_call`, `junc_len` and `clone` columns, 
respectively. For any given clone, each value in these columns should be identical.

To allow for cases where heavy and light chains are used, this function returns three
sequence columns for heavy chains (sequence), light chain (lsequence, empty if none 
available), and concatenated heavy+light chain (hlsequence). These contain sequences
in alignment with germline, lgermline, and hlgermline slots, respectively. The sequence
column used for build trees is specified in the `phylo_seq` slot. Importantly, 
this column is also the sequence column that also has uninformative columns removed
by `cleanAlignment`. It is highly likely we will change this system to a single 
`sequence` and `germline` slot in the near future.

The airrClone object also contains vectors `locus`, `region`, and 
`numbers`, which contain the locus, IMGT region, and IMGT number for each position
in the sequence column specified in `phylo_seq`. If IMGT-gapped sequences are not 
supplied, this will likely result in an error. Specify `use_regions=FALSE` if not
using IMGT-gapped sequences



Examples
-------------------

```R
data(ExampleAirr)
airr_clone <- makeAirrClone(ExampleAirr[ExampleAirr$clone_id=="3184",])
```



See also
-------------------

Returns an [airrClone](airrClone-class.md). See [formatClones](formatClones.md) to generate an 
ordered list of airrClone objects.






