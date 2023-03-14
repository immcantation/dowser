**airrClone-class** - *S4 class defining a clone in Dowser*

Description
--------------------

`airrClone` defines a common data structure for perform lineage recontruction
from AIRR data, based heavily on alakazam::ChangeoClone.






Slots
-------------------



`data`
:   data.frame containing sequences and annotations. Contains the
columns `sequence_id` and `sequence`, as well as any additional 
sequence-specific annotation columns

`clone`
:   string defining the clone identifier

`germline`
:   string containing the heavy chain germline sequence for the clone

`lgermline`
:   string containing the light chain germline sequence for the clone

`hlgermline`
:   string containing the combined germline sequence for the clone

`v_gene`
:   string defining the V segment gene call

`j_gene`
:   string defining the J segment gene call

`junc_len`
:   numeric junction length (nucleotide count)

`locus`
:   index showing which locus represented at each site

`region`
:   index showing FWR/CDR region for each site

`phylo_seq`
:   sequence column used for phylogenetic tree building

`numbers`
:   index (usually IMGT) number of each site in `phylo_seq`




See also
-------------------

See [formatClones](formatClones.md) for use.






