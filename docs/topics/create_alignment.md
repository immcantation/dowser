**create_alignment** - *Takes an airr clone object and returns BEAST2 Alignment xml of the sequences*

Description
--------------------

Takes an airr clone object and returns BEAST2 Alignment xml of the sequences


Usage
--------------------
```
create_alignment(clone, id, include_germline_as_tip)
```

Arguments
-------------------

clone
:   an `airrClone` object

id
:   unique identifer for this analysis

include_germline_as_tip
:   include the germline as a tip in the alignment?




Value
-------------------

String of BEAST2 Alignment and TaxonSet xml









