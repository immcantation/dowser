**create_traitset** - *Takes an airr clone object and returns BEAST2 XML for a trait/traitSet from a column*

Description
--------------------

Takes an airr clone object and returns BEAST2 XML for a trait/traitSet from a column


Usage
--------------------
```
create_traitset(
clone,
trait_name,
column,
id,
trait_data_type = NULL,
isSet = FALSE,
include_germline_as_tip = FALSE
)
```

Arguments
-------------------

clone
:   an `airrClone` object

trait_name
:   name of the trait

column
:   column in the clone data to use for the trait

id
:   unique identifer for this analysis

trait_data_type
:   optional data type for the trait

isSet
:   is this a traitSet (TRUE) or a trait (FALSE)?

include_germline_as_tip
:   include the germline as a tip




Value
-------------------

String of XML of the trait or traitSet









