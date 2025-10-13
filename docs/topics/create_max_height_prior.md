**create_max_height_prior** - *Takes an airr clone object and returns BEAST2 XML to set a maximum height prior*

Description
--------------------

Takes an airr clone object and returns BEAST2 XML to set a maximum height prior


Usage
--------------------
```
create_max_height_prior(clone, id, max_start_date)
```

Arguments
-------------------

clone
:   an `airrClone` object

id
:   unique identifer for this analysis

max_start_date
:   max start date to use for prior, in forward time




Value
-------------------

String of XML setting the MRCA prior of the observed sequences









