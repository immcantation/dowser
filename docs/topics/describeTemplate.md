**describeTemplate** - *Describe an xml-templates file: what it needs and what it supports*

Description
--------------------

Reports every free-form `${VAR}` placeholder in a template, with
descriptions pulled from an embedded `dowser-template-info` comment
block if the template author provided one, plus a built-in checklist of
structural elements Dowser recognizes (MRCA, ROOTFREQS, OPERATORS,
ROOT_TYPE_PROBABILITIES, a `nodeTypes` parameter, and whether
the template uses `GermlineRootTree`). The structural checklist
works on any template, documented or not, since that behavior is defined
by Dowser's own code rather than by the template author -- several of
these are silent no-ops rather than errors when an argument is supplied
but the template doesn't support it, which is exactly what this function
is meant to catch before a long MCMC run rather than after. If the
template's `dowser-template-info` block has an overall
`@description` entry, that's printed first, before the structural
checklist and variable list.


Usage
--------------------
```
describeTemplate(template, print = TRUE, var = NULL)
```

Arguments
-------------------

template
:   path to an xml-templates file

print
:   print a human-readable summary? Default TRUE. Set
FALSE if you only want the returned object.

var
:   optional: a single variable name (bare, no
`${}`), free-form or Dowser-special (e.g.
`"EMP_EQFREQS"`). If supplied, instead of the
full dump, this looks up just that one variable and
prints/returns a focused answer: whether it's in
the template at all, and either its
`dowser-template-info` description (free-form)
or its built-in structural effect (Dowser-special).




Value
-------------------

invisibly: if `var` is NULL, a list with
`$description` (character scalar, or NULL if the template
has no `@description`), `$vars` (data.frame: name,
tier, documented, description), and `$structural`
(data.frame: element, present, effect). If `var` is
supplied, a list with `$name`, `$kind` ("free",
"structural", or "not_found"), `$present`, and whichever
of `$documented`/`$tier`/`$description` (free)
or `$effect` (structural) applies.


Details
-------------------

See the xml-templates variant reference at https://tyche.readthedocs.io
for the `dowser-template-info` format if you're documenting your own
template.



Examples
-------------------

```R
### Not run:
describeTemplate("TypeLinked/TypeLinkedExpectedOccupancy_EstTraitClockRates_EmpFreq.xml")

```

**Error in describeTemplate("TypeLinked/TypeLinkedExpectedOccupancy_EstTraitClockRates_EmpFreq.xml")**: could not find function "describeTemplate"
```R
# describeTemplate("TypeLinked/TypeLinkedExpectedOccupancy_EstTraitClockRates_EmpFreq.xml",
# var = "TRAIT_RATE_MEAN_1")

```








