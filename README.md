
This is work in progress (barely) to provide a consistent interface to
different variant annotations such as from [snpEff ANN field](http://snpeff.sourceforge.net/) and the [VEP CSQ field](http://www.ensembl.org/info/docs/tools/vep/index.html).

This will be used in [gemini](http://gemini.rtfd.org/) but should also be of
general utility.

Design
======

We will have an effect base-class and then a sub-class for each tool (one
for snpeff and one for VEP).

`Effect` objects will be orderable (via \_\_le\_\_ ) and should have an equality test.  Then we can use [functools.total_ordering](https://docs.python.org/2/library/functools.html#functools.total_ordering) to provide the other comparison operators.

Given 2 effects, `a` and `b`: `a < b == True` iff the *severity* of `b` is greater than `a`.

We will have a classmethod: `Effect.top_severity([eff1, ... effn]) that will return the single highest serverity if that exists or
a list of the ties for highest
