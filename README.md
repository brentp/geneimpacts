Given multiple snpEff or VEP or BCFTools consequence annotations for a single variant, get an orderable python object for each annotation.

[![Build Status](https://travis-ci.org/brentp/geneimpacts.svg?branch=master)](https://travis-ci.org/brentp/geneimpacts)

This is to provide a consistent interface to
different variant annotations such as from [snpEff ANN field](http://snpeff.sourceforge.net/) and the [VEP CSQ field](http://www.ensembl.org/info/docs/tools/vep/index.html).
and the [BCFTools consequence field](http://biorxiv.org/content/early/2016/12/01/090811)

This will be used in [gemini](http://gemini.rtfd.org/) but should also be of
general utility.

Design
======

There is an effect base-class and then a sub-class for `snpEff`, one for `VEP`, and one for `BCFT`

`Effect` objects are orderable (via \_\_le\_\_ ) and should have an \_\_eq\_\_ method so that we can use [functools.total_ordering](https://docs.python.org/2/library/functools.html#functools.total_ordering) to provide the other comparison operators.

Given 2 effects objects, `a` and `b`: `a < b == True` iff the *severity* of `b` is greater than `a`.

We will have a classmethod: `Effect.top_severity([eff1, ... effn]) that will return the single highest serverity if that exists or
a list of the ties for highest

Rules for severity:
===================

Given 2 annotations, *a* and *b*
*a* is more severe than *b* if:

1. *b* is a pseudogene and *a* is not
2. *a* is coding and *b* is not
3. *a* has higher severity than *b* ( see below)
4. polyphen, then sift
5. ??? transcript length? (we dont have access to this).

severity
--------

Severity is based on the [impacts from VEP](http://uswest.ensembl.org/info/docs/tools/vep/script/vep_other.html#pick)
and the [impacts from snpEff](http://snpeff.sourceforge.net/VCFannotationformat_v1.0.pdf). We reduce from the 4 categories HIGH, MEDIUM, LOW, MODIFIER to 3 by renaming MEDIUM to MED and renaming MODIFIER to LOW.
