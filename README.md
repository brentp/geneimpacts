
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

Rules for severity:
===================

Given 2 annotations, *a* and *b*
*a* is more severe than *b* if:

1. *b* is a pseudogene and *a* is not
2. *a* is coding and *b* is not
3. *a* has higher severity than *b* ( see below)
4. polyphen/sift?


severity
--------

Here are the impacts from both snpEff and VEP (as of snpEff 4+, they both use Sequence ontology terms)
from 20 million variants ordered by frequency:

We need to decide a severity ordering for these.

```
85964122 intron_variant
13478317 upstream_gene_variant
13429780 downstream_gene_variant
11645981 nc_transcript_variant
9610719 intergenic_region
7745274 intergenic_variant
7356090 sequence_feature
4529431 feature_truncation
4428300 feature_elongation
3588296 NMD_transcript_variant
1147256 non_coding_exon_variant
 765406 3_prime_UTR_variant
 407180 missense_variant
 339328 synonymous_variant
 236790 5_prime_UTR_variant
 180396 splice_region_variant
  21729 TF_binding_site_variant
  17067 5_prime_UTR_premature_start_codon_gain_variant
  16092 frameshift_variant
   8689 splice_acceptor_variant
   8503 splice_donor_variant
   7312 stop_gained
   6245 inframe_deletion
   5265 inframe_insertion
   3832 intragenic_variant
   3154 disruptive_inframe_deletion
   2017 disruptive_inframe_insertion
   1010 mature_miRNA_variant
    734 coding_sequence_variant
    679 stop_lost
    626 initiator_codon_variant
    479 start_lost
    465 transcript
    323 stop_retained_variant
     50 incomplete_terminal_codon_variant
     16 exon_loss_variant
     12 5_prime_UTR_truncation
      5 transcript_ablation
      2 non_canonical_start_codon
      1 3_prime_UTR_truncation
      1 exon_loss
```
