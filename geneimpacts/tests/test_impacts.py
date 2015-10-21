import sys
import os
import gzip
from geneimpacts import SnpEff, VEP, Effect, OldSnpEff


HERE = os.path.dirname(__file__)

def test_snpeff():

    ann = SnpEff("C|splice_donor_variant&splice_region_variant&splice_region_variant&intron_variant|HIGH|DDX11L1|ENSG00000223972|transcript|ENST00000518655|transcribed_unprocessed_pseudogene|3/3|n.734+2_734+3delAG||||||")

    assert ann.gene == "DDX11L1"
    assert ann.transcript == "ENST00000518655"
    assert ann.biotype == "transcribed_unprocessed_pseudogene", ann.biotype
    assert ann.consequences == 'splice_donor_variant&splice_region_variant&splice_region_variant&intron_variant'.split('&')
    assert ann.severity == 3
    assert ann.impact_severity == "HIGH"
    assert ann.aa_change is None
    assert ann.exon == '3/3', ann.exon
    assert not ann.coding
    assert ann.is_pseudogene

def test_vep():

    ann = VEP('missense_variant|tTt/tGt|F/C|ENSG00000186092|OR4F5|ENST00000335137|1/1|possibly_damaging(0.568)|deleterious(0)|113/305|protein_coding')
    assert ann.gene == 'OR4F5'
    assert ann.transcript == 'ENST00000335137'
    assert ann.aa_change == "F/C", ann.aa_change
    assert ann.consequences == ['missense_variant']
    assert ann.coding
    assert ann.biotype == "protein_coding"
    assert ann.severity == 2
    assert ann.impact_severity == "MED", ann.impact_severity
    assert not ann.is_pseudogene
    assert ann.polyphen_value == 0.568, ann.polyphen
    assert ann.polyphen_class == "possibly_damaging", ann.polyphen
    assert ann.sift_value == 0.0, ann.sift
    assert ann.sift_class == "deleterious", ann.sift


def test_veps():

    f = os.path.join(HERE, "vep-csqs.txt.gz")
    with gzip.open(f) as veps:
        for csq in (VEP(l.strip()) for l in veps):
            assert csq.severity in (1, 2, 3)
            assert csq.is_pseudogene in (True, False)
            assert csq.coding in (True, False)
            assert isinstance(csq.polyphen_value, float) or csq.polyphen_value is None
            csq.gene
            assert isinstance(csq.sift_value, float) or csq.sift_value is None

def test_snpeffs():
    f = os.path.join(HERE, "snpeff-anns.txt.gz")
    with gzip.open(f) as anns:
        for csq in (SnpEff(l.strip()) for l in anns):
            assert csq.severity in (1, 2, 3)
            assert csq.is_pseudogene in (True, False)
            assert csq.coding in (True, False)
            assert csq.polyphen_value is None

EFFECTS = [VEP("upstream_gene_variant|||ENSG00000223972|DDX11L1|ENST00000456328|||||processed_transcript"),
           VEP("downstream_gene_variant|||ENSG00000227232|WASH7P|ENST00000488147|||||unprocessed_pseudogene"),
           VEP("non_coding_exon_variant&nc_transcript_variant|||ENSG00000223972|DDX11L1|ENST00000456328|2/3||||processed_transcript"),
           VEP("non_coding_exon_variant&nc_transcript_variant|||ENSG00000223972|DDX11L1|ENST00000456328|2/3||||processed_transcript"),
           VEP("splice_region_variant&non_coding_exon_variant&nc_transcript_variant|||ENSG00000223972|DDX11L1|ENST00000456328|2/3||||processed_transcript"),
           VEP("splice_region_variant&non_coding_exon_variant&nc_transcript_variant|||ENSG00000223972|DDX11L1|ENST00000456328|2/3||||processed_transcript"),
           VEP("splice_region_variant&non_coding_exon_variant&nc_transcript_variant|||ENSG00000223972|DDX11L1|ENST00000456328|2/3||||processed_transcript"),
           VEP("intron_variant&nc_transcript_variant|||ENSG00000223972|DDX11L1|ENST00000450305|||||transcribed_unprocessed_pseudogene"),
           VEP("intron_variant&nc_transcript_variant|||ENSG00000223972|DDX11L1|ENST00000450305|||||transcribed_unprocessed_pseudogene"),
           VEP('missense_variant|tTt/tGt|F/C|ENSG00000186092|OR4F5|ENST00000335137|1/1|possibly_damaging(0.568)|deleterious(0)|113/305|protein_coding'),
           VEP("non_coding_exon_variant&nc_transcript_variant&feature_elongation|||ENSG00000223972|DDX11L1|ENST00000456328|3/3||||processed_transcript"),
           ]


def test_order():

    effects = sorted(EFFECTS)
    assert effects[-1].impact_severity == "MED"
    assert effects[0].impact_severity == "LOW"

def test_o2():

    keys = [x.strip() for x in "Effect | Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_change| Amino_Acid_length | Gene_Name | Gene_BioType | Coding | Transcript | Exon  | ERRORS | WARNINGS".split("|")]

    effects = [OldSnpEff(v, keys) for v in "DOWNSTREAM(MODIFIER|||||RP5-902P8.10|processed_transcript|NON_CODING|ENST00000434139|),DOWNSTREAM(MODIFIER|||||RP5-902P8.10|processed_transcript|NON_CODING|ENST00000453732|),INTRON(MODIFIER||||138|SCNN1D|protein_coding|CODING|ENST00000470022|3),INTRON(MODIFIER||||638|SCNN1D|protein_coding|CODING|ENST00000338555|3),INTRON(MODIFIER||||638|SCNN1D|protein_coding|CODING|ENST00000400928|2),INTRON(MODIFIER||||669|SCNN1D|protein_coding|CODING|ENST00000379110|6),INTRON(MODIFIER||||704|SCNN1D|protein_coding|CODING|ENST00000325425|2),INTRON(MODIFIER||||802|SCNN1D|protein_coding|CODING|ENST00000379116|5),INTRON(MODIFIER|||||SCNN1D|nonsense_mediated_decay|CODING|ENST00000379101|5),INTRON(MODIFIER|||||SCNN1D|processed_transcript|CODING|ENST00000467651|3)".split(",")]

    effects = sorted(effects)
    assert effects[-1].gene == "SCNN1D"


    effects = sorted([OldSnpEff(v, keys) for v in "DOWNSTREAM(MODIFIER||||85|FAM138A|protein_coding|CODING|ENST00000417324|),DOWNSTREAM(MODIFIER|||||FAM138A|processed_transcript|CODING|ENST00000461467|),DOWNSTREAM(MODIFIER|||||MIR1302-10|miRNA|NON_CODING|ENST00000408384|),EXON(MODIFIER|||||MIR1302-10|antisense|NON_CODING|ENST00000469289|1),INTRON(MODIFIER|||||MIR1302-10|antisense|NON_CODING|ENST00000473358|1),UPSTREAM(MODIFIER|||||WASH7P|unprocessed_pseudogene|NON_CODING|ENST00000423562|),UPSTREAM(MODIFIER|||||WASH7P|unprocessed_pseudogene|NON_CODING|ENST00000430492|),UPSTREAM(MODIFIER|||||WASH7P|unprocessed_pseudogene|NON_CODING|ENST00000438504|),UPSTREAM(MODIFIER|||||WASH7P|unprocessed_pseudogene|NON_CODING|ENST00000488147|),UPSTREAM(MODIFIER|||||WASH7P|unprocessed_pseudogene|NON_CODING|ENST00000538476|)".split(",")])
    s = "\n".join(e.effect_string for e in effects[::-1])

    # reversed so that most significant is first
    assert s == """\
DOWNSTREAM(MODIFIER||||85|FAM138A|protein_coding|CODING|ENST00000417324|
DOWNSTREAM(MODIFIER|||||FAM138A|processed_transcript|CODING|ENST00000461467|
INTRON(MODIFIER|||||MIR1302-10|antisense|NON_CODING|ENST00000473358|1
EXON(MODIFIER|||||MIR1302-10|antisense|NON_CODING|ENST00000469289|1
DOWNSTREAM(MODIFIER|||||MIR1302-10|miRNA|NON_CODING|ENST00000408384|
UPSTREAM(MODIFIER|||||WASH7P|unprocessed_pseudogene|NON_CODING|ENST00000423562|
UPSTREAM(MODIFIER|||||WASH7P|unprocessed_pseudogene|NON_CODING|ENST00000430492|
UPSTREAM(MODIFIER|||||WASH7P|unprocessed_pseudogene|NON_CODING|ENST00000438504|
UPSTREAM(MODIFIER|||||WASH7P|unprocessed_pseudogene|NON_CODING|ENST00000488147|
UPSTREAM(MODIFIER|||||WASH7P|unprocessed_pseudogene|NON_CODING|ENST00000538476|"""

def test_highest():
    effects = sorted(EFFECTS)

    top = Effect.top_severity(effects)
    assert isinstance(top, Effect)
    assert top.impact_severity == "MED"


    effects.append(effects[-1])

    top = Effect.top_severity(effects)
    assert isinstance(top, list)
    assert top[0].impact_severity == "MED"

def test_splice():

    e = VEP('splice_acceptor_variant&intron_variant&feature_truncation|||ENSG00000221978|CCNL2|ENST00000408918||||-/226|protein_coding|1')
    assert (e.is_coding, e.is_exonic, e.is_splicing) == (False, False, True)

    e = VEP('intron_variant&feature_elongation|||ENSG00000187634|SAMD11|ENST00000341065||||-/589|protein_coding|1')
    assert (e.is_coding, e.is_exonic, e.is_splicing) == (False, False, False)

def test_eff_splice():

    keys = [x.strip() for x in "Effect | Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_change| Amino_Acid_length | Gene_Name | Gene_BioType |  Coding | Transcript | Exon  | ERRORS | WARNINGS".split("|")]
    e = OldSnpEff("SPLICE_SITE_REGION+SYNONYMOUS_CODING(LOW|SILENT|acG/acA|T245|1134|ANKS1A|protein_coding|CODING|ENST00000360359|5|A)", keys)
    print e.effects
    assert e.aa_change == "T245"
    assert e.is_coding, e.is_coding
    # note that we choose synonymous coding over splice_site_region
    assert not e.is_splicing
