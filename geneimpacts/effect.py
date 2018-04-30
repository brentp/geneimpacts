from __future__ import print_function
import sys

from functools import total_ordering
import re
import itertools as it
try:
    izip = it.izip
except AttributeError:
    izip = zip
    basestring = str

old_snpeff_effect_so = {'CDS': 'coding_sequence_variant',
               'CODON_CHANGE': 'coding_sequence_variant',
               'CODON_CHANGE_PLUS_CODON_DELETION': 'disruptive_inframe_deletion',
               'CODON_CHANGE_PLUS_CODON_INSERTION': 'disruptive_inframe_insertion',
               'CODON_DELETION': 'inframe_deletion',
               'CODON_INSERTION': 'inframe_insertion',
               'DOWNSTREAM': 'downstream_gene_variant',
               'EXON': 'exon_variant',
               'EXON_DELETED': 'exon_loss_variant',
               'FRAME_SHIFT': 'frameshift_variant',
               'GENE': 'gene_variant',
               'INTERGENIC': 'intergenic_variant',
               'INTERGENIC_REGION': 'intergenic_region',
               'INTERGENIC_CONSERVED': 'conserved_intergenic_variant',
               'INTRAGENIC': 'intragenic_variant',
               'INTRON': 'intron_variant',
               'INTRON_CONSERVED': 'conserved_intron_variant',
               'NON_SYNONYMOUS_CODING': 'missense_variant',
               'RARE_AMINO_ACID': 'rare_amino_acid_variant',
               'SPLICE_SITE_ACCEPTOR': 'splice_acceptor_variant',
               'SPLICE_SITE_DONOR': 'splice_donor_variant',
               'SPLICE_SITE_REGION': 'splice_region_variant',
               #'START_GAINED': '5_prime_UTR_premature_start_codon_gain_variant',
               'START_GAINED': '5_prime_UTR_premature_start_codon_variant',
               'START_LOST': 'start_lost',
               'STOP_GAINED': 'stop_gained',
               'STOP_LOST': 'stop_lost',
               'SYNONYMOUS_CODING': 'synonymous_variant',
               'SYNONYMOUS_START': 'start_retained_variant',
               'SYNONYMOUS_STOP': 'stop_retained_variant',
               'TRANSCRIPT': 'transcript_variant',
               'UPSTREAM': 'upstream_gene_variant',
               'UTR_3_DELETED': '3_prime_UTR_truncation_+_exon_loss_variant',
               'UTR_3_PRIME': '3_prime_UTR_variant',
               'UTR_5_DELETED': '5_prime_UTR_truncation_+_exon_loss_variant',
               'UTR_5_PRIME': '5_prime_UTR_variant',
               'NON_SYNONYMOUS_START': 'initiator_codon_variant',
               'NONE': 'None',
               'CHROMOSOME_LARGE_DELETION': 'chromosomal_deletion'}

old_snpeff_lookup = {'CDS': 'LOW',
     'CHROMOSOME_LARGE_DELETION': 'HIGH',
     'CODON_CHANGE': 'MED',
     'CODON_CHANGE_PLUS_CODON_DELETION': 'MED',
     'CODON_CHANGE_PLUS_CODON_INSERTION': 'MED',
     'CODON_DELETION': 'MED',
     'CODON_INSERTION': 'MED',
     'DOWNSTREAM': 'LOW',
     'EXON': 'LOW',
     'EXON_DELETED': 'HIGH',
     'FRAME_SHIFT': 'HIGH',
     'GENE': 'LOW',
     'INTERGENIC': 'LOW',
     'INTERGENIC_CONSERVED': 'LOW',
     'INTRAGENIC': 'LOW',
     'INTRON': 'LOW',
     'INTRON_CONSERVED': 'LOW',
     'NONE': 'LOW',
     'NON_SYNONYMOUS_CODING': 'MED',
     'NON_SYNONYMOUS_START': 'HIGH',
     'RARE_AMINO_ACID': 'HIGH',
     'SPLICE_SITE_ACCEPTOR': 'HIGH',
     'SPLICE_SITE_DONOR': 'HIGH',
     'SPLICE_SITE_REGION': 'MED',
     'START_GAINED': 'LOW',
     'START_LOST': 'HIGH',
     'STOP_GAINED': 'HIGH',
     'STOP_LOST': 'HIGH',
     'SYNONYMOUS_CODING': 'LOW',
     'SYNONYMOUS_START': 'LOW',
     'SYNONYMOUS_STOP': 'LOW',
     'TRANSCRIPT': 'LOW',
     'UPSTREAM': 'LOW',
     'UTR_3_DELETED': 'MED',
     'UTR_3_PRIME': 'LOW',
     'UTR_5_DELETED': 'MED',
     'UTR_5_PRIME': 'LOW'}
# http://uswest.ensembl.org/info/genome/variation/predicted_data.html#consequences
IMPACT_SEVERITY = [
    ('chromosome_number_variation', 'HIGH'), # snpEff
    ('transcript_ablation', 'HIGH'), # VEP
    ('exon_loss_variant', 'HIGH'), # snpEff
    ('exon_loss', 'HIGH'), # snpEff
    ('rare_amino_acid_variant', 'HIGH'),
    ('protein_protein_contact', 'HIGH'), # snpEff
    ('structural_interaction_variant', 'HIGH'), #snpEff
    ('feature_fusion', 'HIGH'), #snpEff
    ('bidirectional_gene_fusion', 'HIGH'), #snpEff
    ('gene_fusion', 'HIGH'), #snpEff
    ('feature_ablation', 'HIGH'), #snpEff, structural varint
    ('splice_acceptor_variant', 'HIGH'), # VEP
    ('splice_donor_variant', 'HIGH'), # VEP
    ('start_retained_variant', 'HIGH'), # new VEP
    ('stop_gained', 'HIGH'), # VEP
    ('frameshift_variant', 'HIGH'), # VEP
    ('stop_lost', 'HIGH'), # VEP
    ('start_lost', 'HIGH'), # VEP
    ('transcript_amplification', 'HIGH'), # VEP


    ('disruptive_inframe_deletion', 'MED'), #snpEff
    ('conservative_inframe_deletion', 'MED'), #snpEff
    ('disruptive_inframe_insertion', 'MED'), #snpEff
    ('conservative_inframe_insertion', 'MED'), #snpEff
    ('duplication', 'MED'), # snpEff, structural variant
    ('inversion', 'MED'), # snpEff, structural variant
    ('exon_region', 'MED'), # snpEff, structural variant
    ('inframe_insertion', 'MED'), # VEP
    ('inframe_deletion', 'MED'), # VEP
    ('missense_variant', 'MED'), # VEP
    ('protein_altering_variant', 'MED'), # VEP
    ('initiator_codon_variant', 'MED'), # snpEff
    ('regulatory_region_ablation', 'MED'), # VEP

    ('5_prime_UTR_truncation', 'MED'), # found in snpEff
    ('splice_region_variant', 'MED'), # VEP changed to have medium priority

    ('3_prime_UTR_truncation', 'LOW'), # found in snpEff
    ('non_canonical_start_codon', 'LOW'), # found in snpEff

    ('synonymous_variant', 'LOW'), # VEP
    ('coding_sequence_variant', 'LOW'), # VEP
    ('incomplete_terminal_codon_variant', 'LOW'), # VEP
    ('stop_retained_variant', 'LOW'), # VEP
    ('mature_miRNA_variant', 'LOW'), # VEP
    ('5_prime_UTR_premature_start_codon_variant', 'LOW'), # snpEff
    ('5_prime_UTR_premature_start_codon_gain_variant', 'LOW'), #snpEff
    ('5_prime_UTR_variant', 'LOW'), # VEP
    ('3_prime_UTR_variant', 'LOW'), # VEP


    ('non_coding_transcript_exon_variant', 'LOW'), # VEP
    ('conserved_intron_variant', 'LOW'), # snpEff
    ('intron_variant', 'LOW'), # VEP
    ('exon_variant', 'LOW'), # snpEff
    ('gene_variant', 'LOW'), # snpEff
    ('NMD_transcript_variant', 'LOW'), # VEP
    ('non_coding_transcript_variant', 'LOW'), # VEP
    ('upstream_gene_variant', 'LOW'), # VEP
    ('downstream_gene_variant', 'LOW'), # VEP
    ('TFBS_ablation', 'LOW'), # VEP
    ('TFBS_amplification', 'LOW'), # VEP
    ('TF_binding_site_variant', 'LOW'), # VEP
    ('regulatory_region_amplification', 'LOW'), # VEP
    ('feature_elongation', 'LOW'), # VEP
    ('miRNA', 'LOW'), # snpEff
    ('transcript_variant', 'LOW'), # snpEff
    ('start_retained', 'LOW'), # snpEff
    ('regulatory_region_variant', 'LOW'), # VEP
    ('feature_truncation', 'LOW'), # VEP
    ('non_coding_exon_variant', 'LOW'),
    ('nc_transcript_variant', 'LOW'),
    ('conserved_intergenic_variant', 'LOW'), # snpEff
    ('intergenic_variant', 'LOW'), # VEP
    ('intergenic_region', 'LOW'), # snpEff
    ('intragenic_variant', 'LOW'), # snpEff
    ('non_coding_transcript_exon_variant', 'LOW'), # snpEff
    ('non_coding_transcript_variant', 'LOW'), # snpEff
    ('transcript', 'LOW'),  # ? snpEff older
    ('sequence_feature', 'LOW'), # snpEff older
    ('non_coding', 'LOW'), # BCSQ


    ('?', 'UNKNOWN'),  # some VEP annotations have '?'
    ('', 'UNKNOWN'),  # some VEP annotations have ''
    ('UNKNOWN', 'UNKNOWN'),  # some snpEFF annotations have 'unknown'
]

# bcftools doesn't add _variant on the end.
for (csq, imp) in list(IMPACT_SEVERITY[::-1]):
    if csq.endswith('_variant'):
        for i, (a, b) in enumerate(IMPACT_SEVERITY):
            if (a, b) == (csq, imp):
                IMPACT_SEVERITY.insert(i, (csq[:-8].lower(), imp))
                break

IMPACT_SEVERITY_ORDER = dict((x[0], i) for i, x in enumerate(IMPACT_SEVERITY[::-1]))
IMPACT_SEVERITY = dict(IMPACT_SEVERITY)

EXONIC_IMPACTS = set(["stop_gained",
                      "exon_variant",
                      "stop_lost",
                      "frameshift_variant",
                      "initiator_codon_variant",
                      "inframe_deletion",
                      "inframe_insertion",
                      "missense_variant",
                      "protein_altering_variant",
                      "incomplete_terminal_codon_variant",
                      "stop_retained_variant",
                      "5_prime_UTR_premature_start_codon_variant",
                      "synonymous_variant",
                      "coding_sequence_variant",
                      "5_prime_UTR_variant",
                      "3_prime_UTR_variant",
                      "transcript_ablation",
                      "transcript_amplification",
                      "feature_elongation",
                      "feature_truncation"])

for im in list(EXONIC_IMPACTS):
    if im.endswith("_variant"):
        EXONIC_IMPACTS.add(im[:-8])
EXONIC_IMPACTS = frozenset(EXONIC_IMPACTS)

def snpeff_aa_length(self):
    try:
        v = self.effects['AA.pos / AA.length']
        if v.strip():
            return int(v.split("/")[1].strip())
    except:
        try:
            return int(self.effects['Amino_Acid_length'])
        except:
            return None

def vep_aa_length(self):
    if not 'Protein_position' in self.effects:
        return None
    try:
        return int(self.effects['Protein_position'])
    except ValueError:
        try:
            return self.effects['Protein_position']
        except KeyError:
            return None

def vep_polyphen_pred(self):
    try:
        return self.effects['PolyPhen'].split('(')[0]
    except (KeyError, IndexError):
        return None

def vep_polyphen_score(self):
    try:
        return float(self.effects['PolyPhen'].split('(')[1][:-1])
    except (KeyError, IndexError):
        return None

def vep_sift_score(self):
    try:
        return float(self.effects['SIFT'].split("(")[1][:-1])
    except (IndexError, KeyError):
        return None

def vep_sift_pred(self):
    try:
        return self.effects['SIFT'].split("(")[0]
    except (IndexError, KeyError):
        return None

snpeff_lookup = {
    'transcript': ['Feature_ID', 'Transcript_ID', 'Transcript'],
    'gene': 'Gene_Name',
    'exon': ['Rank', 'Exon', 'Exon_Rank'],
    'codon_change': ['HGVS.c', 'Codon_Change'],
    'aa_change': ['HGVS.p', 'Amino_Acid_Change', 'Amino_Acid_change'],
    'aa_length': snpeff_aa_length,
    'biotype': ['Transcript_BioType', 'Gene_BioType'],
    'alt': 'Allele',
        }

bcft_lookup = {}

vep_lookup = {
    'transcript': 'Feature',
    'gene': ['SYMBOL', 'HGNC', 'Gene'],
    'exon': 'EXON',
    'codon_change': 'Codons',
    'aa_change': 'Amino_acids',
    'aa_length': vep_aa_length,
    'biotype': 'BIOTYPE',
    'polyphen_pred': vep_polyphen_pred,
    'polyphen_score': vep_polyphen_score,
    'sift_pred': vep_sift_pred,
    'sift_score': vep_sift_score,
    'alt': 'ALLELE',
        }

# lookup here instead of returning ''.
defaults = {'gene': None}

@total_ordering
class Effect(object):
    _top_consequence = None
    lookup = None

    def __init__(self, key, effect_dict, keys, prioritize_canonical):
        raise NotImplemented

    @classmethod
    def new(self, key, effect_dict, keys):
        lookup = {"CSQ": VEP, "ANN": SnpEff, "EFF": OldSnpEff, "BCSQ": BCFT}
        assert key in lookup
        return lookup[key](effect_dict, keys)

    @property
    def is_exonic(self):
        return self.top_consequence in EXONIC_IMPACTS

    def unused(self):
        return []

    @property
    def top_consequence(self):
        # sort by order and return the top
        if self._top_consequence is None:
            self._top_consequence = sorted([(IMPACT_SEVERITY_ORDER.get(c, 0), c) for c in
            self.consequences], reverse=True)[0][1]
        return self._top_consequence

    @property
    def so(self):
        return self.top_consequence

    @property
    def is_coding(self):
        return self.biotype == "protein_coding" and self.is_exonic and ("_UTR_" not in self.top_consequence)

    @property
    def is_splicing(self):
        return "splice" in self.top_consequence

    @property
    def is_lof(self):
        return self.biotype == "protein_coding" and self.impact_severity == "HIGH"

    def __le__(self, other):
        # we sort so that the effects with the highest impacts come last
        # (highest) and so, we:
        # + return true if self has lower impact than other.
        # + return false if self has higher impact than other.
        self_has_lower_impact = True
        self_has_higher_impact = False

        if self.prioritize_canonical:
            scanon, ocanon = self.is_canonical, other.is_canonical
            if scanon and not ocanon:
                return self_has_higher_impact
            elif ocanon and not scanon:
                return self_has_lower_impact

        spg = self.is_pseudogene
        opg = other.is_pseudogene
        if spg and not opg:
            return self_has_lower_impact
        elif opg and not spg:
            return self_has_higher_impact

        sc, oc = self.coding, other.coding
        if sc and not oc:
            # other is not coding. is is splicing?
            # if other is splicing, we have lower impact.
            if not (self.is_splicing or other.is_splicing):
                return self_has_higher_impact
        elif oc and not sc:
            # self. is not coding. is it splicing?
            # if self is splicing it has higher impact
            if not (self.is_splicing or other.is_splicing):
                return self_has_lower_impact

        if self.severity != other.severity:
            return self.severity <= other.severity

        if self.biotype == "protein_coding" and not other.biotype == "protein_coding":
            return False
        elif other.biotype == "protein_coding" and not self.biotype == "protein_coding":
            return True

        if self.biotype == "processed_transcript" and not other.biotype == "processed_transcript":
            return False
        elif other.biotype == "processed_transcript" and not self.biotype == "processed_transcript":
            return True

        # sift higher == more damaing
        if (self.sift_value or 10000) < (other.sift_value or 10000):
            return True

        # polyphen, lower == more damaging
        if (self.polyphen_value or -10000) > (other.polyphen_value or -10000):
            return True

        return max(IMPACT_SEVERITY_ORDER.get(c, 0) for c in self.consequences) <= \
                        max(IMPACT_SEVERITY_ORDER.get(co, 0) for co in other.consequences)

    @classmethod
    def top_severity(cls, effects):
        for i, e in enumerate(effects):
            if isinstance(e, basestring):

                effects[i] = cls(e)

        if len(effects) == 0:
            return None
        if len(effects) == 1:
            return effects[0]
        effects = sorted(effects)
        if effects[-1] > effects[-2]:
            return effects[-1]
        ret = [effects[-1], effects[-2]]
        for i in range(-3, -(len(effects) - 1), -1):
            if effects[-1] > effects[i]: break
            ret.append(effects[i])
        return ret

    def __getitem__(self, key):
        return self.effects[key]

    def __eq__(self, other):
        if not isinstance(other, Effect): return False
        return self.effect_string == other.effect_string

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return "%s(%s-%s, %s)" % (self.__class__.__name__, self.gene,
                self.consequence, self.impact_severity)

    @property
    def effect_severity(self):
        return self.impact_severity

    @property
    def lof(self):
        return self.biotype == "protein_coding" and self.impact_severity == "HIGH"

    @property
    def severity(self, lookup={'HIGH': 3, 'MED': 2, 'LOW': 1, 'UNKNOWN': 0}, sev=IMPACT_SEVERITY):
        # higher is more severe. used for ordering.
        try:
            v = max(lookup[sev[csq]] for csq in self.consequences)
        except KeyError:
            v = 0
        if v == 0:
            excl = []
            for i, c in [(i, c) for i, c in enumerate(self.consequences) if not c in sev]:
                sys.stderr.write("WARNING: unknown severity for '%s' with effect '%s'\n" % (self.effect_string, c))
                sys.stderr.write("Please report this on github with the effect-string above\n")
                excl.append(i)
            if len(excl) == len(self.consequences):
                v = 1
            else:
                v =  max(lookup[sev[csq]] for i, csq in enumerate(self.consequences) if not i in excl)
        return max(v, 1)

    @property
    def impact_severity(self):
        return ['xxx', 'LOW', 'MED', 'HIGH'][self.severity]

    @property
    def consequence(self):
        return self.top_consequence

    @property
    def is_pseudogene(self): #bool
        return self.biotype is not None and 'pseudogene' in self.biotype


    def __getattr__(self, k):
        v = self.lookup.get(k)
        if v is None: return v
        if isinstance(v, basestring):
            ret = self.effects.get(v)
            # if we didnt get value, there may be a column
            # specific value stored in defaults so we look import
            # up.
            if not ret and ret is not False:
                return defaults.get(k, '')
            return ret
        elif isinstance(v, list):
            for key in v:
                try:
                    return self.effects[key]
                except KeyError:
                    continue
            return defaults.get(k, '')
        return v(self)

class BCFT(Effect):
    __slots__ = ('effect_string', 'effects', 'biotype', 'gene', 'transcript', 'aa_change', 'dna_change')
    keys = "consequence,gene,transcript,biotype,strand,amino_acid_change,dna_change".split(",")
    lookup = bcft_lookup

    def __init__(self, effect_string, keys=None, prioritize_canonical=False):
        if keys is not None: self.keys = keys
        self.effect_string = effect_string
        self.effects = dict(izip(self.keys, (x.strip().replace(' ', '_') for x in effect_string.split("|"))))
        self.biotype = self.effects.get('biotype', None)
        self.transcript = self.effects.get('transcript', None)
        self.gene = self.effects.get('gene', None)
        self.aa_change = self.effects.get('amino_acid_change', None)
        self.consequences = self.effects[self.keys[0]].split('&')

    def unused(self, used=frozenset("csq|gene|transcript|biotype|strand|aa_change|dna_change".lower().split("|"))):
        """Return fields that were in the VCF but weren't utilized as part of the standard fields supported here."""
        return [k for k in self.keys if not k.lower() in used]

    @property
    def exonic(self):
        return self.biotype == "protein_coding" and any(csq in EXONIC_IMPACTS for csq in self.consequences)

    @property
    def coding(self):
        # what about start/stop_gained?
        return self.exonic and any(csq[1:] != "_prime_utr" for csq in self.consequences)


class VEP(Effect):
    __slots__ = ('effect_string', 'effects', 'biotype')
    keys = "Consequence|Codons|Amino_acids|Gene|SYMBOL|Feature|EXON|PolyPhen|SIFT|Protein_position|BIOTYPE|CANONICAL".split("|")
    lookup = vep_lookup

    def __init__(self, effect_string, keys=None, checks=True, prioritize_canonical=False):
        if checks:
            assert not "," in effect_string
            assert not "=" in effect_string
        self.effect_string = effect_string
        if keys is not None: self.keys = keys

        self.effect_string = effect_string
        self.effects = dict(izip(self.keys, (x.strip() for x in effect_string.split("|"))))
        self.biotype = self.effects.get('BIOTYPE', None)
        self.prioritize_canonical = prioritize_canonical

    @property
    def consequences(self, _cache={}):
        try:
            # this is a bottleneck so we keep a cache
            return _cache[self.effects['Consequence']]
        except KeyError:
            res = _cache[self.effects['Consequence']] = list(it.chain.from_iterable(x.split("+") for x in self.effects['Consequence'].split('&')))
            return res

    def unused(self, used=frozenset("Consequence|Codons|Amino_acids|Gene|SYMBOL|Feature|EXON|PolyPhen|SIFT|Protein_position|BIOTYPE|CANONICAL".lower().split("|"))):
        """Return fields that were in the VCF but weren't utilized as part of the standard fields supported here."""
        return [k for k in self.keys if not k.lower() in used]

    @property
    def coding(self):
        # what about start/stop_gained?
        return self.exonic and any(csq[1:] != "_prime_UTR_variant" for csq in self.consequences)

    @property
    def exonic(self):
        return self.biotype == "protein_coding" and any(csq in EXONIC_IMPACTS for csq in self.consequences)

    @property
    def is_canonical(self):
        return self.effects.get("CANONICAL", False)

class SnpEff(Effect):
    lookup = snpeff_lookup

    __slots__ = ('effects', 'effect_string', 'biotype')

    keys = [x.strip() for x in 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'.split("|")]

    def __init__(self, effect_string, keys=None, prioritize_canonical=False):
        assert not "," in effect_string
        assert not "=" == effect_string[3]
        self.effect_string = effect_string
        if keys is not None:
            self.keys = keys
        self.effects = dict(izip(self.keys, (x.strip() for x in effect_string.split("|", len(self.keys)))))
        self.biotype = self.effects['Transcript_BioType']

    @property
    def consequences(self):
        return list(it.chain.from_iterable(x.split("+") for x in self.effects['Annotation'].split('&')))

    @property
    def coding(self):
        # TODO: check start_gained and utr
        return self.exonic and not "utr" in self.consequence and not "start_gained" in self.consequence


    @property
    def exonic(self):
        csqs = self.consequence
        if isinstance(csqs, basestring):
            csqs = [csqs]
        return any(csq in EXONIC_IMPACTS for csq in csqs) and self.effects['Transcript_BioType'] == 'protein_coding'


class OldSnpEff(SnpEff):

    keys = [x.strip() for x in "Effect | Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_change| Amino_Acid_length | Gene_Name | Gene_BioType | Coding | Transcript | Exon  | ERRORS | WARNINGS".split("|")]

    def __init__(self, effect_string, keys=None, _patt=re.compile("\||\("),
                 prioritize_canonical=False):
        assert not "," in effect_string
        assert not "=" in effect_string
        effect_string = effect_string.rstrip(")")
        self.effect_string = effect_string
        if keys is not None:
            self.keys = keys
        self.effects = dict(izip(self.keys, (x.strip() for x in _patt.split(effect_string))))

    @property
    def consequence(self):
        if '&' in self.effects['Effect']:
            return self.effects['Effect'].split('&')
        return self.effects['Effect']

    @property
    def consequences(self):
        try:
            return [old_snpeff_effect_so.get(c, old_snpeff_effect_so[c.upper()]) for c in it.chain.from_iterable(x.split("+") for x in
                self.effects['Effect'].split('&'))]
        except KeyError:
            return list(it.chain.from_iterable(x.split("+") for x in self.effects['Effect'].split('&')))

    @property
    def severity(self, lookup={'HIGH': 3, 'MED': 2, 'LOW': 1}):
        # higher is more severe. used for ordering.
        try:
            return max(lookup[old_snpeff_lookup[csq]] for csq in self.consequences)
        except KeyError:
            try:
                #in between
                sevs = [IMPACT_SEVERITY.get(csq, "LOW") for csq in self.consequences]
                return max(lookup[s] for s in sevs)
            except KeyError:
                return Effect.severity.fget(self)

    @property
    def is_lof(self):
        return self.biotype == "protein_coding" and self.impact_severity == "HIGH"
