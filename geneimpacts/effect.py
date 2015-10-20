from functools import total_ordering
import itertools as it

EXONIC_IMPACTS = set(["stop_gained",
                      "stop_lost",
                      "frameshift_variant",
                      "initiator_codon_variant",
                      "inframe_deletion",
                      "inframe_insertion",
                      "missense_variant",
                      "incomplete_terminal_codon_variant",
                      "stop_retained_variant",
                      "synonymous_variant",
                      "coding_sequence_variant",
                      "5_prime_UTR_variant",
                      "3_prime_UTR_variant",
                      "transcript_ablation",
                      "transcript_amplification",
                      "feature_elongation",
                      "feature_truncation"])


IMPACT_SEVERITY = dict([
    ('non_coding_exon_variant', 'LOW'),
    ('incomplete_terminal_codon_variant', 'LOW'),
    ('stop_retained_variant', 'LOW'),
    ('synonymous_variant', 'LOW'),
    ('coding_sequence_variant', 'LOW'),
    ('5_prime_UTR_variant', 'LOW'),
    ('3_prime_UTR_variant', 'LOW'),
    ('intron_variant', 'LOW'),
    ('NMD_transcript_variant', 'LOW'),
    ('nc_transcript_variant', 'LOW'),
    ('upstream_gene_variant', 'LOW'),
    ('downstream_gene_variant', 'LOW'),
    ('intergenic_variant', 'LOW'),
    ('transcript_amplification', 'LOW'),
    ('feature_elongation', 'LOW'),
    ('feature_truncation', 'LOW'),

    ('inframe_deletion', 'MED'),
    ('inframe_insertion', 'MED'),
    ('missense_variant', 'MED'),
    ('splice_region_variant', 'MED'),
    ('mature_miRNA_variant', 'MED'),
    ('regulatory_region_variant', 'MED'),
    ('TF_binding_site_variant', 'MED'),
    ('regulatory_region_ablation', 'MED'),
    ('regulatory_region_amplification', 'MED'),
    ('TFBS_ablation', 'MED'),
    ('TFBS_amplification', 'MED'),

    ('transcript_ablation', 'HIGH'),
    ('splice_acceptor_variant', 'HIGH'),
    ('splice_donor_variant', 'HIGH'),
    ('stop_gained', 'HIGH'),
    ('stop_lost', 'HIGH'),
    ('frameshift_variant', 'HIGH'),
    ('initiator_codon_variant', 'HIGH'),

])

# these are taken from the snpeff manual.
# with MODIFIER => LOW and MEDIUM => MED
IMPACT_SEVERITY_SNPEFF = dict([
    ('chromosome_number_variation', 'HIGH'),
    ('exon_loss_variant', 'HIGH'),
    ('exon_loss', 'HIGH'),
    ('rare_amino_acid_variant', 'HIGH'),
    ('start_lost', 'HIGH'),

    ('3_prime_UTR_truncation', 'MED'),
    ('5_prime_UTR_truncation', 'MED'),
    ('disruptive_inframe_deletion', 'MED'),
    ('disruptive_inframe_insertion', 'MED'),

    ('5_prime_UTR_premature_start_codon_gain_variant', 'LOW'),
    ('conserved_intergenic_variant', 'LOW'),
    ('conserved_intron_variant', 'LOW'),
    ('exon_variant', 'LOW'),
    ('gene_variant', 'LOW'),
    ('intergenic_region', 'LOW'),
    ('intragenic_variant', 'LOW'),
    ('miRNA', 'LOW'),
    ('non_coding_transcript_exon_variant', 'LOW'),
    ('non_coding_transcript_variant', 'LOW'),
    ('start_retained', 'LOW'),
    ('transcript_variant', 'LOW'),
])

# http://uswest.ensembl.org/info/genome/variation/predicted_data.html#consequences
IMPACT_SEVERITY_VEP = dict([
    ('transcript_ablation', 'HIGH'),
    ('splice_acceptor_variant', 'HIGH'),
    ('splice_donor_variant', 'HIGH'),
    ('stop_gained', 'HIGH'),
    ('frameshift_variant', 'HIGH'),
    ('stop_lost', 'HIGH'),
    ('start_lost', 'HIGH'),
    ('transcript_amplification', 'HIGH'),

    ('inframe_insertion', 'MED'),
    ('inframe_deletion', 'MED'),
    ('missense_variant', 'MED'),
    ('protein_altering_variant', 'MED'),
    ('regulatory_region_ablation', 'MED'),
    ('TFBS_ablation', 'MED'),

    ('splice_region_variant', 'LOW'),
    ('incomplete_terminal_codon_variant', 'LOW'),
    ('stop_retained_variant', 'LOW'),
    ('synonymous_variant', 'LOW'),
    ('coding_sequence_variant', 'LOW'),
    ('mature_miRNA_variant', 'LOW'),
    ('5_prime_UTR_variant', 'LOW'),
    ('3_prime_UTR_variant', 'LOW'),
    ('non_coding_transcript_exon_variant', 'LOW'),
    ('intron_variant', 'LOW'),
    ('NMD_transcript_variant', 'LOW'),
    ('non_coding_transcript_variant', 'LOW'),
    ('upstream_gene_variant', 'LOW'),
    ('downstream_gene_variant', 'LOW'),
    ('TFBS_amplification', 'LOW'),
    ('TF_binding_site_variant', 'LOW'),
    ('regulatory_region_amplification', 'LOW'),
    ('feature_elongation', 'LOW'),
    ('regulatory_region_variant', 'LOW'),
    ('feature_truncation', 'LOW'),
    ('intergenic_variant', 'LOW'),
])


# I decided these myself.
IMPACT_SEVERITY_CUSTOM = dict([
    ('sequence_feature', 'LOW'),
    ('transcript', 'LOW'),  # ? snpEff

    # occurs with 'exon_loss' in snpEff
    ('3_prime_UTR_truncation', 'MED'),
    ('3_prime_UTR_truncation+exon_loss', 'MED'),
    ('3_prime_UTR_truncation+exon_loss_variant', 'MED'),
    ('exon_loss', 'MED'),

    ('5_prime_UTR_truncation', 'MED'),
    ('5_prime_UTR_truncation+exon_loss_variant', 'MED'),
    ('non_canonical_start_codon', 'LOW'),
    ('initiator_codon_variant', 'LOW'),
])

IMPACT_SEVERITY.update(IMPACT_SEVERITY_SNPEFF)
IMPACT_SEVERITY.update(IMPACT_SEVERITY_CUSTOM)
IMPACT_SEVERITY.update(IMPACT_SEVERITY_VEP)


@total_ordering
class Effect(object):

    def __init__(self, effect_dict):
        # or maybe arg should be a dict for Effect()
        raise NotImplementedError

    def __le__(self, other):
        if self.is_pseudogene and not other.is_pseudogene:
            return True
        elif other.is_pseudogene and not self.is_pseudogene:
            return False

        if self.coding and not other.coding:
            return False
        elif other.coding and not self.coding:
            return True

        if self.severity != other.severity:
            return self.severity <= other.severity

        # sift higher == more damaing
        if self.sift_value < other.sift_value:
            return True

        # polyphen, lower == more damaging
        if self.polyphen_value > other.polyphen_value:
            return True

        return True
        # TODO: look at transcript length?

    @classmethod
    def top_severity(cls, effects):
        for i, e in enumerate(effects):
            if isinstance(e, basestring):

                effects[i] = cls(e)

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
    def gene(self):
        raise NotImplementedError

    @property
    def transcript(self):
        raise NotImplementedError

    @property
    def exonic(self):
        raise NotImplementedError

    @property
    def coding(self):
        raise NotImplementedError
        return True

    @property
    def lof(self):
        return self.impact_severity == "HIGH" and self.biotype == "protein_coding"
        raise NotImplementedError
        return True

    @property
    def aa_change(self):
        raise NotImplementedError

    @property
    def codon_change(self):
        raise NotImplementedError

    @property
    def severity(self, lookup={'HIGH': 3, 'MED': 2, 'LOW': 1}, sev=IMPACT_SEVERITY):
        # higher is more severe. used for ordering.
        return max(lookup[IMPACT_SEVERITY[csq]] for csq in self.consequences)

    @property
    def impact_severity(self):
        return ['xxx', 'LOW', 'MED', 'HIGH'][self.severity]

    @property
    def consequence(self):
        raise NotImplementedError

    @property
    def biotype(self):
        raise NotImplementedError

    @property
    def is_pseudogene(self): #bool
        return 'pseudogene' in self.biotype

class SnpEff(Effect):

    __slots__ = ('effects', 'effect_string')

    keys = [x.strip() for x in 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'.split("|")]

    def __init__(self, effect_string):
        assert not "," in effect_string
        assert not "=" == effect_string[3]
        self.effect_string = effect_string
        self.effects = dict(it.izip(self.keys, (x.strip() for x in effect_string.split("|", len(self.keys)))))

    @property
    def gene(self):
        return self.effects['Gene_Name'] or None

    @property
    def transcript(self):
        return self.effects['Feature_ID'] or None

    @property
    def exon(self):
        return self.effects['Rank']

    @property
    def consequence(self):
        if '&' in self.effects['Annotation']:
            return self.effects['Annotation'].split('&')
        return self.effects['Annotation']

    @property
    def consequences(self):
        return list(it.chain.from_iterable(x.split("+") for x in self.effects['Annotation'].split('&')))

    @property
    def biotype(self):
        return self.effects['Transcript_BioType']

    @property
    def alt(self):
        return self.effects['Allele']

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

    # not defined in ANN field.
    aa_change = None
    sift = None
    sift_value = None
    sift_class = None
    polyphen = None
    polyphen_value = None
    polyphen_class = None

class VEP(Effect):
    keys = "Consequence|Codons|Amino_acids|Gene|SYMBOL|Feature|EXON|PolyPhen|SIFT|Protein_position|BIOTYPE".split("|")
    def __init__(self, effect_string, keys=None):
        assert not "," in effect_string
        assert not "=" in effect_string
        self.effect_string = effect_string
        if keys is not None: self.keys = keys

        self.effect_string = effect_string
        self.effects = dict(it.izip(self.keys, (x.strip() for x in effect_string.split("|"))))

    @property
    def gene(self):
        return self.effects['SYMBOL'] or self.effects['Gene']

    @property
    def transcript(self):
        return self.effects['Feature']

    @property
    def exon(self):
        return self.effects['EXON']

    @property
    def consequence(self):
        if '&' in self.effects['Consequence']:
            return self.effects['Consequence'].split('&')
        return self.effects['Consequence']

    @property
    def consequences(self):
        return list(it.chain.from_iterable(x.split("+") for x in self.effects['Consequence'].split('&')))

    @property
    def biotype(self):
        return self.effects['BIOTYPE']

    @property
    def alt(self):
        return self.effects.get('ALLELE')

    @property
    def coding(self):
        # what about start/stop_gained?
        return self.exonic and any(csq[1:] != "_prime_UTR_variant" for csq in self.consequences)

    @property
    def exonic(self):
        return any(csq in EXONIC_IMPACTS for csq in self.consequences) and self.effects['BIOTYPE'] == 'protein_coding'

    @property
    def sift(self):
        return self.effects['SIFT']

    @property
    def sift_value(self):
        try:
            return float(self.effects['SIFT'].split("(")[1][:-1])
        except IndexError:
            return None

    @property
    def sift_class(self):
        try:
            return self.effects['SIFT'].split("(")[0]
        except IndexError:
            return None

    @property
    def polyphen(self):
        return self.effects['PolyPhen']

    @property
    def polyphen_value(self):
        try:
            return float(self.effects['PolyPhen'].split('(')[1][:-1])
        except IndexError:
            return None

    @property
    def polyphen_class(self):
        try:
            return self.effects['PolyPhen'].split('(')[0]
        except:
            return None

    @property
    def aa_change(self):
        return self.effects['Amino_acids']

if __name__ == "__main__":

    s = SnpEff("A|stop_gained|HIGH|C1orf170|ENSG00000187642|transcript|ENST00000433179|protein_coding|3/5|c.262C>T|p.Arg88*|262/3064|262/2091|88/696||")
    print s.effects
    print s.gene, s.transcript, s.consequence, s.is_pseudogene
    s = SnpEff("G|splice_donor_variant&intron_variant|HIGH|WASH7P|ENSG00000227232|transcript|ENST00000423562|unprocessed_pseudogene|6/9|n.822+2T>C||||||")
    print s.is_pseudogene
    s = SnpEff("G|missense_variant|MODERATE|OR4F5|ENSG00000186092|transcript|ENST00000335137|protein_coding|1/1|c.338T>G|p.Phe113Cys|338/918|338/918|113/305||")
    print s.coding, s.consequence, s.aa_change
