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


@total_ordering
class Effect(object):

    def __init__(self, effect_dict):
        # or maybe arg should be a dict for Effect()
        raise NotImplementedError

    def __le__(self, other):
        if self.pseudogene and not other.pseudogene:
            return True
        elif other.pseudogene and not self.pseudogene:
            return False
        if self.coding and not other.coding:
            return False
        elif other.coding and not self.coding:
            return True
        return self.severity <= other.severity

    def __eq__(self, other):
        raise NotImplementedError

    def __str__(self):
        raise NotImplementedError

    def __repr__(self):
        raise NotImplementedError

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
        return True

    @property
    def lof(self):
        return True

    @property
    def aa_change(self):
        raise NotImplementedError

    @property
    def codon_change(self):
        raise NotImplementedError

    @property
    def severity(self):
        # higher is more severe. used for ordering.
        return 1

    @property
    def consequence(self):
        raise NotImplementedError

    @property
    def biotype(self):
        raise NotImplementedError

    @property
    def is_pseudogene(self): #bool
        raise NotImplementedError

class SnpEff(Effect):

    keys = [x.strip() for x in 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'.split("|")]

    def __init__(self, effect_string):
        assert not "," in effect_string
        assert not "=" in effect_string
        self.effect_string = effect_string
        self.effects = dict(it.izip(self.keys, (x.strip() for x in effect_string.split("|"))))

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
    def biotype(self):
        return self.effects['Transcript_BioType']

    @property
    def alt(self):
        return self.effects['Allele']

    @property
    def is_pseudogene(self):
        return 'pseudogene' in self.effects['Transcript_BioType']

    @property
    def coding(self):
        # TODO: check start_gained and utr
        return self.exonic and not "utr" in self.consequence and self.consequence != "start_gained"

    @property
    def exonic(self):
        csqs = self.consequence
        if isinstance(csqs, basestring):
            csqs = [csqs]
        return any(csq in EXONIC_IMPACTS for csq in csqs) and self.effects['Transcript_BioType'] == 'protein_coding'

    # not defined in ANN field.
    aa_change = None

    @property
    def sift(self):
        return None

    @property
    def polyphen(self):
        return None

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
        return self.effects['Consequence']

    @property
    def biotype(self):
        return self.effects['BIOTYPE']

    @property
    def alt(self):
        return self.effects.get('ALLELE')

    @property
    def is_pseudogene(self):
        return self.effects['BIOTYPE'] == 'processed_pseudogene'

    @property
    def coding(self):
        # what about start/stop_gained?
        return self.exonic and self.effect_name[1:] != "_prime_UTR_variant"

    def exonic(self):
        return self.consequence in EXONIC_IMPACTS and self.effects['BIOTYPE'] == 'protein_coding'

    @property
    def sift(self):
        return self.effects['SIFT']

    @property
    def polyphen(self):
        return self.effects['PolyPhen']

if __name__ == "__main__":

    s = SnpEff("A|stop_gained|HIGH|C1orf170|ENSG00000187642|transcript|ENST00000433179|protein_coding|3/5|c.262C>T|p.Arg88*|262/3064|262/2091|88/696||")
    print s.effects
    print s.gene, s.transcript, s.consequence, s.is_pseudogene
    s = SnpEff("G|splice_donor_variant&intron_variant|HIGH|WASH7P|ENSG00000227232|transcript|ENST00000423562|unprocessed_pseudogene|6/9|n.822+2T>C||||||")
    print s.is_pseudogene
    s = SnpEff("G|missense_variant|MODERATE|OR4F5|ENSG00000186092|transcript|ENST00000335137|protein_coding|1/1|c.338T>G|p.Phe113Cys|338/918|338/918|113/305||")
    print s.coding, s.consequence, s.aa_change
