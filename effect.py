from functools import total_ordering
import itertools as it

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
    def exon(self):
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
    def consequence(self):
        if '&' in self.effects['Annotation']:
            return self.effects['Annotation'].split('&')
        return self.effects['Annotation']

    @property
    def alt(self):
        return self.effects['Allele']

    @property
    def is_pseudogene(self):
        return 'pseudogene' in self.effects['Transcript_BioType']

    @property
    def coding(self):
        return self.effects['Transcript_BioType'] == 'protein_coding'

    # not defined in ANN field.
    aa_change = None

if __name__ == "__main__":

    s = SnpEff("A|stop_gained|HIGH|C1orf170|ENSG00000187642|transcript|ENST00000433179|protein_coding|3/5|c.262C>T|p.Arg88*|262/3064|262/2091|88/696||")
    print s.effects
    print s.gene, s.transcript, s.consequence, s.is_pseudogene
    s = SnpEff("G|splice_donor_variant&intron_variant|HIGH|WASH7P|ENSG00000227232|transcript|ENST00000423562|unprocessed_pseudogene|6/9|n.822+2T>C||||||")
    print s.is_pseudogene
    s = SnpEff("G|missense_variant|MODERATE|OR4F5|ENSG00000186092|transcript|ENST00000335137|protein_coding|1/1|c.338T>G|p.Phe113Cys|338/918|338/918|113/305||")
    print s.coding, s.consequence, s.aa_change
