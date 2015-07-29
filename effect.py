from functools import total_ordering
import itertools as it

@total_ordering
class Effect(object):

    def __init__(self, effect_dict):
        # or maybe arg should be a dict for Effect()
        raise NotImplemented

    def __le__(self, other):
        return self.severity <= other.severity

    def __eq__(self, other):
        raise NotImplemented

    def __str__(self):
        raise NotImplemented

    def __repr__(self):
        raise NotImplemented

    @property
    def gene(self):
        raise NotImplemented

    @property
    def transcript(self):
        raise NotImplemented

    @property
    def exon(self):
        raise NotImplemented

    @property
    def coding(self):
        return True

    @property
    def lof(self):
        return True

    @property
    def aa_change(self):
        raise NotImplemented

    @property
    def codon_change(self):
        raise NotImplemented

    @property
    def severity(self):
        # higher is more severe. used for ordering.
        return 1

    @property
    def consequence(self):
        raise NotImplemented

    @property
    def biotype(self):
        raise NotImplemented

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
        return self.effects['Annotation']


if __name__ == "__main__":

    s = SnpEff("A|stop_gained|HIGH|C1orf170|ENSG00000187642|transcript|ENST00000433179|protein_coding|3/5|c.262C>T|p.Arg88*|262/3064|262/2091|88/696||")
    print s.effects
    print s.gene, s.transcript, s.consequence
