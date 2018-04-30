import sys
import os
import gzip
from geneimpacts import SnpEff, VEP, Effect, OldSnpEff, BCFT


HERE = os.path.dirname(__file__)

def test_bug():
    e = sorted([VEP('missense_variant|tTt/tGt|F/C|ENSG00000186092|OR4F5|ENST00000335137|1/1|possibly_damaging(0.568)|deleterious(0)|113/305|protein_coding'),
          VEP("splice_region_variant&non_coding_exon_variant&nc_transcript_variant|||ENSG00000223972|DDX11L1|ENST00000456328|2/3||||processed_transcript")])
    assert e[-1].so == 'missense_variant', e[-1].so

def test_snpeff():

    ann = SnpEff("C|splice_donor_variant&splice_region_variant&splice_region_variant&intron_variant|HIGH|DDX11L1|ENSG00000223972|transcript|ENST00000518655|transcribed_unprocessed_pseudogene|3/3|n.734+2_734+3delAG||||||")

    assert ann.gene == "DDX11L1"
    assert ann.transcript == "ENST00000518655"
    assert ann.biotype == "transcribed_unprocessed_pseudogene", ann.biotype
    assert ann.consequences == 'splice_donor_variant&splice_region_variant&splice_region_variant&intron_variant'.split('&')
    assert ann.severity == 3
    assert ann.impact_severity == "HIGH"
    assert ann.aa_change == ""
    assert ann.exon == '3/3', ann.exon
    assert not ann.coding
    assert ann.is_pseudogene

def test_unused():
    extra = ['YYY']
    keys = VEP.keys + extra
    ann = VEP('missense_variant|tTt/tGt|F/C|ENSG00000186092|OR4F5|ENST00000335137|1/1|possibly_damaging(0.568)|deleterious(0)|113/305|protein_coding|xval|yval', keys=keys)
    assert ann.unused() == extra, ann.unused()

    assert ann.effects['YYY'] == 'yval'



def test_vep():

    ann = VEP('missense_variant|tTt/tGt|F/C|ENSG00000186092|OR4F5|ENST00000335137|1/1|possibly_damaging(0.568)|deleterious(0)|113/305|protein_coding|')
    assert ann.gene == 'OR4F5'
    assert ann.transcript == 'ENST00000335137'
    assert ann.aa_change == "F/C", ann.aa_change
    assert ann.consequences == ['missense_variant']
    assert ann.coding
    assert ann.biotype == "protein_coding"
    assert ann.severity == 2
    assert ann.impact_severity == "MED", ann.impact_severity
    assert not ann.is_pseudogene
    assert ann.polyphen_score == 0.568, ann.polyphen
    assert ann.polyphen_pred == "possibly_damaging", ann.polyphen
    assert ann.sift_score == 0.0, ann.sift
    assert ann.sift_pred == "deleterious", ann.sift
    assert not ann.canonical

def test_vep_canonical():

    ann = VEP('missense_variant|tTt/tGt|F/C|ENSG00000186092|OR4F5|ENST00000335137|1/1|possibly_damaging(0.568)|deleterious(0)|113/305|protein_coding|*', prioritize_canonical=True)
    assert ann.gene == 'OR4F5'
    assert ann.transcript == 'ENST00000335137'
    assert ann.aa_change == "F/C", ann.aa_change
    assert ann.consequences == ['missense_variant']
    assert ann.coding
    assert ann.biotype == "protein_coding"
    assert ann.severity == 2
    assert ann.impact_severity == "MED", ann.impact_severity
    assert not ann.is_pseudogene
    assert ann.polyphen_score == 0.568, ann.polyphen
    assert ann.polyphen_pred == "possibly_damaging", ann.polyphen
    assert ann.sift_score == 0.0, ann.sift
    assert ann.sift_pred == "deleterious", ann.sift
    assert ann.is_canonical

def test_bcfts():
    f = os.path.join(HERE, "bcfts.txt.gz")
    with gzip.open(f, "rt") as fh:
        for csq in (BCFT(l.rstrip()) for l in fh):
            assert csq.severity in (1, 2, 3)
            assert csq.is_pseudogene in (True, False)
            assert csq.coding in (True, False), (csq.coding, csq)
            assert csq.is_exonic in (True, False)


def test_veps():

    f = os.path.join(HERE, "vep-csqs.txt.gz")
    with gzip.open(f, "rt") as veps:
        for csq in (VEP(l.strip()) for l in veps):
            assert csq.severity in (1, 2, 3)
            assert csq.is_pseudogene in (True, False)
            assert csq.coding in (True, False)
            assert isinstance(csq.polyphen_value, float) or csq.polyphen_value is None
            csq.gene
            assert isinstance(csq.sift_value, float) or csq.sift_value is None

def test_snpeffs():
    f = os.path.join(HERE, "snpeff-anns.txt.gz")
    with gzip.open(f, "rt") as anns:
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

def test_canonical_order():
    effects = EFFECTS[:]
    effects.append(VEP("intron_variant&nc_transcript_variant|||ENSG00000223972|DDX11L1|ENST00000450305|||||transcribed_unprocessed_pseudogene|*", prioritize_canonical=True))
    effects = sorted(effects)
    assert effects[-1].is_canonical
    assert effects[0].impact_severity == "LOW"
    assert not effects[0].is_canonical

def test_o2():

    keys = [x.strip() for x in "Effect | Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_change| Amino_Acid_length | Gene_Name | Gene_BioType | Coding | Transcript | Exon  | ERRORS | WARNINGS".split("|")]

    effects = [OldSnpEff(v, keys) for v in "DOWNSTREAM(MODIFIER|||||RP5-902P8.10|processed_transcript|NON_CODING|ENST00000434139|),DOWNSTREAM(MODIFIER|||||RP5-902P8.10|processed_transcript|NON_CODING|ENST00000453732|),INTRON(MODIFIER||||138|SCNN1D|protein_coding|CODING|ENST00000470022|3),INTRON(MODIFIER||||638|SCNN1D|protein_coding|CODING|ENST00000338555|3),INTRON(MODIFIER||||638|SCNN1D|protein_coding|CODING|ENST00000400928|2),INTRON(MODIFIER||||669|SCNN1D|protein_coding|CODING|ENST00000379110|6),INTRON(MODIFIER||||704|SCNN1D|protein_coding|CODING|ENST00000325425|2),INTRON(MODIFIER||||802|SCNN1D|protein_coding|CODING|ENST00000379116|5),INTRON(MODIFIER|||||SCNN1D|nonsense_mediated_decay|CODING|ENST00000379101|5),INTRON(MODIFIER|||||SCNN1D|processed_transcript|CODING|ENST00000467651|3)".split(",")]

    effects = sorted(effects)
    assert effects[-1].gene == "SCNN1D", effects[-1].gene

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
    assert top.impact_severity == "MED"
    assert top.so == "missense_variant"
    #assert top[0].


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
    assert e.aa_change == "T245"
    # note that we choose splice_site_region over synonymous coding
    assert e.is_splicing, e.is_splicing

    assert not e.is_coding

    e = OldSnpEff("intergenic_region(MODIFIER|||n.null_nulldelAAGGAAGG|||||||A",
            keys)
    assert e.consequences != []

def test_regr():
    keys = [x.strip() for x in 'Effect | Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_change| Amino_Acid_length | Gene_Name | Transcript_BioType | Gene_Coding | Transcript_ID | Exon_Rank  | Genotype_Number  | ERRORS | WARNINGS'.split("|")]
    v = OldSnpEff('SPLICE_SITE_REGION+SYNONYMOUS_CODING(LOW|SILENT|acG/acA|T245|1134|ANKS1A|protein_coding|CODING|ENST00000360359|5|A)', keys)
    assert v.consequences == ['splice_region_variant', 'synonymous_variant'], v.consequences
    assert v.severity == 2, v.severity
    assert v.aa_change == 'T245'
    v = OldSnpEff('UPSTREAM(MODIFIER||2771|||PSMB1|processed_transcript|CODING|ENST00000462957||C)', keys)
    assert v.consequences == ['upstream_gene_variant'], v.consequences
    assert v.severity == 1, v.severity

    v = OldSnpEff('NEXT_PROT[maturation_peptide](LOW||||241|PSMB1|protein_coding|CODING|||C)', keys)
    assert v.consequences == ['NEXT_PROT[maturation_peptide]'], v.consequences
    assert v.severity == 1, v.severity

    assert v <= v

def test_aa_change():

    eff = OldSnpEff('NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|Agc/Ggc|S418G|696|C1orf170|protein_coding|CODING|ENST00000433179|3|C)')
    assert eff.aa_change == 'S418G'
    ann = SnpEff('C|missense_variant|MODERATE|C1orf170|ENSG00000187642|transcript|ENST00000433179|protein_coding|3/5|c.1252A>G|p.Ser418Gly|1252/3064|1252/2091|418/696||')
    assert ann.aa_change == 'p.Ser418Gly'

def test_old():
    keys = [x.strip() for x in 'Effect | Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_change| Amino_Acid_length | Gene_Name | Transcript_BioType | Gene_Coding | Transcript_ID | Exon_Rank  | Genotype_Number  | ERRORS | WARNINGS'.split("|")]
    v = OldSnpEff('SPLICE_SITE_REGION+SYNONYMOUS_CODING(LOW|SILENT|acG/acA|T245|1134|ANKS1A|protein_coding|CODING|ENST00000360359|5|A)', keys)
    assert v.so == "splice_region_variant", v.so
    v = OldSnpEff('SYNONYMOUS_CODING+SPLICE_SITE_REGION(LOW|SILENT|acG/acA|T245|1134|ANKS1A|protein_coding|CODING|ENST00000360359|5|A)', keys)
    assert v.so == "splice_region_variant", v.so
    assert v.aa_length == 1134, v.aa_length
    assert v.exon == "5", v.exon
    assert v.codon_change == "acG/acA", v.codon_change
    assert v.transcript == "ENST00000360359", v.transcript

def test_old2():
    keys = [x.strip() for x in 'Effect | Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_change| Amino_Acid_length | Gene_Name | Transcript_BioType | Gene_Coding | Transcript_ID | Exon_Rank  | Genotype_Number  | ERRORS | WARNINGS'.split("|")]
    v = OldSnpEff('SPLICE_SITE_REGION+NON_SYNONYMOUS_CODING(LOW|SILENT|acG/acA|T245|1134|ANKS1A|protein_coding|CODING|ENST00000360359|5|A)', keys)
    assert v.so == "missense_variant", v.so

def test_weird_vep():
    keys = "Consequence|Codons|Amino_acids|Gene|SYMBOL|Feature|EXON|PolyPhen|SIFT|Protein_position|BIOTYPE|CANONICAL|CCDS|RadialSVM_score|RadialSVM_pred|LR_score|LR_pred|CADD_raw|CADD_phred|Reliability_index".split("|")

    csqs = ["?|||117581|TWIST2|NM_001271893.1|1/1||||protein_coding|YES||||||||,non_coding_transcript_exon_variant&non_coding_transcript_variant|||117581|TWIST2|NM_001271893.1_dupl8|1/1||||mRNA|||||||||",
            "non_coding_transcript_exon_variant&non_coding_transcript_variant|||117581|TWIST2|NM_001271893.1_dupl8|1/1||||mRNA|||||||||,?|||117581|TWIST2|NM_001271893.1|1/1||||protein_coding|YES||||||||",
"?|||115286|SLC25A26|NM_173471.3|1/1||||protein_coding|YES||||||||",

    "|||ENSG00000138190|EXOC6|ENST00000260762||||-/804|protein_coding,|||ENSG00000138190|EXOC6|ENST00000371547||||-/820|protein_coding,|||ENSG00000138190|EXOC6|ENST00000443748||||-/701|protein_coding,NMD_transcript_variant|||ENSG00000138190|EXOC6|ENST00000495132||||-/404|nonsense_mediated_decay,|||ENSG00000138190|EXOC6|ENST00000371552||||-/799|protein_coding",
    "|||ENSG00000013503|POLR3B|ENST00000539066||||-/1075|protein_coding,nc_transcript_variant|||ENSG00000013503|POLR3B|ENST00000549195|||||processed_transcript,|||ENSG00000013503|POLR3B|ENST00000549569||||-/170|protein_coding,|||ENSG00000013503|POLR3B|ENST00000228347||||-/1133|",
    "|||ENSG00000147202|DIAPH2|ENST00000373054||||-/1097|protein_coding,|||ENSG00000147202|DIAPH2|ENST00000355827||||-/1096|protein_coding,|||ENSG00000147202|DIAPH2|ENST00000324765||||-/1101|protein_coding,|||ENSG00000147202|DIAPH2|ENST00000373049||||-/1096|protein_coding,|||ENSG00000147202|DIAPH2|ENST00000373061||||-/1101|protein_coding",

]
    import sys
    for cs in csqs:
        for c in cs.split(","):
            v = VEP(c, keys)
            assert v.impact_severity in ('LOW', 'MEDIUM', 'HIGH')

def test_empty_snpeff():

    keys = [x.strip() for x in 'Effect | Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_change| Amino_Acid_length | Gene_Name | Transcript_BioType | Gene_Coding | Transcript_ID | Exon_Rank  | Genotype_Number  | ERRORS | WARNINGS'.split("|")]

    eff = "(MODIFIER||||||||||A|ERROR_CHROMOSOME_NOT_FOUND)"
    v = OldSnpEff(eff, keys)
    assert v.impact_severity == "LOW", v.impact_severity

def test_protein_contact():
    ann = SnpEff('C|protein_protein_contact|HIGH|C1orf170|ENSG00000187642|transcript|ENST00000433179|protein_coding|3/5|c.1252A>G|p.Ser418Gly|1252/3064|1252/2091|418/696||')
    assert ann.impact_severity == "HIGH"

def test_gemini_issue812():
    ann = VEP('protein_altering_variant|caGCAGCAGCAGCAGCAACAGCAG/caA|QQQQQQQQ/Q|ENSG00000204842|ATXN2|ENST00000608853|1/25|||14-21/1153|protein_coding|', keys="Consequence|Codons|Amino_acids|Gene|SYMBOL|Feature|EXON|PolyPhen|SIFT|Protein_position|BIOTYPE|CANONICAL".split("|"))
    assert ann.is_coding

def test_bug_vcf2db_21():
	ann = VEP('synonymous_variant|tcA/tcG|S|ENSG00000186092|OR4F5|ENST00000335137|1/1|||60/305|protein_coding||Low_complexity_(Seg):seg&Transmembrane_helices:TMhelix&Prints_domain:PR00237&Superfamily_domains:SSF81321&Gene3D:1.20.1070.10&hmmpanther:PTHR26451&hmmpanther:PTHR26451:SF72&PROSITE_profiles:PS50262||||ENST00000335137.3:c.180A>G|ENST00000335137.3:c.180A>G(p.%3D)|||-0.817044|0.039', keys="Consequence|Codons|Amino_acids|Gene|SYMBOL|Feature|EXON|PolyPhen|SIFT|Protein_position|BIOTYPE|CANONICAL|DOMAINS|CLIN_SIG".split("|"))

	assert ann.codon_change == "tcA/tcG", ann.codon_change

def test_32():
    keys = "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|REFSEQ_MATCH|SOURCE|GIVEN_REF|USED_REF|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|AA_AF|EA_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|MaxEntScan_alt|MaxEntScan_diff|MaxEntScan_ref|SpliceRegion".split("|")
    s = "-|frameshift_variant&start_lost&start_retained_variant|HIGH|HRNR|ENSG00000197915|Transcript|ENST00000368801|protein_coding|2/3||ENST00000368801.2:c.1del|ENSP00000357791.2:p.Met1?|77/9623|1/8553|1/2850|M/X|Atg/tg|rs34061715&COSM111478||-1||deletion|HGNC|HGNC:20846|YES|1|P1|CCDS30859.1|ENSP00000357791|Q86YZ3||UPI00001D7CAD||Ensembl|T|T||||||0.874|0.7337|0.8818|0.9544|0.9592|0.8875|||0.9028|0.7227|0.8276|0.9554|0.9063|0.9541|0.9411|0.9142|0.9069|0.9592|EUR||0&1|0&1|||||||||,-|intron_variant&non_coding_transcript_variant|MODIFIER|FLG-AS1|ENSG00000237975|Transcript|ENST00000420707|antisense_RNA||1/8|ENST00000420707.5:n.159-25632del|||||||rs34061715&COSM111478||1||deletion|HGNC|HGNC:27913||5||||||||Ensembl|T|T|||||10|0.874|0.7337|0.8818|0.9544|0.9592|0.8875|||0.9028|0.7227|0.8276|0.9554|0.9063|0.9541|0.9411|0.9142|0.9069|0.9592|EUR||0&1|0&1|||||||||,-|intron_variant&non_coding_transcript_variant|MODIFIER|FLG-AS1|ENSG00000237975|Transcript|ENST00000593011|antisense_RNA||1/3|ENST00000593011.5:n.296+54843del|||||||rs34061715&COSM111478||1||deletion|HGNC|HGNC:27913||4||||||||Ensembl|T|T|||||10|0.874|0.7337|0.8818|0.9544|0.9592|0.8875|||0.9028|0.7227|0.8276|0.9554|0.9063|0.9541|0.9411|0.9142|0.9069|0.9592|EUR||0&1|0&1|||||||||,-|frameshift_variant&start_lost&start_retained_variant|HIGH|HRNR|388697|Transcript|NM_001009931.2|protein_coding|2/3||NM_001009931.2:c.1del|NP_001009931.1:p.Met1?|80/9632|1/8553|1/2850|M/X|Atg/tg|rs34061715&COSM111478||-1||deletion|EntrezGene|HGNC:20846|YES||||NP_001009931.1||||rseq_mrna_match|RefSeq|T|T||||||0.874|0.7337|0.8818|0.9544|0.9592|0.8875|||0.9028|0.7227|0.8276|0.9554|0.9063|0.9541|0.9411|0.9142|0.9069|0.9592|EUR||0&1|0&1|||||||||".split(",")
    for e in s:
        eff = VEP(e, keys=keys)
        if not "intron" in e.lower():
            assert eff.impact_severity == "HIGH", (eff.impact_severity, e)
