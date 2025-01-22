import enum


######################################
## Variation selection (annotation) ##
######################################
MAX_SOMATIC_VARIANTS = 500_000
MAX_SOMATIC_VARIANTS_GNOMAD_FILTER = 0.01


####################################
## Filtering exclusion thresholds ##
####################################
# NOTE(SW): rename or organise as data classes; delinate somatic and (possibly) germline constants;
# choice will become more clear once the germline workflow is implemented
# Variant filtering
MIN_AD = 4
MIN_AD_DIFFICULT_REGIONS = 6
MIN_AF = 0.10
MAX_GNOMAD_AF = 0.01
MIN_PON_GERMLINE_AF = 0.20
PON_HIT_THRESHOLD = 5


########################################
## Filtering rescue thresholds, tiers ##
########################################
# Variant rescue
MIN_COSMIC_COUNT_RESCUE = 10
MIN_TCGA_PANCANCER_COUNT_RESCUE = 5
MIN_ICGC_PCAWG_COUNT_RESCUE = 5
CLINVAR_CLINSIGS_RESCUE = {
    'conflicting_interpretations_of_pathogenicity',
    'likely_pathogenic',
    'pathogenic',
    'uncertain_significance',
}
PCGR_TIERS_RESCUE = {
    'TIER_1',
    'TIER_2',
}


################################
## Hypermutated report filter ##
################################
PCGR_TIERS_FILTERING = (
    'TIER_1',
    'TIER_2',
    'TIER_3',
    'TIER_4',
    'NONCODING',
)

VEP_IMPACTS_FILTER = (
    'intergenic',
    'intronic',
    'downstream',
    'upstream',
    'impacts_other',
)

GENOMIC_REGIONS_FILTERING = (
    'difficult',
    'none',
    'giab_conf',
)

HOTSPOT_FIELDS_FILTERING = (
    'SAGE_HOTSPOT',
    'hotspot',
    'PCGR_MUTATION_HOTSPOT',
)

RETAIN_FIELDS_FILTERING = (
    'PANEL',
    *HOTSPOT_FIELDS_FILTERING,
)


##################################################
## VCF FILTER tags and FORMAT, INFO annotations ##
##################################################
# Symbol definitions
class VcfFilter(enum.Enum):

    SAGE_LOWCONF = 'SAGE_lowconf'

    MAX_VARIANTS_NON_PASS = 'max_variants_non_pass'
    MAX_VARIANTS_GNOMAD = 'max_variants_gnomad'
    MAX_VARIANTS_NON_CANCER_GENES = 'max_variants_non_cancer_genes'

    MIN_AF = 'min_AF'
    MIN_AD = 'min_AD'
    MIN_AD_DIFFICULT = 'min_AD_difficult_region'
    MIN_AD_NON_GIAB = 'min_AD_non_GIAB'
    PON = 'PON'
    ENCODE = 'ENCODE'
    GNOMAD_COMMON = 'gnomAD_common'

    PCGR_COUNT_LIMIT = 'PCGR_count_limit'

    @property
    def namespace(self):
        return 'FILTER'


class VcfInfo(enum.Enum):

    TUMOR_AF = 'TUMOR_AF'
    TUMOR_DP = 'TUMOR_DP'
    NORMAL_AF = 'NORMAL_AF'
    NORMAL_DP = 'NORMAL_DP'

    SAGE_HOTSPOT = 'SAGE_HOTSPOT'
    SAGE_NOVEL = 'SAGE_NOVEL'
    SAGE_RESCUE = 'SAGE_RESCUE'

    PCGR_TIER = 'PCGR_TIER'
    PCGR_CSQ = 'PCGR_CSQ'
    PCGR_MUTATION_HOTSPOT = 'PCGR_MUTATION_HOTSPOT'
    PCGR_CLINVAR_CLNSIG = 'PCGR_CLINVAR_CLNSIG'
    PCGR_COSMIC_COUNT = 'PCGR_COSMIC_COUNT'
    PCGR_TCGA_PANCANCER_COUNT = 'PCGR_TCGA_PANCANCER_COUNT'
    PCGR_ICGC_PCAWG_COUNT = 'PCGR_ICGC_PCAWG_COUNT'

    CPSR_FINAL_CLASSIFICATION = 'CPSR_FINAL_CLASSIFICATION'
    CPSR_PATHOGENICITY_SCORE = 'CPSR_PATHOGENICITY_SCORE'
    CPSR_CLINVAR_CLASSIFICATION = 'CPSR_CLINVAR_CLASSIFICATION'
    CPSR_CSQ = 'CPSR_CSQ'

    PON_COUNT = 'PON_COUNT'

    HMF_HOTSPOT = 'HMF_HOTSPOT'

    GIAB_CONF = 'GIAB_CONF'
    ENCODE = 'ENCODE'

    DIFFICULT_BAD_PROMOTER = 'DIFFICULT_bad_promoter'
    DIFFICULT_GC15 = 'DIFFICULT_gc15'
    DIFFICULT_GC70TO75 = 'DIFFICULT_gc70to75'
    DIFFICULT_GC75TO80 = 'DIFFICULT_gc75to80'
    DIFFICULT_GC80TO85 = 'DIFFICULT_gc80to85'
    DIFFICULT_GC80 = 'DIFFICULT_gc80'
    DIFFICULT_LOW_COMPLEXITY_DITR = 'DIFFICULT_low_complexity_ditr'
    DIFFICULT_LOW_COMPLEXITY_QUADTR = 'DIFFICULT_low_complexity_quadtr'
    DIFFICULT_LOW_COMPLEXITY_TANDEMREPEATS = 'DIFFICULT_low_complexity_tandemrepeats'
    DIFFICULT_LOW_COMPLEXITY_TRITR = 'DIFFICULT_low_complexity_tritr'
    DIFFICULT_MAPPABILITY_NONUNIQUE = 'DIFFICULT_mappability_nonunique'
    DIFFICULT_SEGDUP = 'DIFFICULT_segdup'

    GNOMAD_AF = 'gnomAD_AF'

    PCGR_TIER_RESCUE = 'PCGR_TIER_RESCUE'
    SAGE_HOTSPOT_RESCUE = 'SAGE_HOTSPOT_RESCUE'
    CLINICAL_POTENTIAL_RESCUE = 'CLINICAL_POTENTIAL_RESCUE'

    ANN = 'ANN'

    RESCUED_FILTERS_EXISTING = 'RESCUED_FILTERS_EXISTING'
    RESCUED_FILTERS_PENDING = 'RESCUED_FILTERS_PENDING'

    PANEL = 'PANEL'

    @property
    def namespace(self):
        return 'INFO'


class VcfFormat(enum.Enum):

    SAGE_AD = 'SAGE_AD'
    SAGE_AF = 'SAGE_AF'
    SAGE_DP = 'SAGE_DP'
    SAGE_SB = 'SAGE_SB'

    @property
    def namespace(self):
        return 'FORMAT'


# Header definitions
VCF_HEADER_ENTRIES = {

    # FILTER
    VcfFilter.SAGE_LOWCONF: {
        'Description': 'Called but filtered by SAGE',
    },

    VcfFilter.MAX_VARIANTS_NON_PASS: {
        'Description': 'Non-pass variant removed when max variant count is exceeded',
    },
    VcfFilter.MAX_VARIANTS_GNOMAD: {
        'Description': 'Non-hotspot population variant removed when max variant count is exceeded',
    },
    VcfFilter.MAX_VARIANTS_NON_CANCER_GENES: {
        'Description': 'Non-cancer gene variant removed when max variant count is exceeded',
    },

    VcfFilter.MIN_AF: {
        'Description': f'Somatic variant AF < {(MIN_AF * 100):.2f}%',
    },
    VcfFilter.MIN_AD: {
        'Description': f'Somatic variant AD < {MIN_AD}',
    },
    VcfFilter.MIN_AD_DIFFICULT: {
        'Description': (
            f'Somatic variant AD < {MIN_AD_DIFFICULT_REGIONS} and in a difficult to '
            'call region (segmental duplications, GC 0-15%, GC 70-100%, LCR, low mappability mask)'
        ),
    },
    VcfFilter.MIN_AD_NON_GIAB: {
        'Description': (
            f'Somatic variant AD < {MIN_AD_DIFFICULT_REGIONS} and not in a GIAB high '
            'confidence region'
        ),
    },
    VcfFilter.PON: {
        'Description': (
            f'Somatic variant that is in {PON_HIT_THRESHOLD} or more samples from the '
            'panel of normals'
        ),
    },
    VcfFilter.ENCODE: {
        'Description': 'ENCODE blocklist region https://github.com/Boyle-Lab/Blacklist',
    },
    VcfFilter.GNOMAD_COMMON: {
        'Description': f'gnomAD AF >= {MAX_GNOMAD_AF}',
    },

    VcfFilter.PCGR_COUNT_LIMIT: {
        'Description': 'Manually filtered to meet PCGR 500,000 variant limit',
    },

    # INFO
    VcfInfo.TUMOR_AF: {
        'Number': '1',
        'Type': 'Float',
        'Description': 'Tumor sample AF',
    },
    VcfInfo.TUMOR_DP: {
        'Number': '1',
        'Type': 'Integer',
        'Description': 'Tumor sample DP',
    },
    VcfInfo.NORMAL_AF: {
        'Number': '1',
        'Type': 'Float',
        'Description': 'Normal sample AF',
    },
    VcfInfo.NORMAL_DP: {
        'Number': '1',
        'Type': 'Integer',
        'Description': 'Normal sample DP',
    },

    VcfInfo.SAGE_HOTSPOT: {
        'Number': '0',
        'Type': 'Flag',
        'Description': 'Variant called by SAGE in a provided hotspot location',
    },
    VcfInfo.SAGE_NOVEL: {
        'Number': '0',
        'Type': 'Flag',
        'Description': 'Novel variant called by SAGE (AD, AF, DP, and SAGE_SB are taken used from the SAGE call)',
    },
    VcfInfo.SAGE_RESCUE: {
        'Number': '0',
        'Type': 'Flag',
        'Description': 'Variant rescued by a matching SAGE call',
    },

    VcfInfo.PCGR_TIER: {
        'Number': '1',
        'Type': 'String',
        'Description': (
            'Tier reported by PCGR with the following meaning: TIER_1: strong clinical '
            'significance; TIER_2: potential clinical significance; TIER_3: uncertain clinical '
            'significance; TIER_4: other coding variants; NONCODING: other non-coding variants'
        ),
    },
    VcfInfo.PCGR_CSQ: {
        'Number': '.',
        'Type': 'String',
        'Description': (
            'Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|'
            'Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|'
            'CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|'
            'DISTANCE|STRAND|FLAGS|PICK|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|'
            'MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|'
            'UNIPROT_ISOFORM|RefSeq|DOMAINS|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|'
            'gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|'
            'gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|CLIN_SIG|SOMATIC|PHENO|CHECK_REF|'
            'NearestExonJB'
        ),
    },
    VcfInfo.PCGR_MUTATION_HOTSPOT: {
        'Number': '.',
        'Type': 'String',
        'Description': 'Known cancer mutation hotspot, as found in cancerhotspots.org_v2, Gene|Codon|Q-value',
    },
    VcfInfo.PCGR_CLINVAR_CLNSIG: {
        'Number': '.',
        'Type': 'String',
        'Description': 'ClinVar clinical significance',
    },
    VcfInfo.PCGR_COSMIC_COUNT: {
        'Number': '1',
        'Type': 'Integer',
        'Description': 'Count of COSMIC mutations',
    },
    VcfInfo.PCGR_TCGA_PANCANCER_COUNT: {
        'Number': '1',
        'Type': 'Integer',
        'Description': 'Raw variant count across all tumor types',
    },
    VcfInfo.PCGR_ICGC_PCAWG_COUNT: {
        'Number': '1',
        'Type': 'Integer',
        'Description': 'Count of ICGC PCAWG hits',
    },

    VcfInfo.CPSR_FINAL_CLASSIFICATION: {
        'Number': '1',
        'Type': 'String',
        'Description': (
            'Final variant classification based on the combination of CLINVAR_CLASSIFICTION (for '
            'ClinVar-classified variants), and CPSR_CLASSIFICATION (for novel variants)'
        ),
    },
    VcfInfo.CPSR_PATHOGENICITY_SCORE: {
        'Number': '1',
        'Type': 'Float',
        'Description': 'Aggregated CPSR pathogenicity score',
    },
    VcfInfo.CPSR_CLINVAR_CLASSIFICATION: {
        'Number': '1',
        'Type': 'String',
        'Description': 'Clinical significance of variant on a five-tiered scale',
    },
    VcfInfo.CPSR_CSQ: {
        'Number': '.',
        'Type': 'String',
        'Description': (
            'Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|'
            'Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|'
            'Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|'
            'FLAGS|PICK|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|'
            'MANE_PLUS_CLINICAL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|UNIPROT_ISOFORM|RefSeq|'
            'DOMAINS|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|gnomAD_AF|gnomAD_AFR_AF|'
            'gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|'
            'gnomAD_SAS_AF|CLIN_SIG|SOMATIC|PHENO|CHECK_REF|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|'
            'MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|NearestExonJB|LoF|LoF_filter|LoF_flags|'
            'LoF_info'
        ),
    },

    VcfInfo.SAGE_HOTSPOT_RESCUE: {
        'Number': '0',
        'Type': 'Flag',
        'Description': '',
    },
    VcfInfo.PCGR_TIER_RESCUE: {
        'Number': '0',
        'Type': 'Flag',
        'Description': '',
    },
    VcfInfo.CLINICAL_POTENTIAL_RESCUE: {
        'Number': '0',
        'Type': 'Flag',
        'Description': '',
    },

    VcfInfo.ANN: {
        'Number': '.',
        'Type': 'String',
        'Description': (
            'Functional annotations: \'Allele | Annotation | Annotation_Impact | Gene_Name | '
            'Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | '
            'cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | '
            'ERRORS / WARNINGS / INFO\''
        ),
    },

    VcfInfo.RESCUED_FILTERS_EXISTING: {
        'Number': '1',
        'Type': 'String',
        'Description': 'Filters existing prior to variant rescue',
    },

    VcfInfo.RESCUED_FILTERS_PENDING: {
        'Number': '1',
        'Type': 'String',
        'Description': 'Filters pending prior to variant rescue',
    },

    VcfInfo.PANEL: {
        'Number': '0',
        'Type': 'Flag',
        'Description': 'UMCCR somatic panel CDS (2,000 bp padding)',
    },


    # FORMAT
    VcfFormat.SAGE_AD: {
        'Number': 'R',
        'Type': 'Integer',
        'Description': 'Allelic depths for the ref and alt alleles in the order listed',
    },
    VcfFormat.SAGE_AF: {
        'Number': '1',
        'Type': 'Float',
        'Description': 'Allelic frequency calculated from read context counts as (Full + Partial + Core + Realigned + Alt) / Coverage',
    },
    VcfFormat.SAGE_DP: {
        'Number': '1',
        'Type': 'Integer',
        'Description': 'Approximate read depth (reads with MQ=255 or with bad mates are filtered)',
    },
    VcfFormat.SAGE_SB: {
        'Number': '1',
        'Type': 'Float',
        'Description': 'Strand bias - percentage of first-in-pair reads',
    },
}


#####################
##      Other      ##
#####################
CONTIGS_MAIN = [f'chr{e}' for e in [*range(1, 23), 'X', 'Y']]
