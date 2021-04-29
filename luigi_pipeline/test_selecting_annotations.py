from typing import List

import logging
from functools import wraps
from datetime import datetime
from collections import defaultdict
from contextlib import contextmanager
from inspect import getmembers, ismethod

import hail as hl

# requires --pyfiles "../hail_scripts"
from hail_scripts.v02.utils.computed_fields import variant_id
from hail_scripts.v02.utils.computed_fields import vep


logger = logging.getLogger(__name__)

inp_path = f'gs://cpg-seqr-testing/NA12878_trio/mfranklin/NA12878_trio.vep-annotated.chr20.mt'
output_path = f'gs://cpg-seqr-testing/NA12878_trio/mfranklin/NA12878_trio.test01.mt'
genome_version = "38"

ref_path = "gs://seqr-reference-data/GRCh38/all_reference_data/combined_reference_data_grch38.ht"
clinvar_path = "gs://seqr-reference-data/GRCh38/clinvar/clinvar.GRCh38.2020-06-15.ht"


class RowAnnotationOmit(Exception):
    pass


class RowAnnotation:
    def __init__(self, fn, name=None, requirements: List[str]=None):
        self.fn = fn
        self.name = name or fn.__name__
        self.requirements = requirements

    def __repr__(self):
        requires = None
        if self.requirements:
            requires = f' (requires: {", ".join(self.requirements)})'
        return f"{self.name}{requires}"

    # probably shouldn't let this happen
    # def __call__(self, *args, **kwargs):
    #     return self.fn(*args, **kwargs)


def row_annotation(name=None, fn_require=None):
    """
    Function decorator for methods in a subclass of BaseMTSchema.
    Allows the function to be treated like an row_annotation with annotation name and value.

        @row_annotation()
        def a(self):
            return 'a_val'

        @row_annotation(name='b', fn_require=a)
        def b_1(self):
            return 'b_val'

    Will generate a mt with rows of {a: 'a_val', 'b': 'b_val'} if the function is called.
    TODO: Consider changing fn_require to be a list of requirements.

    When calling the function with annotation already set in the MT, the default behavior is to
    skip unless an overwrite=True is passed into the call.

    :param name: name in the final MT. If not provided, uses the function name.
    :param fn_require: method name strings in class that are dependencies.
    :return:
    """
    def mt_prop_wrapper(func):
        requirements = None
        if fn_require:
            requirements = [fn_require.name]

        c = RowAnnotation(func, name=name, requirements=requirements)
        return c

    return mt_prop_wrapper


class SeqrVariantsAndGenotypesSchema:

    def __init__(self, mt, *, ref_data, clinvar_data, hgmd_data):
        self.mt = mt
        self.mt_instance_meta = {
            'row_annotations': defaultdict(lambda: {
                'annotated': 0,
                'result': {},
            })
        }

        self._ref_data = ref_data
        self._selected_ref_data = ref_data[mt.row_key] if ref_data else None
        self._clinvar_data = clinvar_data
        self._hgmd_data = hgmd_data

    def all_annotation_fns(self) -> List[RowAnnotation]:
        """
        Get all row_annotation decorated methods using introspection.
        :return: list of all annotation functions
        """
        return [a[1] for a in getmembers(self, lambda x: isinstance(x, RowAnnotation))]

    def annotate_all(self, overwrite=False):
        """
        Iterate over all annotation functions and call them on the instance.
        :return: instance object
        """
        called_annotations = set()
        rounds: List[List[RowAnnotation]] = [self.all_annotation_fns()]
        print(f'Trying to apply {len(rounds)} row annotations')
        while len(rounds) > 0:
            rnd = rounds.pop(0)
            print(f'Starting round with {len(rnd)} annotations')
            # add callers that you can't run yet to this list
            next_round = []
            annotations_to_apply = {}
            for annotation in rnd:
                # apply each atn_fn here
                instance_metadata = self.mt_instance_meta['row_annotations'][annotation.name]
                if instance_metadata['annotated'] > 0:
                    # already called
                    continue

                # MT already has annotation, so only continue if overwrite requested.
                if annotation.name in self.mt.rows()._fields:
                    logger.warning(
                        'MT using schema class %s already has "%s" annotation.' % (self.__class__.__name__, annotation.name))
                    if not overwrite:
                        continue
                    logger.info('Overwriting matrix table annotation %s' % annotation.name)

                if annotation.requirements and any(r not in called_annotations for r in annotation.requirements):
                    # skip, and do this next time
                    next_round.append(annotation)
                    continue

                try:
                    # evaluate the function
                    func_ret = annotation.fn(self)
                    # pass
                # Do not annotate when RowAnnotationOmit raised.
                except RowAnnotationOmit:
                    print(f'Received RowAnnotationOmit for "{annotation.name}"')
                    continue

                annotations_to_apply[annotation.name] = func_ret

                instance_metadata['annotated'] += 1
                instance_metadata['result'] = func_ret

            # update the mt
            print('Applying annotations: ' + ', '.join(annotations_to_apply.keys()))
            self.mt = self.mt.annotate_rows(**annotations_to_apply)
            if self._ref_data:
                self._selected_ref_data = self._ref_data[self.mt.row_key]

            called_annotations = called_annotations.union(set(annotations_to_apply.keys()))

            if len(next_round) > 0:
                rounds.append(next_round)

        return self

    def select_annotated_mt(self):
        """
        Returns a matrix table with an annotated rows where each row annotation is a previously called
        annotation (e.g. with the corresponding method or all in case of `annotate_all`).
        :return: a matrix table
        """
        # Selection field is the annotation name of any function that has been called.
        select_fields = []
        for annotation in self.all_annotation_fns():
            if annotation.name in self.mt_instance_meta['row_annotations']:
                inst_fn_metadata = self.mt_instance_meta['row_annotations'][annotation.name]
            else:
                continue

            # Not called.
            if inst_fn_metadata['annotated'] <= 0:
                continue

            select_fields.append(annotation.name)
        return self.mt.select_rows(*select_fields)

    @row_annotation(name='AC')
    def ac(self):
        return self.mt.info.AC[self.mt.a_index - 1]

    @row_annotation(name='AF')
    def af(self):
        return self.mt.info.AF[self.mt.a_index - 1]

    @row_annotation(name='AN')
    def an(self):
        return self.mt.info.AN

    @row_annotation()
    def vep(self):
        return self.mt.vep

    @row_annotation()
    def rsid(self):
        return self.mt.rsid

    @row_annotation()
    def filters(self):
        return self.mt.filters

    @row_annotation()
    def aIndex(self):
        return self.mt.a_index

    @row_annotation()
    def wasSplit(self):
        return self.mt.was_split

    # # Genotypes

    @row_annotation()
    def genotypes(self):
        return hl.agg.collect(hl.struct(**self._genotype_fields()))

    @row_annotation(fn_require=genotypes)
    def samples_no_call(self):
        return self._genotype_filter_samples(lambda g: g.num_alt == -1)

    @row_annotation(fn_require=genotypes)
    def samples_num_alt(self, start=1, end=3, step=1):
        return hl.struct(**{
            '%i' % i: self._genotype_filter_samples(lambda g: g.num_alt == i)
            for i in range(start, end, step)
        })

    @row_annotation(fn_require=genotypes)
    def samples_gq(self, start=0, end=95, step=5):
        # struct of x_to_y to a set of samples in range of x and y for gq.
        return hl.struct(**{
            '%i_to_%i' % (i, i+step): self._genotype_filter_samples(lambda g: ((g.gq >= i) & (g.gq < i+step)))
            for i in range(start, end, step)
        })

    @row_annotation(fn_require=genotypes)
    def samples_ab(self, start=0, end=45, step=5):
        # struct of x_to_y to a set of samples in range of x and y for ab.
        return hl.struct(**{
            '%i_to_%i' % (i, i+step): self._genotype_filter_samples(
                lambda g: ((g.num_alt == 1) & ((g.ab*100) >= i) & ((g.ab*100) < i+step))
            )
            for i in range(start, end, step)
        })


    def _genotype_filter_samples(self, filter):
        # Filter on the genotypes.
        return hl.set(self.mt.genotypes.filter(filter).map(lambda g: g.sample_id))

    def _genotype_fields(self):
        # Convert the mt genotype entries into num_alt, gq, ab, dp, and sample_id.
        is_called = hl.is_defined(self.mt.GT)
        return {
            'num_alt': hl.cond(is_called, self.mt.GT.n_alt_alleles(), -1),
            'gq': hl.cond(is_called, self.mt.GQ, hl.null(hl.tint)),
            'ab': hl.bind(
                lambda total: hl.cond((is_called) & (total != 0) & (hl.len(self.mt.AD) > 1),
                                      hl.float(self.mt.AD[1] / total),
                                      hl.null(hl.tfloat)),
                hl.sum(self.mt.AD)
            ),
            'dp': hl.cond(is_called, hl.int(hl.min(self.mt.DP, 32000)), hl.null(hl.tfloat)),
            'sample_id': self.mt.s
        }

    # # variant IDS
    #
    @row_annotation()
    def originalAltAlleles(self):
        # TODO: This assumes we annotate `locus_old` in this code because `split_multi_hts` drops the proper `old_locus`.
        # If we can get it to not drop it, we should revert this to `old_locus`
        return variant_id.get_expr_for_variant_ids(self.mt.locus_old, self.mt.alleles_old)

    @row_annotation(name='sortedTranscriptConsequences', fn_require=vep)
    def sorted_transcript_consequences(self):
        return vep.get_expr_for_vep_sorted_transcript_consequences_array(self.mt.vep)

    @row_annotation(name='docId')
    def doc_id(self, length=512):
        return variant_id.get_expr_for_variant_id(self.mt, length)

    @row_annotation(name='variantId')
    def variant_id(self):
        return variant_id.get_expr_for_variant_id(self.mt)

    @row_annotation()
    def contig(self):
        return variant_id.get_expr_for_contig(self.mt.locus)

    @row_annotation()
    def pos(self):
        return variant_id.get_expr_for_start_pos(self.mt)

    @row_annotation()
    def start(self):
        return variant_id.get_expr_for_start_pos(self.mt)

    @row_annotation()
    def end(self):
        return variant_id.get_expr_for_end_pos(self.mt)

    @row_annotation()
    def ref(self):
        return variant_id.get_expr_for_ref_allele(self.mt)

    @row_annotation()
    def alt(self):
        return variant_id.get_expr_for_alt_allele(self.mt)

    @row_annotation()
    def xpos(self):
        return variant_id.get_expr_for_xpos(self.mt.locus)

    @row_annotation()
    def xstart(self):
        return variant_id.get_expr_for_xpos(self.mt.locus)

    @row_annotation()
    def xstop(self):
        return variant_id.get_expr_for_xpos(self.mt.locus) + hl.len(variant_id.get_expr_for_ref_allele(self.mt)) - 1

    @row_annotation()
    def rg37_locus(self):
        if self.mt.locus.dtype.reference_genome.name != "GRCh38":
            raise RowAnnotationOmit
        return self.mt.rg37_locus

    # VEP
    @row_annotation(fn_require=sorted_transcript_consequences)
    def domains(self):
        return vep.get_expr_for_vep_protein_domains_set_from_sorted(
            self.mt.sortedTranscriptConsequences)

    @row_annotation(name='transcriptConsequenceTerms', fn_require=sorted_transcript_consequences)
    def transcript_consequence_terms(self):
        return vep.get_expr_for_vep_consequence_terms_set(self.mt.sortedTranscriptConsequences)

    @row_annotation(name='transcriptIds', fn_require=sorted_transcript_consequences)
    def transcript_ids(self):
        return vep.get_expr_for_vep_transcript_ids_set(self.mt.sortedTranscriptConsequences)

    @row_annotation(name='mainTranscript', fn_require=sorted_transcript_consequences)
    def main_transcript(self):
        return vep.get_expr_for_worst_transcript_consequence_annotations_struct(
            self.mt.sortedTranscriptConsequences)

    @row_annotation(name='geneIds', fn_require=sorted_transcript_consequences)
    def gene_ids(self):
        return vep.get_expr_for_vep_gene_ids_set(self.mt.sortedTranscriptConsequences)

    @row_annotation(name='codingGeneIds', fn_require=sorted_transcript_consequences)
    def coding_gene_ids(self):
        return vep.get_expr_for_vep_gene_ids_set(self.mt.sortedTranscriptConsequences, only_coding_genes=True)

    # REF DATA
    @row_annotation()
    def cadd(self):
        return self._selected_ref_data.cadd

    @row_annotation()
    def dbnsfp(self):
        return self._selected_ref_data.dbnsfp

    @row_annotation()
    def geno2mp(self):
        return self._selected_ref_data.geno2mp

    @row_annotation()
    def gnomad_exomes(self):
        return self._selected_ref_data.gnomad_exomes

    @row_annotation()
    def gnomad_exome_coverage(self):
        return self._selected_ref_data.gnomad_exome_coverage

    @row_annotation()
    def gnomad_genomes(self):
        return self._selected_ref_data.gnomad_genomes

    @row_annotation()
    def gnomad_genome_coverage(self):
        return self._selected_ref_data.gnomad_genome_coverage

    @row_annotation()
    def eigen(self):
        return self._selected_ref_data.eigen

    @row_annotation()
    def exac(self):
        return self._selected_ref_data.exac

    @row_annotation()
    def g1k(self):
        return self._selected_ref_data.g1k

    @row_annotation()
    def mpc(self):
        return self._selected_ref_data.mpc

    @row_annotation()
    def primate_ai(self):
        return self._selected_ref_data.primate_ai

    @row_annotation()
    def splice_ai(self):
        return self._selected_ref_data.splice_ai

    @row_annotation()
    def topmed(self):
        return self._selected_ref_data.topmed

    @row_annotation()
    def hgmd(self):
        if self._hgmd_data is None:
            raise RowAnnotationOmit
        return hl.struct(**{'accession': self._hgmd_data[self.mt.row_key].rsid,
                            'class': self._hgmd_data[self.mt.row_key].info.CLASS})

    @row_annotation()
    def clinvar(self):
        clinvar_ref = self._clinvar_data[self.mt.row_key]
        return hl.struct(**{'allele_id': clinvar_ref.info.ALLELEID,
                            'clinical_significance': hl.delimit(clinvar_ref.info.CLNSIG),
                            'gold_stars': clinvar_ref.gold_stars})


@contextmanager
def time_function(name):
    s = datetime.now()
    logger.warning(f"\nStarting {name}")
    yield
    logger.warning(f"\nTIME {(datetime.now() - s).total_seconds()} to {name}")


with time_function("reading reference tables"):
    ref_data = None
    ref_data = hl.read_table(ref_path)
    clinvar = hl.read_table(clinvar_path)
    # hgmd is optional.
    # hgmd = hl.read_table(self.hgmd_ht_path) if self.hgmd_ht_path else None
    hgmd = None

with time_function("Opening VCF"):
    mt = None
    mt = hl.read_matrix_table(inp_path)

with time_function("generating annotation to evaluate"):
    s = SeqrVariantsAndGenotypesSchema(
        mt=mt,
        ref_data=ref_data,
        clinvar_data=clinvar,
        hgmd_data=hgmd
    )
    annotated_mt = s.annotate_all(overwrite=True)
    mt = annotated_mt.select_annotated_mt()

with time_function("annotate_globals to evaluate"):
    mt = mt.annotate_globals(
        genomeVersion=genome_version,
        sampleType="WGS",
        hail_version=hl.version()
    )

# with time_function("evaluate"):
#     mt._force_count_rows()

with time_function("write"):
    mt.write(output_path, overwrite=True)


