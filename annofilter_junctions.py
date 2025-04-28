#!/usr/bin/env python3

"""
Python module to annotate splice junctions with features from a GTF file.
Additionally, the module can take a STAR SJ.out.tab file and remove entries that
fail user-specified filters.
"""

__author__ = 'Rob Moccia'
__version__ = '0.1'

import argparse
import os
import transcriptome
import logging
import pybedtools
from collections import defaultdict

def parse_bed12_row(row):
    field_names = ['seqname', 'start', 'end', 'name', 'score', 'strand',
        'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes',
        'blockStarts']
    return {k: v for k,v in zip(field_names, row.strip().split('\t'))}

def group_bed_by_name(bed):
    """
    Helper function to group rows of a bed file by
    the name field. Returns a dictionary with keys
    corresponding to unique names in the bed file
    and a list of corresponding rows as values.
    """
    grouped = defaultdict(lambda: defaultdict(list))
    for line in bed:
        gtf = pybedtools.create_interval_from_list(line[8:])
        gtf.attrs['type'] = line[7]
        side = line[6]
        grouped[line[3]][side].append(gtf)
    return grouped

def get_junction_id(junction):
    """
    Helper function that takes a junction in BED coordinate format and returns a
    junction ID using the GTF coordinate system. The junction ID is a tuple
    of (seqname, start, end) to match the keys of junction dictionaries returned
    by the Transcript and Transcriptome classes.
    It is assumed that the input junction is defined by intronic coordinates.

    Keyword arguments:
    junction (dict) with required keys `seqname`, `start`, and `end`
    """
    return (junction['seqname'], int(junction['start']) + 1, int(junction['end']))

def get_donor_id(junction):
    """
    Helper function that takes a junction in BED coordinate format and returns an
    donor ID using the GTF coordinate system. The donor ID is a tuple
    of (seqname, start, end) or (seqname, end, start), depending on strand, matching the
    keys of donor dictionaries returned by the Transcript and Transcriptome classes.
    It is assumed that the input junction is defined by intronic coordinates.

    Keyword arguments:
    junction -- a dictionary with required keys `seqname`, `start`, and `end`
    """
    if junction['strand'] == '-':
        donor_id = (junction['seqname'], int(junction['end']) + 1, int(junction['end']))
    else:
        donor_id = (junction['seqname'], int(junction['start']), int(junction['start']) + 1)
    return donor_id

def get_acceptor_id(junction):
    """
    Helper function that takes a junction in BED coordinate format and returns an
    acceptor ID using the GTF coordinate system. The acceptor ID is a tuple
    of (seqname, start, end) or (seqname, end, start), depending on strand, matching the
    keys of acceptor dictionaries returned by the Transcript and Transcriptome classes.
    It is assumed that the input junction is defined by intronic coordinates.

    Keyword arguments:
    junction (dict) with required keys `seqname`, `start`, and `end`
    """
    if junction['strand'] == '-':
        acceptor_id = (junction['seqname'], int(junction['start']) + 1, int(junction['start']))
    else:
        acceptor_id = (junction['seqname'], int(junction['end']), int(junction['end']) + 1)
    return acceptor_id

def get_junction_info(junction_id, transcriptome, strand):
    """
    Extract annotations from a Transcriptome class object by junction key

    Keyword args:
    junction_id (tuple) matching the pattern used for keys in transcriptome.junctions
    transcriptome (Transcriptome) a Transcriptome class object
    strand (str) the strand of the junction, either '+', '-', or '.'
    """
    assert(junction_id in transcriptome.junctions)
    transcript_id = transcriptome.junctions[junction_id].copy()
    transcripts = [transcriptome.transcripts[tid] for tid in transcript_id]
    transcript_biotype = [transcript.attributes.get('transcript_biotype') for transcript in transcripts]
    tsl = [transcript.attributes.get('transcript_support_level') for transcript in transcripts]
    gene_id = [transcript.gene_id for transcript in transcripts]
    gene_name = [transcript.attributes.get('gene_name') for transcript in transcripts]
    gene_biotype = [transcript.attributes.get('gene_biotype') for transcript in transcripts]
    strand_match = [strand == transcript.strand
        for transcript in transcripts]
    return {
        'transcripts': transcripts,
        'transcript_id': transcript_id,
        'transcript_biotype': transcript_biotype,
        'tsl': tsl,
        'gene_id': gene_id,
        'gene_name': gene_name,
        'gene_biotype': gene_biotype,
        'strand_match': strand_match
    }

def get_ss_info(ss_id, transcriptome, strand, ss_type):
    """
    Extract annotations from a Transcriptome class object by donor or acceptor key

    Keyword args:
    ss_id (tuple) matching the pattern used for keys in transcriptome.donors and
                  transcriptome.acceptors
    transcriptome (Transcriptome) a Transcriptome class object
    strand (str) the strand of the junction, either '+', '-', or '.'
    """
    if ss_type == 'donor':
        ss_dict = transcriptome.donors
    elif ss_type == 'acceptor':
        ss_dict = transcriptome.acceptors
    else:
        raise TypeError('ss_type must be either donor or acceptor')

    transcript_id = ss_dict[ss_id].copy()
    transcripts = [transcriptome.transcripts[tid] for tid in transcript_id]
    transcript_biotype = [transcript.attributes.get('transcript_biotype') for transcript in transcripts]
    tsl = [transcript.attributes.get('transcript_support_level') for transcript in transcripts]
    gene_id = [transcript.gene_id for transcript in transcripts]
    gene_name = [transcript.attributes.get('gene_name') for transcript in transcripts]
    gene_biotype = [transcript.attributes.get('gene_biotype') for transcript in transcripts]
    strand_match = [strand == transcript.strand for transcript
        in transcripts]
    return {
        'transcripts': transcripts,
        'transcript_id': transcript_id,
        'transcript_biotype': transcript_biotype,
        'tsl': tsl,
        'gene_id': gene_id,
        'gene_name': gene_name,
        'gene_biotype': gene_biotype,
        'strand_match': strand_match
    }

def is_known_junction(junction_id, strand, transcriptome):
    """
    Check if a junction is annotated in the transcriptome. If so,
    return list of SSOverlap objects describing matches. If not,
    return None.

    Keyword args:
    junction_id: tuple junction_id (seqname, donor pos, acceptor pos)
    strand: strand of the junction
    transcriptome: a Transcriptome class object
    """
    if junction_id in transcriptome.junctions.keys():
        match = []
        junc_annot = get_junction_info(junction_id, transcriptome,
            strand)

        for idx in range(len(junc_annot['transcript_id'])):
            match.append(SSOverlap(
                seqname = junc_annot['transcripts'][idx].seqname,
                start = junc_annot['transcripts'][idx].start,
                end = junc_annot['transcripts'][idx].end,
                strand = junc_annot['transcripts'][idx].strand,
                transcript_id = junc_annot['transcript_id'][idx],
                transcript_biotype = junc_annot['transcript_biotype'][idx],
                tsl = junc_annot['tsl'][idx],
                gene_id = junc_annot['gene_id'][idx],
                gene_name = junc_annot['gene_name'][idx],
                gene_biotype = junc_annot['gene_biotype'][idx],
                strand_match = junc_annot['strand_match'][idx]))
        # match = split_by_biotype(junc_annot)
        return match

def is_known_donor(donor_id, strand, transcriptome):
    """
    Check if a donor is annotated in the transcriptome. If so,
    return list of SSOverlap objects describing matches. If not,
    return None.

    Keyword args:
    donor_id: tuple donor_id (seqname, donor pos, acceptor pos)
    strand: strand of the splice site
    transcriptome: a Transcriptome class object
    """
    if donor_id in transcriptome.donors.keys():
        match = []
        annot = get_ss_info(donor_id, transcriptome, strand,
            'donor')

        for idx in range(len(annot['transcript_id'])):
            match.append(SSOverlap(
                seqname = annot['transcripts'][idx].seqname,
                start = annot['transcripts'][idx].start,
                end = annot['transcripts'][idx].end,
                strand = annot['transcripts'][idx].strand,
                transcript_id = annot['transcript_id'][idx],
                transcript_biotype = annot['transcript_biotype'][idx],
                tsl = annot['tsl'][idx],
                gene_id = annot['gene_id'][idx],
                gene_name = annot['gene_name'][idx],
                gene_biotype = annot['gene_biotype'][idx],
                strand_match = annot['strand_match'][idx]))
        return match

def is_known_acceptor(acceptor_id, strand, transcriptome):
    """
    Check if a acceptor is annotated in the transcriptome. If so,
    return list of SSOverlap objects describing matches. If not,
    return None.

    Keyword args:
    acceptor_id: tuple acceptor_id (seqname, acceptor pos, acceptor pos)
    strand: strand of the splice site
    transcriptome: a Transcriptome class object
    """
    if acceptor_id in transcriptome.acceptors.keys():
        match = []
        annot = get_ss_info(acceptor_id, transcriptome, strand,
            'acceptor')

        for idx in range(len(annot['transcript_id'])):
            match.append(SSOverlap(
                seqname = annot['transcripts'][idx].seqname,
                start = annot['transcripts'][idx].start,
                end = annot['transcripts'][idx].end,
                strand = annot['transcripts'][idx].strand,
                transcript_id = annot['transcript_id'][idx],
                transcript_biotype = annot['transcript_biotype'][idx],
                tsl = annot['tsl'][idx],
                gene_id = annot['gene_id'][idx],
                gene_name = annot['gene_name'][idx],
                gene_biotype = annot['gene_biotype'][idx],
                strand_match = annot['strand_match'][idx]))
        return match

def reformat_transcripts(transcript_dict):
    """
    Helper function that takes a dictionary of dictionaries with gene_id
    as top-level keys and minimally a key 'transcripts' in the second level
    dict that contains a list of key-value pairs (i.e., dict) describing the
    transcript, and outputs a dictionary with the latter key-value pairs flattened
    into key: [list of values]. The ith entry of each list corresponds to the same
    transcript. 
    """
    gene_id, gene_name, gene_biotype, strand, strand_match, transcript_id, transcript_biotype, tsl = (
        [] for _ in range(8))
    for gid in transcript_dict:
        for transcript in transcript_dict[gid]['transcripts']:
            if transcript['transcript_id'] not in transcript_id:
                gene_id.append(gid)
                gene_name.append(transcript['gene_name'])
                gene_biotype.append(transcript['gene_biotype'])
                strand.append(transcript['strand'])
                strand_match.append(transcript['strand_match'])
                transcript_id.append(transcript['transcript_id'])
                transcript_biotype.append(transcript['transcript_biotype'])
                tsl.append(transcript['tsl'])
    return {
        'gene_id': gene_id,
        'gene_name': gene_name,
        'gene_biotype': gene_biotype,
        'strand': strand,
        'strand_match': strand_match,
        'transcript_id': transcript_id,
        'transcript_biotype': transcript_biotype,
        'tsl': tsl
    }


# TODO: this class is the start for a potential feature to flag junctions that might
# be introducing a frameshift; the plan is to replace SSOverlap with this although this
# feature might be pushing beyond the intended scope of this program
# test code for class development:
# gtf = '/home/moccir/Share/genomes/annotations/ensembl/release-98/mus_musculus/Mus_musculus.GRCm38.98.gtf'
# transcriptome_model = transcriptome.Transcriptome(gtf)
# # test_bed_line = '1\t3999484\t4007737\tportcullis_pass_5\t0.940\t-\t3999617\t4007655\t255,0,0\t2\t133,82\t0,8171'
# # test_bed_line2 = '1\t5084416\t5089115\tportcullis_pass_92\t1.000\t+\t5084563\t5089008\t255,0,0  2\t147,107\t0,4592'
# # junction = parse_bed12_row(test_bed_line2)
# # junction['start'] = junction.pop('thickStart')
# # junction['end'] = junction.pop('thickEnd')
# # test_annotated_junction = AnnotatedJunction(
# #     **junction)
# # # expected phase
# # [(x.overlap((int(junction['start']), int(junction['end']))), 
# #     x.phase) for x in sorted(transcriptome.transcripts['ENSMUST00000192847'].annotations['CDS'])] 
# class SpliceSite(transcriptome.GtfFeature):
#     """
#     PLACEHOLDER FOR FUTURE FEATURE: NOT WORKING
#     Class to describe a genomic feature intersecting with a splice site
#     """

#     def __init__(self, **kwargs):
#         super().__init__(**kwargs)
#         # self.expected_phase = expected_phase
#         self.features = []

#     # @property
#     # def expected_phase(self):
#     #   return self._expected_phase
    
#     # @expected_phase.setter
#     # def expected_phase(self, value):
#     #   valid_phase = [0, 1, 2, '.', '0', '1', '2']
#     #   if value not in valid_phase:
#     #       raise ValueError(f'expected_phase must be in {valid_phase}')
#     #   self._expected_phase = str(value)

#     @property
#     def features(self):
#         return self._features

#     @features.setter
#     def features(self, value):
#         if isinstance(value, (list, set, tuple)):
#             self._features = value.split()
#         elif isinstance(value, str):
#                 self._features = value
#         else:
#             raise TypeError(f'{value} is not of type str, list, tuple, or set')

#     def add_features(self, values):
#         if isinstance(values, str):
#             values = values.split()
#         if isinstance(values, (list, set, tuple)):
#             for val in values:
#                 if val not in self.features:
#                     self.features.append(val)

#     def phase_diff(self, pos, ss_type):
#         """
#         Return the difference between the expected phase of the
#         feature and the phase of the splice donor/acceptor whose
#         last exonic base is at `pos`. Phase is typically only
#         described for CDS features so self.start and self.end
#         must refer to those coordinates.
#         e.g., If expected_phase is 1 (i.e., 1 base left over
#         contributing to the next codon) and a splice site donor
#         intersecting with this feature ends 5 bases prior to the
#         expected end, the phase_diff is 1. In other words, there
#         is 1 base too few to contribute to the next codon if the
#         reading frame is to be maintained.
#                               .
#         expected: --- ATT ACA G|GT--
#         actual:   --- ATT ACA G|GT--
#                        .
#         """
#         if ss_type == 'donor':
#             if self.end is None:
#                 return None
#             phase_diff = abs(self.end - pos) % 3
#         elif ss_type == 'acceptor':
#             if self.start is None:
#                 return None
#             phase_diff = abs(pos - self.start) % 3
#         else:
#             raise ValueError('ss_type must be "acceptor" or "donor"')
#         return phase_diff


class SSOverlap:
    """
    Class to describe a genomic feature intersecting with a splice site
    """

    def __init__(self, gene_id, gene_biotype, strand_match, transcript_id, transcript_biotype,
            tsl, expected_phase = '.', gene_name=None, seqname=None, start=None,
            end=None, strand=None, **kwargs):
        self.seqname = seqname
        self.start = start
        self.end = end
        self.strand = strand
        self.gene_id = gene_id
        if gene_name is not None:
            self.gene_name = gene_name
        else:
            self.gene_name = gene_id
        self.gene_biotype = gene_biotype
        self.strand_match = strand_match
        self.transcript_id = transcript_id
        self.transcript_biotype = transcript_biotype
        self.tsl = tsl
        self.expected_phase = expected_phase
        self._features = []

    def __repr__(self):
        attr_repr = ', '.join(f"{k}={v!r}" for k, v in self.__dict__.items())
        return f'{self.__class__.__name__} ({attr_repr})'

    def __str__(self):
        return (
            f"{','.join(self.gene_id)}\t"
            f"{','.join(self.gene_name)}\t"
            f"{','.join(self.gene_biotype)}\t"
            f"{','.join(self.strand_match)}\t"
            f"{','.join(self.transcript_id)}\t"
            f"{','.join(self.transcript_biotype)}\t"
            f"{','.join(self.tsl)}\t")

    @property
    def start(self):
        return self._start

    @start.setter
    def start(self, value):
        if value is not None and not isinstance(value, int):
            raise TypeError('start must be of type int or None')
        if value is not None and value < 0:
            raise ValueError('start coordinate must be a positive integer')
        self._start = value

    @property
    def end(self):
        return self._end

    @end.setter
    def end(self, value):
        if value is not None and not isinstance(value, int):
            raise TypeError('end must be of type int or None')
        if value is not None and value < 0:
            raise ValueError('end coordinate must be a positive integer')
        self._end = value

    @property
    def strand(self):
        return self._strand

    @strand.setter
    def strand(self, value):
        valid_strand_args = ['+', '-', '.', None]
        if value not in valid_strand_args:
            raise ValueError(f'strand must be one of {valid_strand_args}')
        self._strand = value

    @property
    def gene_id(self):
        return self._gene_id

    @gene_id.setter
    def gene_id(self, value):
        if isinstance(value, str):
            self._gene_id = value
        else:
            raise TypeError('gene_id argument must be of type str')

    @property
    def gene_name(self):
        return self._gene_name

    @gene_name.setter
    def gene_name(self, value):
        if isinstance(value, str):
            self._gene_name = value
        else:
            raise TypeError('gene_name argument must be of type str')

    @property
    def gene_biotype(self):
        return self._gene_biotype

    @gene_biotype.setter
    def gene_biotype(self, value):
        if isinstance(value, str):
            self._gene_biotype = value
        else:
            raise TypeError('gene_biotype argument must be of type str')

    @property
    def strand_match(self):
        return self._strand_match

    @strand_match.setter
    def strand_match(self, value):
        if isinstance(value, bool):
            self._strand_match = value
        else:
            raise TypeError('strand_match argument must be of type bool')

    @property
    def transcript_id(self):
        return self._transcript_id

    @transcript_id.setter
    def transcript_id(self, value):
        if isinstance(value, str):
            self._transcript_id = value
        else:
            raise TypeError('transcript_id argument must be of type str')

    @property
    def transcript_biotype(self):
        return self._transcript_biotype

    @transcript_biotype.setter
    def transcript_biotype(self, value):
        if isinstance(value, str):
            self._transcript_biotype = value
        else:
            raise TypeError('transcript_biotype argument must be of type str')

    @property
    def tsl(self):
        return self._tsl

    @tsl.setter
    def tsl(self, value):
        if value is None or isinstance(value, (int, str)):
            self._tsl = value
        else:
            raise TypeError(f'{value} is not of type int, str, or None')

    @property
    def expected_phase(self):
        return self._expected_phase
    
    @expected_phase.setter
    def expected_phase(self, value):
        valid_phase = [0, 1, 2, '.', '0', '1', '2']
        if value not in valid_phase:
            raise ValueError(f'expected_phase must be in {valid_phase}')
        self._expected_phase = str(value)

    @property
    def features(self):
        return self._features

    def add_features(self, values):
        if isinstance(values, str):
            values = values.split()
        if isinstance(values, (list, set, tuple)):
            for val in values:
                if val not in self.features:
                    self._features.append(val)

    def phase_diff(self, pos, ss_type):
        """
        Return the difference between the expected phase of the
        feature and the phase of the splice donor/acceptor whose
        last exonic base is at `pos`. Phase is typically only
        described for CDS features so self.start and self.end
        must refer to those coordinates.
        e.g., If expected_phase is 1 (i.e., 1 base left over
        contributing to the next codon) and a splice site donor
        intersecting with this feature ends 5 bases prior to the
        expected end, the phase_diff is 1. In other words, there
        is 1 base too few to contribute to the next codon if the
        reading frame is to be maintained.
                              .
        expected: --- ATT ACA G|GT--
        actual:   --- ATT ACA G|GT--
                       .
        """
        if ss_type == 'donor':
            if self.end is None:
                return None
            phase_diff = abs(self.end - pos) % 3
        elif ss_type == 'acceptor':
            if self.start is None:
                return None
            phase_diff = abs(pos - self.start) % 3
        else:
            raise ValueError('ss_type must be "acceptor" or "donor"')
        return phase_diff

    def to_dict(self):
        return {
            'gene_id': self.gene_id,
            'gene_name': self.gene_name,
            'gene_biotype': self.gene_biotype,
            'strand_match': self.strand_match,
            'transcript_id': self.transcript_id,
            'transcript_biotype': self.transcript_biotype,
            'tsl': self.tsl
        }


class AnnotatedJunction:
    """
    Class to manage information about junction overlaps with a transcriptome
    """

    def __init__(self, seqname=None, start=None, end=None, name=None, score=None,
            strand=None, category=None, donor_known=None, acceptor_known=None,
            **kwargs):
        self.seqname = seqname
        self.start = start
        self.end = end
        self.name = name
        self.score = score
        self.strand = strand
        self._inferred_strand = '.'
        self.category = category
        self.donor_known = donor_known
        self.donor = []
        self.donor_rss = False # possible recursive splice site
        self.acceptor_known = acceptor_known
        self.acceptor = []
        self.acceptor_rss = False # possible recursive splice site
        self.junction_annotation = None
        self.primary_gene_count = None

    def __repr__(self):
        attr_repr = ', '.join(f"{k}={v!r}" for k, v in self.__dict__.items())
        return f'{self.__class__.__name__} ({attr_repr})'


    def __str__(self):
        """
        String representation of junction annotation as a tsv row.
        """
        # make sure any recently added annotations are considered
        self.update_junction_annotation()

        # if no annotation, an empty dictionary will be returned
        primary_strings = defaultdict()
        secondary_strings = defaultdict()

        # NOTE: This is implemented somewhat defensively to protect against duplicate
        # gene or transcript entries. This should no longer be possible since updating
        # the code. However, this works and isn't really hurting anything beyond aesthetics
        # so no time is yet being invested to simplify it.
        if self.junction_annotation:
            primary_gene_idx = [self.junction_annotation['primary']['gene_id'].index(item)
                for item in set(self.junction_annotation['primary']['gene_id'])]
            for k,v in self.junction_annotation['primary'].items():
                if k.startswith('gene') or k.startswith('strand'):
                    string_rep = ','.join([str(item) if item is not None
                        else 'NA' for e,item in enumerate(v) if e
                        in primary_gene_idx])
                else:
                    string_rep = ','.join([str(item) if item is not None
                        else 'NA' for item in v])
                if string_rep == '':
                    primary_strings[k] = 'NA'
                else:
                    primary_strings[k] = string_rep

            secondary_gene_idx = [self.junction_annotation['secondary']['gene_id'].index(item)
                for item in set(self.junction_annotation['secondary']['gene_id'])]
            for k,v in self.junction_annotation['secondary'].items():
                if k.startswith('gene') or k.startswith('strand'):
                    string_rep = ','.join([str(item) if item is not None
                        else 'NA' for e,item in enumerate(v) if e
                        in secondary_gene_idx])
                else:
                    string_rep = ','.join([str(item) if item is not None
                        else 'NA' for item in v])
                if string_rep == '':
                    secondary_strings[k] = 'NA'
                else:
                    secondary_strings[k] = string_rep

        flag_str = ','.join(sorted(list(self.flags)))

        donor_str = self._donor_str()
        acceptor_str = self._acceptor_str()

        return (
            f'{self._bed6_str()}\t'
            f'{self._inferred_strand}\t'
            f'{self.consensus_strand}\t'
            f'{self.category}\t'
            f"{primary_strings.get('gene_id') or 'NA'}\t"
            f"{primary_strings.get('gene_name') or 'NA'}\t"
            f"{primary_strings.get('gene_biotype') or 'NA'}\t"
            f"{primary_strings.get('strand_match') or 'NA'}\t"
            f"{primary_strings.get('transcript_id') or 'NA'}\t"
            f"{primary_strings.get('tsl')}\t"
            f'{self.donor_known}\t'
            f"{donor_str['gene_id']}\t"
            f"{donor_str['gene_name']}\t"
            f"{donor_str['gene_biotype']}\t"
            f"{donor_str['strand_match']}\t"
            f"{donor_str['transcript_id']}\t"
            f"{donor_str['transcript_biotype']}\t"
            f"{donor_str['tsl']}\t"
            f'{self.acceptor_known}\t'
            f"{acceptor_str['gene_id']}\t"
            f"{acceptor_str['gene_name']}\t"
            f"{acceptor_str['gene_biotype']}\t"
            f"{acceptor_str['strand_match']}\t"
            f"{acceptor_str['transcript_id']}\t"
            f"{acceptor_str['transcript_biotype']}\t"
            f"{acceptor_str['tsl']}\t"
            f"{secondary_strings.get('gene_id') or 'NA'}\t"
            f"{secondary_strings.get('gene_name') or 'NA'}\t"
            f"{secondary_strings.get('gene_biotype') or 'NA'}\t"
            f"{secondary_strings.get('strand_match') or 'NA'}\t"
            f"{secondary_strings.get('transcript_id') or 'NA'}\t"
            f"{secondary_strings.get('tsl')}\t"
            f'{flag_str}\n')

    @property
    def donor(self):
        return self._donor

    @donor.setter
    def donor(self, value):
        if not isinstance(value, (tuple, set, list)):
            raise TypeError(
            f'value must be of type tuple, set, or list;'
            f'if trying to add a single donor use .add_donor()')
        if all([isinstance(item, SSOverlap) for item in value]):
            self._donor = value
        else:
            raise TypeError('all entries in donor must be type SSOverlap')

    @property
    def acceptor(self):
        return self._acceptor

    @acceptor.setter
    def acceptor(self, value):
        if not isinstance(value, (tuple, set, list)):
            raise TypeError(
            f'value must be of type tuple, set, or list;'
            f'if trying to add a single acceptor use .add_acceptor()')
        if all([isinstance(item, SSOverlap) for item in value]):
            self._acceptor = value
        else:
            raise TypeError('all entries in acceptor must be type SSOverlap')

    @property
    def donor_rss(self):
        return self._donor_rss

    @donor_rss.setter
    def donor_rss(self, value):
        if isinstance(value, bool):
            self._donor_rss = value
        else:
            raise TypeError('donor_rss must be True or False')

    @property
    def acceptor_rss(self):
        return self._acceptor_rss

    @acceptor_rss.setter
    def acceptor_rss(self, value):
        if isinstance(value, bool):
            self._acceptor_rss = value
        else:
            raise TypeError('acceptor_rss must be True or False')

    @property
    def gene_id(self):
        if self.junction_annotation:
            annot = self.junction_annotation.get('primary')
            return annot['gene_id']
        else:
            return None

    @property
    def gene_name(self):
        if self.junction_annotation:
            annot = self.junction_annotation.get('primary')
            return annot['gene_name']
        else:
            return None

    @property
    def gene_biotype(self):
        if self.junction_annotation:
            annot = self.junction_annotation.get('primary')
            return annot['gene_biotype']
        else:
            return None

    @property
    def strand_match(self):
        if self.junction_annotation:
            annot = self.junction_annotation.get('primary')
            return annot['strand_match']
        else:
            return None

    @property
    def transcript_id(self):
        if self.junction_annotation:
            annot = self.junction_annotation.get('primary')
            return annot['transcript_id']
        else:
            return None

    @property
    def transcript_biotype(self):
        if self.junction_annotation:
            annot = self.junction_annotation.get('primary')
            return annot['transcript_biotype']
        else:
            return None

    @property
    def tsl(self):
        if self.junction_annotation:
            annot = self.junction_annotation.get('primary')
            return annot['tsl']
        else:
            return None

    @property
    def consensus_strand(self):
        """
        Declare a consensus strand for the junction by comparing assigned strand
        (i.e., self.strand) to the inferred strand based on feature overlaps.
        Inferred strand is only used if self.strand == '.' and self.strand takes
        precedence in all other cases.
        """
        if self.strand == '.':
            return self._inferred_strand
        else:
            return self.strand

    @property
    def chimeric(self):
        donor_gid = set(self.get_donor_feature('gene_id'))
        acceptor_gid = set(self.get_acceptor_feature('gene_id'))
        donor_gid.discard(None)
        donor_gid.discard('NA')
        acceptor_gid.discard(None)
        acceptor_gid.discard('NA')
        if donor_gid and acceptor_gid:
            chimeric = len(donor_gid.intersection(
                acceptor_gid)) == 0
        else:
            chimeric = False
        return chimeric

    @property
    def flags(self):
        """
        List of various warnings and other annotations
        """
        flags = set()
        if self.consensus_strand == '.':
            flags.add('strand_ambiguous')
        if self.strand == '.' and self.consensus_strand != '.':
            flags.add('strand_inferred')
        if self._inferred_strand != '.' and self.strand != '.' and self.strand != self._inferred_strand:
            flags.add('inferred_strand_disagrees')
        if self.chimeric:
            flags.add('chimeric')
        if len(self.donor) == 0:
            flags.add('unannotated_donor')
        if len(self.acceptor) == 0:
            flags.add('unannotated_acceptor')
        if self.donor_rss:
            flags.add('donor_rss')
        if self.acceptor_rss:
            flags.add('acceptor_rss')
        if self.primary_gene_count is None:
            logger.warning(f'{self.seqname}_{self.start}:{self.end} {self.name} not analyzed')
            flags.add('junction_not_analyzed')
        elif self.primary_gene_count == 0:
            flags.add('no_primary_gene')
        elif self.primary_gene_count > 1:
            flags.add('multiple_primary_genes')
        if self.primary_biotype_count > 1:
            flags.add('multiple_primary_biotypes')
        return flags

    def update_inferred_strand(self, junction):
        """
        Take a junction in format {strand: {gene: [{transcript}]}} and update
        self._inferred_strand. Use the simple rule that if either '+' or '-' is
        the only key in junction (other than '.'), then that is the inferred strand.
        Otherwise, '.' is the inferred strand (i.e., cannot be inferred).
        """
        possible_strands = set(junction.keys())
        possible_strands.discard('.')
        if len(possible_strands) == 1:
            inferred_strand = possible_strands.pop()
        else:
            inferred_strand = '.'
        self._inferred_strand = inferred_strand
        logger.debug(
            f'update_inferred_strand on {self.name} {self.seqname}_{self.start}_{self.end};'
            f' possibilities = {junction.keys()}; chose {inferred_strand}')

    def add_donor(self, site):
        """
        Add an SSOverlap object to the donor list
        """
        if isinstance(site, SSOverlap):
            self.donor.append(site)
        else:
            raise TypeError('donor must be of type SSOverlap')

    def add_acceptor(self, site):
        """
        Add an SSOverlap object to the acceptor list
        """
        if isinstance(site, SSOverlap):
            self.acceptor.append(site)
        else:
            raise TypeError('acceptor must be of type SSOverlap')

    def get_donor_feature(self, feature):
        return [getattr(item, feature) for item in self.donor
            if isinstance(item, SSOverlap)]

    def get_acceptor_feature(self, feature):
        return [getattr(item, feature) for item in self.acceptor
            if isinstance(item, SSOverlap)]

    def _bed6_fields(self):
        return [self.seqname, self.start, self.end, self.name, self.score, self.strand]

    def _bed6_str(self):
        return '\t'.join(map(str, self._bed6_fields()))

    def left_ss_bed(self):
        fields = self._bed6_fields()
        start = int(fields[1])
        fields[1] = str(start - 1)
        fields[2] = str(start)
        return fields

    def right_ss_bed(self):
        fields = self._bed6_fields()
        end = int(fields[2])
        fields[1] = str(end)
        fields[2] = str(end + 1)
        return fields

    def update_attr(self, overwrite=False, **kwargs):
        """
        Update values for any attributes named as arguments to this function.
        If the attribute is currently set to None it will always be updated.
        If the attribute already has a value, then it will only be overwritten
        if overwrite is set to True.
        """
        # get the values passed by user
        updated_attributes = {k:v for k,v in kwargs.items() if v is not None
            and k not in ['self', 'overwrite']}

        # if attempting to update start or end coords, need to validate new values
        change_start = 'start' in updated_attributes
        change_end = 'end' in updated_attributes
        if change_start and change_end:
            self.validate_start_end(updated_attributes['start'], updated_attributes['end'])
        elif change_start:
            self.validate_start_end(updated_attributes['start'], self.end)
        elif change_end:
            self.validate_start_end(self.start, updated_attributes['end'])

        for attribute, value in updated_attributes.items():
            if value is None:
                continue
            if overwrite or getattr(self, attribute) is None:
                if attribute in ['start', 'end']:
                    value = int(value) if value is not None else None
                setattr(self, attribute, value)
            elif getattr(self, attribute) != value:
                logger.warning(
                    f'AnnotatedJunction {self.name}: '
                    f'ignored request to change {attribute} to {value}; '
                    f'already set to {getattr(self, attribute)} '
                    f'and overwrite=False')

    def _realign_splice_sites(self):
        """
        Returns donors and acceptors in correct orientation. Compares consensus_strand
        to strand and decide if donor and acceptor should be reversed. Return type is a
        dictionary with keys `donor` and `acceptor` containing the SSOverlap objects.
        """
        # inferred_strand only overrules strand when strand is unspecified
        # also, if strand is unspecified, it is treated like + by default
        # therefore, only flip donors and acceptors when strand is unspecified
        # but inferred strand suggests -
        if self.strand == '.' and self.consensus_strand == '-':
            return {
                'donor': self.acceptor,
                'acceptor': self.donor
                }
        else:
            return {
                'donor': self.donor,
                'acceptor': self.acceptor
                }

    def _check_strand_match(self, query_strand):
        """
        Determine if strands of other features passed as query_strand match the
        strand of the junction. Use self.consensus_strand to resolve unspecified
        junction strands.
        """
        if isinstance(query_strand, str):
            query_strand = query_strand.split()
        elif not isinstance(query_strand, (list, set, tuple)):
            raise TypeError('query_strand must be type str, list, set, or tuple')
        if not all([strand in ['+', '-', '.'] for strand in query_strand]):
            raise ValueError(
                'at least one element of query_strand is not a valid strand value')
        return [strand == self.consensus_strand for strand in query_strand]

    def _donor_str(self):
        """
        Merge all donor SSOverlap objects into a single string representation.
        """
        realigned = self._realign_splice_sites()
        donor_list = realigned['donor']

        strand_match = self._check_strand_match([donor.strand for donor
            in donor_list]) or ['NA']

        return {
            'gene_id': ','.join([str(donor.gene_id)
                for donor in donor_list] if donor_list else ['NA']),
            'gene_name': ','.join([str(donor.gene_name)
                for donor in donor_list] if donor_list else ['NA']),
            'gene_biotype': ','.join([str(donor.gene_biotype)
                for donor in donor_list] if donor_list else ['NA']),
            'strand_match': ','.join(str(match) for match in strand_match),
            'transcript_id': ','.join([str(donor.transcript_id)
                for donor in donor_list] if donor_list else ['NA']),
            'transcript_biotype': ','.join([str(donor.transcript_biotype)
                for donor in donor_list] if donor_list else ['NA']),
            'tsl': ','.join([str(donor.tsl) for donor
                in donor_list] if donor_list else ['NA'])
        }

    def _acceptor_str(self):
        """
        Merge all acceptor SSOverlap objects into a single string representation.
        """
        # NOTE: this basically duplicates _donor_str and should arguably be simplified
        # to a single function that also takes 'donor' or 'acceptor' as an argument.
        realigned = self._realign_splice_sites()
        acceptor_list = realigned['acceptor']

        strand_match = self._check_strand_match([acceptor.strand for acceptor
            in acceptor_list]) or ['NA']

        return {
            'gene_id': ','.join([str(acceptor.gene_id)
                for acceptor in acceptor_list] if acceptor_list else ['NA']),
            'gene_name': ','.join([str(acceptor.gene_name)
                for acceptor in acceptor_list] if acceptor_list else ['NA']),
            'gene_biotype': ','.join([str(acceptor.gene_biotype)
                for acceptor in acceptor_list] if acceptor_list else ['NA']),
            'strand_match': ','.join(str(match) for match in strand_match),
            'transcript_id': ','.join([str(acceptor.transcript_id)
                for acceptor in acceptor_list] if acceptor_list else ['NA']),
            'transcript_biotype': ','.join([str(acceptor.transcript_biotype)
                for acceptor in acceptor_list] if acceptor_list else ['NA']),
            'tsl': ','.join([str(acceptor.tsl) for acceptor
                in acceptor_list] if acceptor_list else ['NA'])
        }

    def _get_common_transcripts(self):
        """
        Return transcript IDs shared by donor and acceptor annotations, if any.
        """
        donor_set = set(filter(None, self.get_donor_feature('transcript_id')))
        acceptor_set = set(filter(None, self.get_acceptor_feature('transcript_id')))
        return list(donor_set.intersection(acceptor_set))

    def _merge_common_transcripts(self, common_id):
        """
        Merge donor and acceptor into a junction by shared transcript_id.
        Because they are shared, all information can be extracted from either
        the donor or the acceptor. Therefore, a call to self._realign_splice_sites()
        would not change the result. Arbitrarily, the information is retrieved from the donor.

        common_id: list of transcript IDs
        
        Return: a dict of dicts structured as:
            {strand: {
                gene_id: {
                    gene_name: ..,
                    gene_biotype: ..,
                    transcripts: [{
                        gene_id: ..,
                        gene_name: ..,
                        gene_biotype: ..,
                        strand: ..,
                        transcript_id: ..,
                        transcript_biotype: ..,
                        tsl: ..
                        }]
                    }
                }
            }
        """
        gene_id = []
        gene_name = []
        gene_biotype = []
        strand = []
        transcript_id = []
        transcript_biotype = []
        tsl = []

        transcript_idx = [self.get_donor_feature('transcript_id').index(id) for id in common_id]

        gene_id.extend([self.get_donor_feature('gene_id')[idx] for idx
            in transcript_idx])
        gene_name.extend([self.get_donor_feature('gene_name')[idx] for idx
            in transcript_idx])
        gene_biotype.extend([self.get_donor_feature('gene_biotype')[idx]
            for idx in transcript_idx])
        strand.extend([self.get_donor_feature('strand')[idx]
            for idx in transcript_idx])
        transcript_id.extend([self.get_donor_feature('transcript_id')[idx]
            for idx in transcript_idx])
        transcript_biotype.extend([self.get_donor_feature('transcript_biotype')[idx]
            for idx in transcript_idx])
        tsl.extend([self.get_donor_feature('tsl')[idx] for idx in transcript_idx])

        common_transcripts = defaultdict(dict)
        for query_gene in set(gene_id):
            gene_index = gene_id.index(query_gene)
            transcripts = []
            for i in [e for e,gid in enumerate(gene_id) if gid == query_gene]:
                transcripts.append(
                    {
                        # gene level info is redundant but useful for sanity checks and debugging
                        'gene_id': gene_id[i],
                        'gene_name': gene_name[i],
                        'gene_biotype': gene_biotype[i],
                        'strand': strand[i],
                        'transcript_id': transcript_id[i],
                        'transcript_biotype': transcript_biotype[i],
                        'tsl': tsl[i]
                    }
                )
            gene_dict = {
                'gene_name': gene_name[gene_index],
                'gene_biotype': gene_biotype[gene_index],
                'transcripts': transcripts
            }
            common_transcripts[strand[i]][query_gene] = gene_dict
        return common_transcripts

    def _get_common_genes(self):
        """
        Return gene IDs shared by donor and acceptor annotations, if any.
        """
        donor_set = set(filter(None, self.get_donor_feature('gene_id')))
        acceptor_set = set(filter(None, self.get_acceptor_feature('gene_id')))
        return list(donor_set.intersection(acceptor_set))

    def _merge_common_genes(self, common_id):
        """
        Merge donor and acceptor into a junction by shared gene_id. Because
        they are shared, all information can be extracted from either
        the donor or the acceptor. Therefore, a call to self._realign_splice_sites()
        would not change the result. Arbitrarily, the information is retrieved from the donor.
        Only intended to be called if there are no common transcripts (i.e., _get_common_transcripts
        found no matches) so there will always be no existing transcript info for these genes.
        Instead, unique transcript identifiers of the form novel_{gene_id}_{unique_integer}
        are generated and returned as transcript_id.

        common_id: list of gene IDs
        
        Return: a dict of dicts structured as:
            {strand: {
                gene_id: {
                    gene_name: ..,
                    gene_biotype: ..,
                    transcripts: [{
                        gene_id: ..,
                        gene_name: ..,
                        gene_biotype: ..,
                        strand: ..,
                        transcript_id: ..,
                        transcript_biotype: ..,
                        tsl: ..
                        }]
                    }
                }
            }

        Note that transcript info is uninformative here because there are no known common
        transcripts between the donors and acceptors. The splice junction is predicted to
        be part of some novel transcript.
        """
        gene_id = []
        gene_name = []
        gene_biotype = []
        strand = []

        gene_idx = [self.get_donor_feature('gene_id').index(id) for id in common_id]

        gene_id.extend([self.get_donor_feature('gene_id')[idx] for idx
            in gene_idx])
        gene_name.extend([self.get_donor_feature('gene_name')[idx] for idx
            in gene_idx])
        gene_biotype.extend([self.get_donor_feature('gene_biotype')[idx]
            for idx in gene_idx])
        strand.extend([self.get_donor_feature('strand')[idx]
            for idx in gene_idx])
        # all transcripts are novel by definition because this function is
        # only called if there are no overlaps with annotated transcripts
        transcript_id = [f'novel_{gid}_{count}' if gid != 'NA' else 'NA' for count,gid
            in enumerate([x for x in gene_id])]
        transcript_biotype = ['unknown' if gid != 'NA' else 'NA' for gid
            in [x for x in gene_id]]
        tsl = [None for _ in [x for x in gene_id]]

        common_genes = defaultdict(dict)
        for i in range(len(gene_idx)):
             common_genes[strand[i]][gene_id[i]] = {
                'gene_name': gene_name[i],
                'gene_biotype': gene_biotype[i],
                'transcripts': [{
                    'gene_id': gene_id[i],
                    'gene_name': gene_name[i],
                    'gene_biotype': gene_biotype[i],
                    'strand': strand[i],
                    'transcript_id': transcript_id[i],
                    'transcript_biotype': transcript_biotype[i],
                    'tsl': tsl[i]
                    }]
                }
        return common_genes

    def _merge_donor_acceptor(self):
        """
        Concatenate donors and acceptor feature annotations into a single result. This is
        only intended to be used after _get_common_transcripts() and _get_common_genes()
        have both failed to find matching annotations across the donors and acceptors. It
        is the final catch-all for any feature the junction (or part of the junction if either
        the donor or acceptor does not overlap any annotated gene) might be part of. Features
        are grouped by strand.

        Return: a dict of dicts structured as:
            {strand: {
                gene_id: {
                    gene_name: ..,
                    gene_biotype: ..,
                    transcripts: [{
                        gene_id: ..,
                        gene_name: ..,
                        gene_biotype: ..,
                        strand: ..,
                        transcript_id: ..,
                        transcript_biotype: ..,
                        tsl: ..
                        }]
                    }
                }
            }
        """
        gene_id = []
        gene_name = []
        gene_biotype = []
        transcript_id = []
        transcript_biotype = []
        tsl = []
        strand = []

        if (donor := self.donor):
            donor_transcript_id = [getattr(item, 'transcript_id') for item
                in donor if isinstance(item, SSOverlap)]
            donor_idx = [donor_transcript_id.index(txid) for txid
                in sorted(set(donor_transcript_id))]
            donor_gene_id = [getattr(item, 'gene_id') for item
                in donor if isinstance(item, SSOverlap)]
            gene_id.extend([donor_gene_id[idx] for idx in donor_idx])
            donor_gene_name = [getattr(item, 'gene_name') for item
                in donor if isinstance(item, SSOverlap)]
            gene_name.extend([donor_gene_name[idx] for idx in donor_idx])
            donor_gene_biotype = [getattr(item, 'gene_biotype') for item
                in donor if isinstance(item, SSOverlap)]
            gene_biotype.extend([donor_gene_biotype[idx] for idx in donor_idx])
            transcript_id.extend([donor_transcript_id[idx] for idx in donor_idx])
            donor_transcript_biotype = [getattr(item, 'transcript_biotype') for item
                in donor if isinstance(item, SSOverlap)]
            transcript_biotype.extend([donor_transcript_biotype[idx] for idx in donor_idx])
            donor_tsl = [getattr(item, 'tsl') for item
                in donor if isinstance(item, SSOverlap)]
            tsl.extend([donor_tsl[idx] for idx in donor_idx])
            donor_strand = [getattr(item, 'strand') for item
                in donor if isinstance(item, SSOverlap)]
            strand.extend([donor_strand[idx] for idx in donor_idx])

        if (acceptor := self.acceptor):
            acceptor_transcript_id = [getattr(item, 'transcript_id') for item
                in acceptor if isinstance(item, SSOverlap)]
            acceptor_idx = [acceptor_transcript_id.index(txid) for txid
                in sorted(set(acceptor_transcript_id))]
            acceptor_gene_id = [getattr(item, 'gene_id') for item
                in acceptor if isinstance(item, SSOverlap)]
            gene_id.extend([acceptor_gene_id[idx] for idx in acceptor_idx])
            acceptor_gene_name = [getattr(item, 'gene_name') for item
                in acceptor if isinstance(item, SSOverlap)]
            gene_name.extend([acceptor_gene_name[idx] for idx in acceptor_idx])
            acceptor_gene_biotype = [getattr(item, 'gene_biotype') for item
                in acceptor if isinstance(item, SSOverlap)]
            gene_biotype.extend([acceptor_gene_biotype[idx] for idx in acceptor_idx])
            transcript_id.extend([acceptor_transcript_id[idx] for idx in acceptor_idx])
            acceptor_transcript_biotype = [getattr(item, 'transcript_biotype') for item
                in acceptor if isinstance(item, SSOverlap)]
            transcript_biotype.extend([acceptor_transcript_biotype[idx] for idx in acceptor_idx])
            acceptor_tsl = [getattr(item, 'tsl') for item
                in acceptor if isinstance(item, SSOverlap)]
            tsl.extend([acceptor_tsl[idx] for idx in acceptor_idx])
            acceptor_strand = [getattr(item, 'strand') for item
                in acceptor if isinstance(item, SSOverlap)]
            strand.extend([acceptor_strand[idx] for idx in acceptor_idx])

        merged = defaultdict(dict)
        for query_gene in set(gene_id):
            gene_index = gene_id.index(query_gene)
            transcripts = []
            for i in [e for e,gid in enumerate(gene_id) if gid == query_gene]:
                transcripts.append(
                    {
                        # gene level info is redundant but useful for sanity checks and debugging
                        'gene_id': gene_id[i],
                        'gene_name': gene_name[i],
                        'gene_biotype': gene_biotype[i],
                        'strand': strand[i],
                        'transcript_id': transcript_id[i],
                        'transcript_biotype': transcript_biotype[i],
                        'tsl': tsl[i]
                    }
                )
            gene_dict = {
                'gene_name': gene_name[gene_index],
                'gene_biotype': gene_biotype[gene_index],
                'transcripts': transcripts
            }
            merged[strand[i]][query_gene] = gene_dict
        return merged

    def _solve_junction(self):
        """
        Find the overlapping annotations between the current donor and acceptor
        annotations. This function wraps other internal methods that merge donor
        and acceptor annotations. Return type is therefore of the same format as
        the _merge_* functions:

            {strand: {
                gene_id: {
                    gene_name: ..,
                    gene_biotype: ..,
                    transcripts: [{
                        gene_id: ..,
                        gene_name: ..,
                        gene_biotype: ..,
                        strand: ..,
                        transcript_id: ..,
                        transcript_biotype: ..,
                        tsl: ..
                        }]
                    }
                }
            }
        """
        if not self.donor and not self.acceptor:
            return {}

        # gather list of any common transcripts among the donor/acceptor annotations
        if (common_txs := self._get_common_transcripts()):
            logger.debug(
                f'_solve_junction {self.name}: '
                f'common transcript(s) {common_txs}')
            return self._merge_common_transcripts(common_txs)

        # if no transcripts, gather list of any common genes among the donor/acceptor annotations
        elif (common_genes := self._get_common_genes()):
            logger.debug(
                f'_solve_junction {self.name}: '
                f'common gene(s) {common_genes}')
            return self._merge_common_genes(common_genes)

        # # if no shared genes, report all genes/transcripts potentially involved
        else:
            return self._merge_donor_acceptor()

    def _resolve_primary_gene(self, junction):
        """
        Helper function to prioritize annotations on the consensus strand
        as well as certain gene biotypes more likely to truly be spliced.
        Takes a dictionary of junction annotation information following the
        standard output format of the _merge_* functions (and, therefore, the
        _solve_junction function that calls them as well). The full data structure
        is described in those docstrings but summarized here as: 
        {strand: {
            gene_id: {
                gene_biotype: ..,
                transcripts: [{transcript_id: .., strand: .., etc.}]}
                }
            }
        The second level gene dictionary uses `gene_id` as keys indexing another
        dictionary containing the required key `gene_biotype` as well as other information
        like `transcripts` that is not directly used by this function. Returns a list
        of gene IDs representing all 'spliceable' gene biotypes on the consensus strand.
        """
        SPLICEABLE = [
            'protein_coding',
            'lncRNA',
            'polymorphic_pseudogene',
            'processed_pseudogene',
            'transcribed_unprocessed_pseudogene',
            'translated_unprocessed_pseudogene',
            'transcribed_processed_pseudogene',
            'unprocessed_transcribed_pseudogene',
            'transcribed_unitary_pseudogene',
            'unprocessed_translated_pseudogene',
            'unitary_pseudogene',
            'unprocessed_pseudogene'
        ]

        biotype_dict = defaultdict(list)

        # use dict.get() to prevent adding a new key to defaultdict if it doesn't exist
        if (candidates := junction.get(self.consensus_strand, None)):
            for gene_id, meta in candidates.items():
                biotype_dict[meta['gene_biotype']].append(gene_id)
        else:
            self.primary_gene_count = 0
            self.primary_biotype_count = 0
            return None

        gene_ids = set()
        biotype_count = 0
        for biotype in biotype_dict:
            if biotype in SPLICEABLE:
                gene_ids.update(biotype_dict[biotype])
                biotype_count += 1

        self.primary_gene_count = len(gene_ids)
        self.primary_biotype_count = biotype_count
        return gene_ids

    def update_junction_annotation(self):
        """
        Combine donor and acceptor annotations into a single
        junction annotation. Splits into primary and secondary
        annotations if multiple genes overlap the junction.
        Creates a dictionary with keys 'primary' (annotations
        for the primary gene(s)) and 'secondary' (annotations for all
        other genes, even if not on the consensus strand) and stores it
        in self.junction_annotation. Each key holds a dictionary of
        identical length arrays where the ith index contains the entry for the
        ith gene/transcript annotation. If there are no overlapping annotations,
        the junction_annotation attribute is set to an empty dictionary.
        """
        junction = self._solve_junction()

        if not junction:
            self.primary_gene_count = 0
            self.primary_biotype_count = 0
            self.junction_annotation = {}

        self.update_inferred_strand(junction)

        # check for strand matches against newly determined consensus strand
        for strand in junction:
            strand_match = (self.consensus_strand != '.') and (strand == self.consensus_strand)
            for gene in junction[strand].values():
                for transcript in gene['transcripts']:
                    transcript['strand_match'] = strand_match

        primary_gene_id = self._resolve_primary_gene(junction)

        primary_genes = dict()
        secondary_genes = dict()
        if primary_gene_id:
            correct_strand = junction.pop(self.consensus_strand)
            for gid in primary_gene_id:
                primary_genes[gid] = correct_strand.pop(gid)
            if len(correct_strand) > 0:
                secondary_genes.update(correct_strand)
        # anything left in junction is a secondary gene (e.g., wrong strand, wrong biotype, etc.)
        if len(junction) > 0:
            for key in junction.keys():
                current_strand = junction.get(key)
                secondary_genes.update(current_strand)

        self.junction_annotation = {
            'primary': reformat_transcripts(primary_genes),
            'secondary': reformat_transcripts(secondary_genes)
        }

def make_output_filenames(infile):
    """
    Create output file names based on input filename
    """
    basename = os.path.splitext(infile)[0]
    outfilename = basename + '_annotated.tsv'
    logfilename = basename + '_annotated.log'
    return (outfilename, logfilename)

def annotate_knowns(junctions, transcriptome_model):
    """
    Add annotations for all junctions, donor, and acceptors that are known based
    on the transcriptome model.
    """
    counts = {'DA': 0, 'NDA': 0, 'D': 0, 'A': 0, 'N': 0,
        'potential_rss_donor': 0, 'potential_rss_acceptor': 0,
        'strand_unspecified': 0, 'total': 0}
    annotated_junctions = defaultdict(AnnotatedJunction)
    query_strings = []

    with open(junctions, 'r') as f:
        for line in f:
            # bed can have track lines that need to be skipped
            if not line.startswith('track'):
                counts['total'] += 1
                junction = parse_bed12_row(line)
                # Portcullis specific: the bed file puts the junctions
                # coords in thickStart and thickEnd rather than start and end
                junction['start'] = junction.pop('thickStart')
                junction['end'] = junction.pop('thickEnd')

                annotated_junctions[junction['name']] = AnnotatedJunction(
                    **junction)

                # generate id for the junction
                junction_id = get_junction_id(junction)

                ##### CHECK FOR KNOWN JUNCTIONS: category DA #####
                if (known_junction := is_known_junction(
                        junction_id, junction['strand'], transcriptome_model)):

                    annotated_junctions[junction['name']].update_attr(
                        category = 'DA',
                        donor_known = True,
                        acceptor_known = True)
                    annotated_junctions[junction['name']].donor = known_junction
                    annotated_junctions[junction['name']].acceptor = known_junction

                    counts['DA'] += 1

                else:
                    '''
                    Donor and acceptor could still be known, but not used together
                    in transcriptome annotation.
                    Four conditions need to be checked:
                        1. donor_id is a known donor
                        2. acceptor_id is a known acceptor
                        3. donor_id is a known acceptor
                        4. acceptor_id is a known donor
                    '''
                    # generate ids for the donor and acceptor
                    donor_id = get_donor_id(junction)
                    acceptor_id = get_acceptor_id(junction)

                    donor_known = is_known_donor(donor_id, junction['strand'], transcriptome_model)
                    acceptor_known = is_known_acceptor(acceptor_id, junction['strand'], transcriptome_model)
                    donor_match_acceptor = is_known_acceptor(donor_id, junction['strand'], transcriptome_model)
                    acceptor_match_donor = is_known_donor(acceptor_id, junction['strand'], transcriptome_model)

                    '''
                    Transcriptomes are complex so various patterns will be
                    encountered:

                    1. donor_known, acceptor_known: category NDA,
                         known donor and acceptor but not used together
                         in the same transcript in reference transcriptome
                    2. donor_known, acceptor_known, donor_match_acceptor:
                         category NDA, but donor might be a recursive splice
                         site
                    3. donor_known, acceptor_known, acceptor_match_donor:
                         category NDA, but acceptor might be a recursive splice
                         site
                    4. all four are true: category NDA and both might be recursive
                         splice sites
                    5. donor_known, all others false: category D
                    6. donor_known, donor_match_acceptor, acceptor_match_donor: 
                         category D, donor and acceptor are possible recursive
                         splice sites
                    7. donor_known, donor_match_acceptor: category D, donor is
                         possible recursive splice site
                    8. donor_known, acceptor_match_donor: category D, acceptor
                         is possible recursive splice site
                    9. acceptor_known, all others false: category A
                    10. acceptor_known, donor_match_acceptor, acceptor_match_donor:
                         category A, donor and acceptor are possible recursive
                         splice sites
                    11. acceptor_known, donor_match_acceptor, category A,
                         donor is possible recursive splice site
                    12. acceptor_known, acceptor_match_donor: category A,
                         acceptor is possible recursive splice site
                    13. all false: category N 
                    14. donor_match_acceptor, all others false: N but
                         donor is possible recursive splice site
                    15. acceptor_match_donor, all others false: N but
                         acceptor is possible recursive splice site
                    16. donor_match_acceptor, acceptor_match_donor: N but
                         donor and acceptor are possible recursive splice sites
                    '''
                    # 1. DONOR KNOWN / ACCEPTOR KNOWN (CATEGORY NDA)
                    if donor_known and acceptor_known:
                        annotated_junctions[junction['name']].update_attr(
                            category = 'NDA',
                            donor_known = True,
                            acceptor_known = True)
                        annotated_junctions[junction['name']].donor = donor_known
                        annotated_junctions[junction['name']].acceptor = acceptor_known
        
                        # update counts
                        counts['NDA'] += 1
                        if junction['strand'] == '.':
                            counts['strand_unspecified'] += 1       

                    # 5. DONOR KNOWN (CATEGORY D)
                    elif donor_known:
                        annotated_junctions[junction['name']].update_attr(
                            category = 'D',
                            donor_known = True,
                            acceptor_known = False)
                        annotated_junctions[junction['name']].donor = donor_known

                        # queue other side for GTF query as acceptor
                        if junction['strand'] == '-':
                            # donor is coming from right so add left
                            ss_bed = annotated_junctions[junction['name']].left_ss_bed()
                            ss_bed.extend(['left', 'acceptor'])
                            query_strings.append('\t'.join(ss_bed))
                        else:
                            ss_bed = annotated_junctions[junction['name']].right_ss_bed()
                            ss_bed.extend(['right', 'acceptor'])
                            query_strings.append('\t'.join(ss_bed))

                        counts['D'] += 1

                    #9. ACCEPTOR KNOWN (CATEGORY A)
                    elif acceptor_known:
                        annotated_junctions[junction['name']].update_attr(
                            category = 'A',
                            donor_known = False,
                            acceptor_known = True)
                        annotated_junctions[junction['name']].acceptor = acceptor_known

                        # queue other side for GTF query as donor
                        if junction['strand'] == '-':
                            # acceptor is on left so add right
                            ss_bed = annotated_junctions[junction['name']].right_ss_bed()
                            ss_bed.extend(['right', 'donor'])
                            query_strings.append('\t'.join(ss_bed))
                        else:
                            ss_bed = annotated_junctions[junction['name']].left_ss_bed()
                            ss_bed.extend(['left', 'donor'])
                            query_strings.append('\t'.join(ss_bed))

                        counts['A'] += 1

                    # 13. NEITHER LEFT ACCEPTOR NOR RIGHT DONOR (CATEGORY N)
                    else:

                        annotated_junctions[junction['name']].update_attr(
                            category = 'N',
                            donor_known = False,
                            acceptor_known = False)

                        # queue both sides for GTF query for both
                        ss_bed_left = annotated_junctions[junction['name']].left_ss_bed()
                        ss_bed_right = annotated_junctions[junction['name']].right_ss_bed()
                        if junction['strand'] == '-':
                            ss_bed_left.extend(['left', 'acceptor'])
                            ss_bed_right.extend(['right', 'donor'])
                        else:
                            ss_bed_left.extend(['left', 'donor'])
                            ss_bed_right.extend(['right', 'acceptor'])
                        query_strings.append('\t'.join(ss_bed_left))
                        query_strings.append('\t'.join(ss_bed_right))

                        counts['N'] += 1

                    # ALL OTHER CATEGORIES 2-4, 6-8, 10-12, 14-16: 
                    # flag potential_rss sites
                    if donor_match_acceptor:
                        counts['potential_rss_donor'] += 1
                        annotated_junctions[junction['name']].donor_rss = True

                    if acceptor_match_donor:
                        counts['potential_rss_acceptor'] += 1
                        annotated_junctions[junction['name']].acceptor_rss = True

                    if junction['strand'] == '.':
                        counts['strand_unspecified'] += 1
    return {
        'annotated_junctions': annotated_junctions,
        'query_strings': query_strings,
        'counts': counts
    }

def query_gtf(query_strings, gtf):
    """
    Overlap unknown splice sites with the GTF file.
    """
    query = pybedtools.BedTool('\n'.join(query_strings),
        from_string=True)
    gtf_bedtool = pybedtools.BedTool(gtf)
    overlaps = query.intersect(gtf_bedtool, wa=True, wb=True)
    return overlaps

def annotate_unknowns(annotated_junctions, grouped_matches):
    """
    Add annotations to the unknown donors and acceptors based on the GTF bed
    overlaps.
    """
    GTF_SORT_ORDER = {
        'transcript': 0,
        'CDS': 1,
        'start_codon': 2,
        'stop_codon': 3,
        'five_prime_utr': 4,
        'three_prime_utr': 5,
        'exon': 6,
        'gene': 7
    }

    for name,junction in grouped_matches.items():
        for side in junction.values():
            collected_overlaps = defaultdict()
            gene_ids = set()
            for interval in sorted(side, key = lambda x: GTF_SORT_ORDER.get(x.fields[2], 999)):
                row = transcriptome.parse_gtf_row(str(interval))
                if interval.fields[2] == 'transcript':
                    try:
                        collected_overlaps[
                            row['attributes']['transcript_id']] = SSOverlap(
                            **row['args'], **row['attributes'],
                            strand_match = (annotated_junctions[name].strand
                                == row['args']['strand']),
                            tsl = row['attributes'].get('transcript_support_level'))
                    except KeyError:
                        logger.error(f"{name}: {row}")
                        raise
                    except TypeError:
                        logger.error(f"{name}: {row}")
                        raise
                    gene_ids.add(row['attributes']['gene_id'])
                elif interval.fields[2] == 'CDS':
                    try:
                        collected_overlaps[row['attributes']['transcript_id']].add_features('CDS')
                        collected_overlaps[row['attributes']['transcript_id']].expected_phase = row['args']['phase']
                    except KeyError:
                        logger.error(f"{name}: CDS without transcript {row}")
                        raise
                elif interval.fields[2] == 'start_codon':
                    collected_overlaps[
                    row['attributes']['transcript_id']].add_features('start_codon')
                elif interval.fields[2] == 'stop_codon':
                    collected_overlaps[
                    row['attributes']['transcript_id']].add_features('stop_codon')
                elif interval.fields[2] == 'five_prime_utr':
                    collected_overlaps[
                    row['attributes']['transcript_id']].add_features('five_prime_utr')
                elif interval.fields[2] == 'three_prime_utr':
                    collected_overlaps[
                    row['attributes']['transcript_id']].add_features('three_prime_utr')
                elif (interval.fields[2] == 'gene'
                        and row['attributes']['gene_id']
                        not in gene_ids):
                    logger.warning(
                        f"{name}: gene without transcripts ({row['attributes']['gene_id']})")
            if row['attributes']['type'] in ['donor', 'both']:
                annotated_junctions[name].donor = [x for x in
                    collected_overlaps.values()]
            elif row['attributes']['type'] in ['acceptor', 'both']:
                annotated_junctions[name].acceptor = [x for x in
                    collected_overlaps.values()]
            else:
                logger.warning(f"{name}: unknown bed overlap type is {row['attributes']['type']}")
    return annotated_junctions

def write_output(annotated_junctions, outfilename, exclude_gene_id=None,
        exclude_gene_name=None, exclude_biotype=None, exclude_seqname=None,
        exclude_flag=None):
    HEADER = ['chrom', 'start', 'end', 'name', 'score', 'strand',
    'inferred_strand', 'consensus_strand', 'category', 'gene_id', 'gene_name',
    'gene_biotype', 'strand_match', 'transcript_id', 'tsl',
    'donor_known', 'donor_gene_id', 'donor_gene_name',
    'donor_gene_biotype', 'donor_strand_match', 'donor_transcript_id',
    'donor_transcript_biotype', 'donor_tsl', 'acceptor_known',
    'acceptor_gene_id', 'acceptor_gene_name', 'acceptor_gene_biotype',
    'acceptor_strand_match', 'acceptor_transcript_id',
    'acceptor_transcript_biotype', 'acceptor_tsl',
    'secondary_gene_id', 'secondary_gene_name',
    'secondary_gene_biotype', 'secondary_strand_match',
    'secondary_transcript_id', 'secondary_tsl', 'flags']

    passing_junctions = set()
    counts = {
        'total_junctions': 0,
        'failed_gene_id': 0,
        'failed_gene_name': 0,
        'failed_biotype': 0,
        'failed_seqname': 0,
        'failed_flags': defaultdict(int),
        'passed': 0,
        'failed': 0,
        'flags': defaultdict(int)
    }
    
    failed_flags = set()
    failed_gene_ids = set()
    failed_gene_names = set()
    failed_biotypes = set()
    failed_seqnames = set()

    with open(outfilename, 'w') as f:
        f.write('\t'.join(HEADER) + '\n')
        for junc in annotated_junctions.values():
            junction_str = str(junc)
            f.write(junction_str)
            counts['total_junctions'] += 1
            for flag in junc.flags:
                if flag.strip() != '':
                    counts['flags'][flag.strip()] += 1
            remove = set()
            if exclude_gene_id and junc.gene_id:
                if (bad_id := any_in_any(exclude_gene_id, junc.gene_id)):
                    counts['failed_gene_id'] += 1
                    remove.update(bad_id)
                    failed_gene_ids.update(bad_id)
            if exclude_gene_name and junc.gene_name:
                if (bad_gene := any_in_any(exclude_gene_name, junc.gene_name)):
                    counts['failed_gene_name'] += 1
                    remove.update(bad_gene)
                    failed_gene_names.update(bad_gene)
            if exclude_biotype and junc.gene_biotype:
                if all([biotype in exclude_biotype for biotype in junc.gene_biotype]):
                    counts['failed_biotype'] += 1
                    remove.update(junc.gene_biotype)
                    failed_biotypes.update(junc.gene_biotype)
            if exclude_seqname and junc.seqname in exclude_seqname:
                counts['failed_seqname'] += 1
                remove.update(junc.seqname)
                failed_seqnames.add(junc.seqname)
            # for flags these must be parsed from the string representation of the AnnotatedJunction
            # object because they are built on-the-fly to ensure they are always up to date with the
            # currently stored donor and acceptor annotations
            # if exclude_flag and (bad_flags := any_in_any(exclude_flag, junction_str.split('\t')[-1].strip().split(','))):
            if exclude_flag and (bad_flags := any_in_any(exclude_flag, junc.flags)):
                for next_flag in bad_flags:
                    counts['failed_flags'][next_flag] += 1
                remove.update(bad_flags)
                failed_flags.update(bad_flags)
            if remove:
                counts['failed'] += 1
                logger.debug(
                    f"{junc.name} removed from STAR: {','.join(remove)}")
            else:
                junction_id = (junc.seqname, int(junc.start) + 1, int(junc.end))
                passing_junctions.add(junction_id)
                counts['passed'] += 1

    return {
        'passing_junctions': passing_junctions,
        'counts': counts,
        'failed_flags': failed_flags,
        'failed_gene_ids': failed_gene_ids,
        'failed_gene_names': failed_gene_names,
        'failed_biotypes': failed_biotypes,
        'failed_seqnames': failed_seqnames
        }

def any_in_any(list1, list2):
    """
    Return any overlapping items in two lists as a set (empty set if no overlaps)
    """
    overlap = set(list1).intersection(set(list2))
    return overlap

def filter_star(star, valid_junctions):
    """
    Filter a STAR SJ.out.tab file to exclude any entries not found in
    valid_junctions.

    Keyword args:
    star (str) path to an SJ.out.tab file from STAR
    valid_junctions (set) tuples containg valid junctions in the form (chr, start, end)
    """
    counts = {'total': 0, 'pass': 0, 'fail': 0}
    outfilename = os.path.splitext(star)[0] + '.filtered.tab'
    with open (outfilename, 'w') as outfile:
        with open(star, 'r') as infile:
            for line in infile:
                counts['total'] += 1
                fields = line.strip().split()
                junction_id = (str(fields[0]), int(fields[1]), int(fields[2]))
                if junction_id in valid_junctions:
                    outfile.write(line)
                    counts['pass'] += 1 
                else:
                    counts['fail'] += 1
    return counts

def main(junctions, gtf, exclude_gene_id, exclude_gene_name, exclude_biotype,
        exclude_seqname, exclude_flag, star, verbosity):
    outfilename, logfilename = make_output_filenames(junctions)
    logging.basicConfig(
        filename = logfilename,
        filemode = 'w',
        level = verbosity.upper(),
        format = '%(asctime)s :: %(levelname)s :: %(message)s',
        datefmt = '%Y-%m-%d %H:%M:%S',
        force = True)
    transcriptome_model = transcriptome.Transcriptome(gtf)
    logger.info('Checking for known splice sites')
    first_pass = annotate_knowns(junctions, transcriptome_model)
    logger.info('Intersecting unknown splice sites with GTF')
    gtf_overlaps = query_gtf(first_pass['query_strings'], gtf)
    logger.info('Grouping GTF overlaps by junction')
    grouped_overlaps = group_bed_by_name(gtf_overlaps)
    logger.info('Adding additional annotations')
    second_pass = annotate_unknowns(first_pass['annotated_junctions'],
        grouped_overlaps)
    logger.info('Writing output')
    final = write_output(annotated_junctions=second_pass, outfilename=outfilename,
        exclude_gene_id=exclude_gene_id, exclude_gene_name=exclude_gene_name,
        exclude_biotype=exclude_biotype, exclude_seqname=exclude_seqname,
        exclude_flag=exclude_flag)
    if star:
        logger.info(f"Filtering {star}")
        logger.info(f"EXCLUDE :: gene_id: {','.join(exclude_gene_id or [''])}")
        logger.info(f"EXCLUDE :: gene_name: {','.join(exclude_gene_name or [''])}")
        logger.info(f"EXCLUDE :: gene_biotype: {','.join(exclude_biotype or [''])}")
        logger.info(f"EXCLUDE :: seqname: {','.join(exclude_seqname or [''])}")
        logger.info(f"EXCLUDE :: flag: {','.join(exclude_flag or [''])}")
        filter_stats = filter_star(star, final['passing_junctions'])
        logger.info(f"filtered STAR output written to {os.path.splitext(star)[0] + '.filtered.tab'}")
        logger.info(f"STAR :: num_input_junctions: {filter_stats['total']}")
        logger.info(f"STAR :: num_kept: {filter_stats['pass']}")
        logger.info(f"STAR :: num_removed: {filter_stats['fail']}")
        logger.info(f"FILTER_STAT :: num_input_junctions: {final['counts']['total_junctions']}")
        logger.info(f"FILTER_STAT :: num_pass: {final['counts']['passed']}")
        logger.info(f"FILTER_STAT :: num_fail_any: {final['counts']['failed']}")
        logger.info(f"FILTER_STAT :: num_fail_gene_id: {final['counts']['failed_gene_id']}")
        logger.info(f"FILTER_STAT :: num_fail_gene_name: {final['counts']['failed_gene_name']}")
        logger.info(f"FILTER_STAT :: num_fail_gene_biotype: {final['counts']['failed_biotype']}")
        logger.info(f"FILTER_STAT :: num_fail_seqname: {final['counts']['failed_seqname']}")
        logger.info(f"FILTER_STAT :: num_fail_flag: {'; '.join([f'{k} - {v}' for k,v in final['counts']['failed_flags'].items()])}")
        logger.debug(f"failed_gene_ids: {final['failed_gene_ids']}")
        logger.debug(f"failed_gene_names: {final['failed_gene_names']}")
        logger.debug(f"failed_biotypes: {final['failed_biotypes']}")
        logger.debug(f"failed_seqnames: {final['failed_seqnames']}")
        logger.debug(f"failed_flags: {final['failed_flags']}")
    logger.info(f"Output written to {outfilename}")
    logger.info(f"Log written to {logfilename}")
    logger.info(f"ANNOTATION_STAT :: num_junctions_annotated: {first_pass['counts']['total']}")
    logger.info(f"ANNOTATION_STAT :: DA: {first_pass['counts']['DA']}")
    logger.info(f"ANNOTATION_STAT :: NDA: {first_pass['counts']['NDA']}")
    logger.info(f"ANNOTATION_STAT :: D: {first_pass['counts']['D']}")
    logger.info(f"ANNOTATION_STAT :: A: {first_pass['counts']['A']}")
    logger.info(f"ANNOTATION_STAT :: N: {first_pass['counts']['N']}")
    logger.info(f"ANNOTATION_STAT :: strand_unspecified: {first_pass['counts']['strand_unspecified']}")
    logger.info(f"ANNOTATION_STAT :: donor_possible_recursive_splice_site: {first_pass['counts']['potential_rss_donor']}")
    logger.info(f"ANNOTATION_STAT :: acceptor_possible_recursive_splice_site: {first_pass['counts']['potential_rss_acceptor']}")
    for flag,num in final['counts']['flags'].items(): 
        logger.info(f"ANNOTATION_STAT :: flag_{flag}: {num}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('junctions', help="Portcullis junctions in bed format")
    parser.add_argument('gtf', help="GTF file of annotations to use")
    parser.add_argument('--star', 
        help='optional STAR SJ.out.tab file to filter based on exclude-* flags',
        required=False)
    parser.add_argument('--exclude-gene-id', action='append', required=False)
    parser.add_argument('--exclude-gene-name', action='append', required=False)
    parser.add_argument('--exclude-biotype', action='append', required=False)
    parser.add_argument('--exclude-seqname', action='append', required=False)
    parser.add_argument('--exclude-flag', action='append', required=False)
    parser.add_argument('-v', '--verbosity', choices=['error', 'warning', 'info', 'debug'],
        default='info', help='logging level')
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    args = parser.parse_args()

    logger = logging.getLogger()
    main(**vars(args))

## START DEBUG CODE ##
    # logger = logging.getLogger()
    # main(junctions='/lustre/workspace/home/moccir/rnaseq/nextflow/test_data/annofilter_debug/inferred_strand/BMPRIIKO_8h_Ab93_0.1nM_S14.portcullis.pass.junctions.bed', 
    #     gtf='/lustre/workspace/home/moccir/rnaseq/nextflow/test_data/annofilter_debug/Homo_sapiens.GRCh38.107.gtf',
    #     exclude_gene_id=None, 
    #     exclude_gene_name=['Myh7'],
    #     exclude_biotype=['transcribed_processed_pseudogene'],
    #     exclude_seqname=None, exclude_flag=['chimeric'],
    #     star='/lustre/workspace/home/moccir/rnaseq/nextflow/test_data/annofilter_debug/inferred_strand/BMPRIIKO_8h_Ab93_0.1nM_S14._pass1_.SJ.out.tab',
    #     verbosity = 'debug')
## END DEBUG CODE
