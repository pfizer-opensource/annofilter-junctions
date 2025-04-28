#!/usr/bin/env python3

"""
Python module to build an object-oriented transcriptome representation from a GTF file.
"""

__author__ = 'Rob Moccia'
__version__ = '0.1'

from collections import defaultdict
from itertools import dropwhile
import re
import logging

logging.basicConfig(
    level = logging.DEBUG,
    format = '%(asctime)s :: %(name)s :: %(levelname)s :: %(message)s',
    datefmt = '%Y-%m-%d %H:%M:%S',
    force = True)
logger = logging.getLogger(__name__)

def parse_gtf_attributes(attributes, delim=' '):
    """
    Function to convert a GTF attribute field into
    a dictionary. By default, this should work on Ensembl
    GTF files which use a space to separate key/value pairs.
    Gencode uses "=" as a delimiter so this would need to be passed
    as the delim argument as appropriate.

    Primarily called inside class instances but has general utility so
    leaving function outside of class definitions
    """

    '''
    put key:value pairs into a set to protect against duplicated entries
    even though this should never happen
    '''
    kv_pairs = set()

    tags = attributes.strip().split(';')

    for tag in tags:
        try:
            k,v = tag.strip().split(delim)
        except ValueError:
            # first trivial possibility is that the attributes field
            # ended with a semicolon leading to a blank space that can
            # be safely ignored
            if tag == '':
                continue
            else:
                # try to handle violations of the format
                # like mixing delimiters -- try to split on "="
                # (note this is terrible practice and really the
                # gtf itself should be fixed)
                try:
                    k,v = tag.strip().split('=')
                except ValueError:
                # Ensembl technically encloses the value in double quotes
                # so it might contain additional spaces leading to additional
                # values with the simple .split() approach
                    try:
                        k = tag.strip().split(delim)[0]
                        v = re.findall(r'"([^"]*)"', tag)[0]
                    except:
                        raise
                except:
                    raise
        kv_pairs.add((k.strip(), v.strip().strip('"')))

    return {k:v for k,v in kv_pairs}

def parse_gtf_row(row):
    """
    Take a row from a GTF file and return a dictionary of dictionaries.
    This can be useful for parsing a GTF file into GtfFeature objects using
    **kwargs.
    The args dictionary contains the first 8 columns as key:value pairs.
    The attributes dictionary contains all parsed attributes found in the 9th
    column.

    Primarily called inside the Transcriptome class but has general utility so
    leaving function outside of class definitions
    """
    field_names = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    args = {k: v for k,v in zip(field_names, row.strip().split('\t'))}
    args['start'] = int(args['start'])
    args['end'] = int(args['end'])
    attributes = parse_gtf_attributes(args.pop('attributes'))
    return {'args': args, 'attributes': attributes}


class GtfFeature:
    """Class representation of a generic GTF file feature"""

    '''
    set defaults for all arguments so that a blank object can be easily created
    the update() method will allow all of this to be filled in later as information
    is encountered while parsing a GTF file
    '''
    def __init__(
            self, seqname=None, source=None, feature=None, start=None, end=None,
            score=None, strand=None, phase='.', attributes=None, gene_id=None,
            transcript_id=None, **kwargs):
        super().__init__(**kwargs)
        self.seqname = seqname
        self.source = source
        self.feature = feature
        self.start = int(start) if start is not None else None
        self.end = int(end) if end is not None else None
        self.score = score
        self.strand = strand
        self.phase = phase
        self.attributes = attributes
        self.gene_id = gene_id
        self.transcript_id = transcript_id

    @property
    def start(self):
        return self._start

    @start.setter
    def start(self, value):
        if value is not None and not isinstance(value, int):
            raise TypeError('start must be of type int or None')
        if value is not None and value < 0:
            raise ValueError('start coordinate must be a positive integer')
        if hasattr(self, '_end'):
            if value is not None and self.end is not None and value > self.end:
                raise ValueError('start coordinate must be <= end coordinate')
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
        if hasattr(self, '_start'):
            if value is not None and self.start is not None and value < self.start:
                raise ValueError('end coordinate must be >= start coordinate')
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
    def attributes(self):
        return self._attributes

    @attributes.setter
    def attributes(self, value):
        """
        Setter function to validate the attributes argument
        Currently not so strict as to require transcript_id as well per GTF2.2 spec.
        Adding this would cause Ensembl GTFs to fail
        """
        if value is not None and (not isinstance(value, dict) or not 'gene_id' in value):
            raise ValueError('attributes must be a dictionary with at least 1 keys: gene_id')
        self._attributes = value

    @property
    def phase(self):
        return self._phase
    
    @phase.setter
    def phase(self, value):
        valid_phase = [0, 1, 2, '.', '0', '1', '2']
        if value not in valid_phase:
            raise ValueError(f'phase must be in {valid_phase}')
        self._phase = str(value)

    def __repr__(self):
        attr_repr = ', '.join(f"{k}={v!r}" for k, v in self.__dict__.items())
        return (f"GtfFeature ({attr_repr})")

    def __str__(self):
        '''
        Return data in gtf format
        '''
        base_str = '\t'.join([str(field) for field in [self.seqname, self.source, self.feature,
            self.start, self.end, self.score, self.strand, self.phase]])
        if self.attributes:
            attribute_str = '; '.join(f'{k} "{v}"' for k, v in self.attributes.items()) + ';'
        else:
            attribute_str = 'MISSING'
        return (f"{base_str}\t{attribute_str}")

    def __eq__(self, other):
        if self.start is None or self.end is None:
            raise ValueError('start and end must be defined to permit comparisons')
        seqnames_match = self.seqname == other.seqname
        strands_match = self.strand == other.strand
        coords_match = self.start == other.start and self.end == other.end
        return seqnames_match and strands_match and coords_match

    def __ne__(self, other):
        return not self == other

    def __gt__(self, other):
        '''
        Define greater than based on seqname and start coordinates to aid sorting.
        Features on reverse strand will be sorted in the order they are encountered
        while traversing the forward strand in the 5' to 3' direction.
        '''
        if self.start is None or self.end is None:
            raise ValueError('start and end must be defined to permit comparisons')
        if self.seqname == other.seqname:
            return self.start > other.start
        else:
            return str(self.seqname) > str(other.seqname)

    def __ge__(self, other):
        return self > other or self == other

    def __lt__(self, other):
        '''
        Define less than based on seqname and start coordinates to aid sorting.
        Features on reverse strand will be sorted in the order they are encountered
        while traversing the forward strand in the 5' to 3' direction.
        '''
        if self.start is None or self.end is None:
            raise ValueError('start and end must be defined to permit comparisons')
        if self.seqname == other.seqname:
            return self.start < other.start
        else:
            return str(self.seqname) < str(other.seqname)

    def __le__(self, other):
        return self > other or self == other

    def contains(self, position):
        """
        Check if position is contained within the feature coordinates, ignoring
        seqname and strand
        """
        return self.start <= position <= self.end

    def overlap(self, coords):
        """
        Report overlapping interval with a set of coordinates, ignoring seqname
        and strand
        """
        overlap_coords = (max(self.start, coords[0]), min(self.end, coords[1]))
        # if there is no overlap, the result will have the larger number in the
        # first index
        if overlap_coords[0] > overlap_coords[1]:
            return None
        else:
            return overlap_coords

    def update(self, overwrite=False, **kwargs):
        """
        Update values for any attributes named as arguments to this function.
        If the attribute is currently set to None it will always be updated.
        If the attribute already has a value, then it will only be overwritten
        if overwrite is set to True.
        """
        # get the values passed by user
        updated_attributes = {k:v for k,v in kwargs.items() if v is not None
            and k not in ['self', 'overwrite']}

        # if start or end is passed, attempt to cast to int as convenience for user
        if 'start' in updated_attributes:
            updated_attributes['start'] = int(updated_attributes['start'])
        if 'end' in updated_attributes:
            updated_attributes['end'] = int(updated_attributes['end'])

        # when trying to update start and end simultaneously, it will fail validation
        # if the new range conflicts with the one that is being changed; this can be
        # prevented by ensuring that 'start' updates first
        if overwrite and ('start' in updated_attributes) and ('end' in updated_attributes):
            # only change if new range is valid; blocks partial changes when validation fails
            # on end coord
            if updated_attributes['start'] > updated_attributes['end']:
                raise ValueError('start coordinate must be <= end coordinate')
            else:
                new_start = updated_attributes.pop('start')
                self.start = new_start

        for attribute, value in updated_attributes.items():
            if value is None:
                continue
            if overwrite or getattr(self, attribute) is None:
                setattr(self, attribute, value)
            elif getattr(self, attribute) != value:
                logger.warning(
                    f'Ignored request to change {attribute} to {value}: '
                    f'{attribute} already set to {getattr(self, attribute)} '
                    f'and overwrite=False')


class Exon(GtfFeature):
    """
    Class representation of an exon feature from a GTF file.

    Extends GtfFeature with the following attributes:
      - exon_features: a dictionary of annotations (i.e., CDS, start/stop
        codons, UTRs)
      - start_phase: number of bases at start before start of next complete
        codon
      - end_phase: number of bases left at end to contribute to next codon

    Extends GtfFeature with the following class methods:
      - _get_feature: intended as an internal function to query for coordinates
        of the specified feature type if contained within exon; wrapped by the
        more specific methods below
      - cds: coordinates of CDS portion of exon (or None)
      - five_prime_utr: coordinates of 5' UTR portion of exon (or None)
      - three_prime_utr: coordinates of 3' UTR portion of exon (or None)
      - start_codon: coordinates of start_codon if exon contains one (or None)
      - stop_codon: coordinates of stop_codon if exon contains one (or None)
    """

    def __init__(self, feature = 'exon', **kwargs):
        super().__init__(feature = feature, **kwargs)
        self.start_phase = None
        self.end_phase = None
        self.exon_features = defaultdict(GtfFeature)

    ##### for possible future use ####
    def _get_feature(self, feature):
        for feature in self.exon_features:
            if feature.feature == feature:
                return (feature.start, feature.end)
            else:
                return None

    def cds(self):
        return self._get_feature('CDS')

    def five_prime_utr(self):
        return self._get_feature('five_prime_utr')

    def three_prime_utr(self):
        return self._get_feature('three_prime_utr')

    def start_codon(self):
        return self._get_feature('start_codon')

    def stop_codon(self):
        return self.get_feature('stop_codon')
    ########################################


class Transcript(GtfFeature):
    """
    Class representation of a transcript feature from a GTF file

    junction_style argument is either `intron` or `exon` and controls how
    junction coordinates are calculated -- `intron` is flanking coords of the
    intron and `exon` uses the last exon position of the donor and acceptor
    """

    def __init__(self, junction_style='intron', **kwargs):
        super().__init__(**kwargs)
        self.junction_style = junction_style
        self._exons = dict()
        self.annotations = {
            'CDS': [],
            'five_prime_utr': [],
            'three_prime_utr': [],
            'start_codon': [],
            'stop_codon': []
        }
        # self.CDS = [] # uppercase because that's how it appears in GTF
        # self.five_prime_utr = []
        # self.three_prime_utr = []
        # self.start_codon = []
        # self.stop_codon = []

    @property
    def junction_style(self):
        return self._junction_style

    @junction_style.setter
    def junction_style(self, value):
        valid_values = ['intron', 'exon']
        if value not in valid_values:
            raise ValueError(f'junction_style must be in {valid_values}')
        self._junction_style = value

    @property
    def exons(self):
        return self._exons

    @exons.setter
    def exons(self, value):
        if not isinstance(value, dict):
            raise ValueError('exons must be type dict')
        for exon in value.items():
            if not isinstance(exon, Exon):
                raise ValueError('every exon in exons must be of type Exon')
        self._exons = value

    def __repr__(self):
        attr_repr = ', '.join(f"{k}={v!r}" for k, v in self.__dict__.items()
            if k not in ['exons', 'junctions', 'annotations'])
        exons_repr = f'exons: {len(self.exons)}'
        # junctions_repr = f'junctions: {len(self.junctions)}'
        annotations_repr = ', '.join("{}={!r}".format(k, len(v)) for k, v in self.annotations.items())
        return (f"{self.__class__.__name__} ({attr_repr}, {exons_repr}, annotations: {{{annotations_repr}}})")

    def add_exon(self, exon):
        """Add a single Exon object to exons dictionary"""
        self.exons[(exon.start, exon.end)] = exon

    def add_annotation(self, annotation):
        self.annotations[annotation.feature].append(annotation)

    def junctions(self):
        """
        Return a dict of all splice junctions. Keys are defined as a tuple of the seqname,
        last base coordinate of the donor exon and the first base coordinate
        of the acceptor exon. This tuple is used as the dictionary key. The value
        is a tuple of the two exons and the strand. The order of the 
        value tuple is always (donor, acceptor, strand). In other words, it will
        seem reversed if the junction is annotated to the minus/reverse strand.
        If stand is unspecified, the order of the tuple is unchanged, but the
        ambiguousness of donor/acceptor will be apparent from the strand.
        """
        junctions = dict()
        sorted_exons = sorted(self.exons.items(), key = lambda x: x[1])
        for exon_idx in range(len(sorted_exons) - 1):
            key_left, exon_left = sorted_exons[exon_idx]
            key_right, exon_right = sorted_exons[exon_idx + 1]
            
            # as an extra safety check, make sure exon and transcript strands match
            if not (self.strand == exon_left.strand) or not (self.strand 
                    == exon_right.strand):
                raise ValueError(
                    f'Strand Mismatch: transcript = {self.strand}, '
                    f'exons = ({exon_left.strand}, {exon_right.strand})')
            if self.junction_style == 'intron':
                key = (self.seqname, exon_left.end + 1, 
                    exon_right.start - 1)
            else:
                key = (self.seqname, exon_left.end, 
                    exon_right.start)

            # arrange values so that donor always comes before acceptor regardless
            # of strand
            if self.strand == '-':
                junctions[key] = (exon_right, exon_left, self.strand)
            else:
                junctions[key] = (exon_left, exon_right, self.strand)
        return junctions


class Gene(GtfFeature):
    """
    Class representation of a gene feature from a GTF file including its
    transcripts
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.transcripts = []

    def add_transcript(self, transcript):
        self.transcripts.append(transcript)


class Transcriptome:
    """
    Class representation of a transcriptome model parsed from a GTF file
    The bed_coords argument will subtract one from the coordinates used for keys
    in the junctions, donors, and acceptors dictionaries to facilitate lookups
    from bed intersections.
    """

    def __init__(self, gtf, species=None, source=None, version=None,
            junction_style='intron'):
        self.genes = defaultdict(Gene)
        self.transcripts = defaultdict(Transcript)
        self.donors = defaultdict(list)
        self.acceptors = defaultdict(list)
        self.junctions = defaultdict(list)
        self.species = species
        self.source = source
        self.version = version
        self.junction_style = junction_style
        self._parse_gtf(gtf)

    @property
    def junction_style(self):
        return self._junction_style

    @junction_style.setter
    def junction_style(self, value):
        valid_values = ['intron', 'exon']
        if value not in valid_values:
            raise ValueError(f'junction_style must be in {valid_values}')
        self._junction_style = value

    def _parse_gtf(self, gtf):
        """
        Parse lines of a GTF and build Transcriptome object from it.
        """
        lines_read = 0
        with open(gtf, 'r') as f:
            logger.info('GTF parsing started')
            for line in dropwhile(lambda s: s.startswith('#'), f):
                lines_read += 1
                parsed = parse_gtf_row(line)
                if parsed['args']['feature'] == 'selenocysteine':
                    continue
                elif parsed['args']['feature'] == 'transcript':
                    self.transcripts[parsed['attributes']['transcript_id']].update(
                        **parsed['args'], attributes = parsed['attributes'],
                        gene_id = parsed['attributes']['gene_id'],
                        transcript_id = parsed['attributes']['transcript_id'])
                    self.genes[parsed['attributes']['gene_id']].add_transcript(
                        parsed['attributes']['transcript_id'])
                elif parsed['args']['feature'] == 'exon':
                    '''
                    Update transcripts dictionary with a placeholder transcript
                    to allow the exon to be attached; the rest of the transcript
                    info will be or already has been grabbed by the transcript
                    feature in the gtf
                    '''
                    exon = Exon(**parsed['args'],
                        attributes = parsed['attributes'],
                        gene_id = parsed['attributes']['gene_id'],
                        transcript_id = parsed['attributes']['transcript_id'])
                    self.transcripts[parsed['attributes']['transcript_id']].update(
                        transcript_id = parsed['attributes']['transcript_id'])
                    self.transcripts[parsed['attributes']['transcript_id']].add_exon(exon)

                    '''
                    For all of the more detailed transcript feature annotations, add them
                    to the proper list in the proper transcript object (which may or may not
                    have been created yet) for further processing. They cannot be processed
                    as they are encountered because they may impact exons that have not yet
                    been encountered.
                    '''
                elif parsed['args']['feature'] in ['CDS', 'five_prime_utr', 'three_prime_utr', 'start_codon', 'stop_codon']:
                    new_feature = GtfFeature(**parsed['args'],
                        attributes = parsed['attributes'],
                        gene_id = parsed['attributes']['gene_id'],
                        transcript_id = parsed['attributes']['transcript_id'])
                    self.transcripts[
                        parsed['attributes']['transcript_id']].add_annotation(
                        new_feature)

                elif parsed['args']['feature'] == 'gene':
                    self.genes[parsed['attributes']['gene_id']].update(
                        **parsed['args'], attributes = parsed['attributes'],
                        gene_id = parsed['attributes']['gene_id'])

                if lines_read > 0 and lines_read % 250000 == 0:
                    logger.info(f'{lines_read} GTF lines read')

            logger.info(f'GTF parsing complete -- {lines_read} lines read')

            logger.info(
                f'Recording known splice junctions: '
                f'junction_style={self.junction_style}')
            for transcript in self.transcripts.values():
                transcript.junction_style = self.junction_style
                junctions = transcript.junctions()
                for junction, details in junctions.items():
                    donor, acceptor, strand = details
                    self.junctions[junction].append(transcript.transcript_id)
                    
                    '''
                    Get donor and acceptor coordinates and put those into their
                    own dictionaries for direct lookups.
                    The expected junction format is (donor, acceptor, strand).
                    Expand range by 1 to include final exon base on each side.
                    Maintain donor/acceptor order such that reverse strand junctions
                    will have coords in decreasing order. If strand is unspecified,
                    treat it as though it was positive.
                    '''
                    if strand == '-':
                        donor_pos = (donor.start, donor.start - 1)
                        acceptor_pos = (acceptor.end + 1, acceptor.end)
                    else:
                        donor_pos = (donor.end, donor.end + 1)
                        acceptor_pos = (acceptor.start - 1, acceptor.start)
                    self.donors[(junction[0], *donor_pos)].append(transcript.transcript_id)
                    self.acceptors[(junction[0], *acceptor_pos)].append(transcript.transcript_id)
            logger.info(
                f'Recorded {len(self.junctions)} known splice junctions: '
                f'{len(self.donors)} donors, {len(self.acceptors)} acceptors')


###### testing code #########
# if __name__ == '__main__':
#     test_line = '1\tensembl\ttranscript\t3102016\t3102125\t.\t+\t.\tgene_id "ENSMUSG00000064842"; gene_version "1"; transcript_id "ENSMUST00000082908"; transcript_version "1"; gene_name "Gm26206"; gene_source "ensembl"; gene_biotype "snRNA"; transcript_name "Gm26206-201"; transcript_source "ensembl"; transcript_biotype "snRNA"; tag "basic"; transcript_support_level "NA";'

#     gtf = '/home/moccir/Share/genomes/annotations/ensembl/release-98/mus_musculus/Mus_musculus.GRCm38.98.gtf'

#     transcriptome = Transcriptome(gtf)
#     mm10 = Transcriptome(gtf, junction_style = 'exon')

#     counts = []
#     for v in mm10.junctions.values():
#         counts.append(len(v))
#     Counter(counts)

#     for k, v in mm10.junctions.items():
#         if len(v) > 42:
#             print(f'{k}: {v}')

#     all_junctions = []
#     for transcript in mm10.transcripts.values():
#         all_junctions.extend(transcript.junctions())
#     len(all_junctions)

#     test_transcript = mm10.transcripts['ENSMUST00000081333']
#     bad_transcript = mm10.transcripts['ENSMUST00000099981']


#     for x in bad_transcript.annotations['CDS']:
#         print(f'{(x.start, x.end)}')
#         for y in bad_transcript.exons.values():
#             overlaps = y.overlap((x.start, x.end))
#             if overlaps:
#                 print(f'{y.attributes["exon_id"]}: {overlaps}')

#     for x in test_transcript.annotations['start_codon']:
#         print(f'{(x.start, x.end)}')
#         for y in test_transcript.exons.values():
#             overlaps = y.overlap((x.start, x.end))
#             print(f'{y.attributes["exon_id"]}: {overlaps}')


    # def time_me():
    #     for x in bad_transcript.annotations['CDS']:
    #         for y in bad_transcript.exons.values():
    #             overlaps = y.overlap((x.start, x.end))
    # %timeit time_me()