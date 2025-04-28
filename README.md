# Annofilter-junctions
This terribly named program annotates splice junctions, including novel splice junctions, and optionally filters them based on various features including gene_name, gene_id, gene_biotype, chromsome/seqname, and various flags (see below) that are applied during the annotation process like chimeric splice junctions. The annotations are taken from a GTF file input. The junctions are read in from the output of Portcullis and a STAR SJ.out.tab file. Outputs include a tsv file of annotated splice junctions (including novel junctions) as well as a new STAR SJ.out.tab file that contains only the junctions passing all filters. This cleaned junctions file can then be used in the second pass of STAR 2-pass mode. The utility is the second pass alignment is improved when false positive novel junctions from the first pass are removed. Additionally, the second pass can be sped up by filtering out genes with very large numbers of novel splice junctions, particularly if these are not of biological interest.

## Getting started
python annofilter_junctions.py --help

## Flags
  - strand_ambiguous
  - strand_inferred
  - inferred_strand_disagrees
  - chimeric
  - unannotated_donor
  - unannotated_acceptor
  - donor_rss
  - acceptor_rss
  - junction_not_analyzed
  - no_primary_gene
  - multiple_primary_genes
  - multiple_primary_biotypes

## Junction category classifications
  - DA: this is a known annotated splice site in the GTF
  - NDA: both the donor and the acceptor are known in the GTF, but their use together as a splice junction is novel
  - D: only the donor is known
  - A: only the acceptor is known
  - N: novel splice junction with neither the donor or acceptor having been seen previously in annotated splice junctions from the GTF file
