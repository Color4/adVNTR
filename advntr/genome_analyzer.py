import logging
import os
from uuid import uuid4

from Bio import SeqIO

from advntr.blast_wrapper import get_blast_matched_ids, make_blast_database
from advntr.profiler import time_usage
from advntr.sam_utils import extract_unmapped_reads_to_fasta_file
from advntr.vntr_finder import VNTRFinder


class GenomeAnalyzer:
    def __init__(self, reference_vntrs, target_vntr_ids, working_directory='./'):
        self.reference_vntrs = reference_vntrs
        self.target_vntr_ids = target_vntr_ids
        self.working_dir = working_directory

        self.vntr_finder = {}
        for ref_vntr in self.reference_vntrs:
            if ref_vntr.id in target_vntr_ids:
                self.vntr_finder[ref_vntr.id] = VNTRFinder(ref_vntr)

    @staticmethod
    def print_genotype(vntr_id, copy_numbers, num_of_reads=None, confidence=None):
        print(vntr_id)
        result = ''
        if copy_numbers is not None:
            result = '/'.join([str(cn) for cn in sorted(copy_numbers)])
        else:
            result = 'None'
        if num_of_reads is not None:
            result += '#%s#%s' % (num_of_reads, confidence)
        print(result)

    @time_usage
    def get_filtered_read_ids(self, read_file, illumina=True):
        db_name = 'blast_db__' + os.path.basename(read_file)
        blast_db_name = self.working_dir + db_name
        empty_db = False
        if not os.path.exists(blast_db_name + '.nsq') and not os.path.exists(blast_db_name + '.nal'):
            empty_db = make_blast_database(read_file, blast_db_name)

        import sys
        ev = sys.maxsize
        word_size = '11'
        identity_cutoff = '100'
        if not illumina:
            identity_cutoff = '70'
        min_keywords = 10 if illumina else 1

        vntr_read_ids = {}
        if not empty_db:
            counts = {}
            for vid in self.target_vntr_ids:
                vntr_read_ids[vid] = []
                keywords = self.vntr_finder[vid].get_keywords_for_filtering(illumina)
                query = '@'.join(keywords)
                query = list(keywords)[0]

                search_id = str(uuid4())
                search_results = get_blast_matched_ids(query, blast_db_name, max_seq='100000', word_size=word_size,
                                                       evalue=ev, search_id=search_id, identity_cutoff=identity_cutoff)
                counts[vid] = {}
                for read_id in search_results:
                    counts[vid][read_id] = counts[vid].get(read_id, 0) + 1

                if vid in counts.keys():
                    blast_ids = set([read_id for read_id, count in counts[vid].items() if count >= min_keywords])
                    from random import random
                    sampling_rate = min(2000.0 / len(blast_ids), 1) if len(blast_ids) > 0 else 1
                    blast_ids = set([e for e in blast_ids if random() <= sampling_rate])
                else:
                    blast_ids = []
                logging.info('blast selected %s reads for %s' % (len(blast_ids), vid))
                vntr_read_ids[vid] = blast_ids

        return vntr_read_ids

    @time_usage
    def get_vntr_filtered_reads_map(self, read_file, illumina=True):
        """
        :param read_file: Fasta file containing unmapped reads of sample
        :param illumina: If this parameter is set, it shows that sample contains short reads
        :return: All filtered reads, and a map from VNTR id to read id
        """
        vntr_read_ids = self.get_filtered_read_ids(read_file, illumina)

        reads = []
        if len(vntr_read_ids.values()) > 0:
            unmapped_reads = SeqIO.parse(read_file, 'fasta')
            for read in unmapped_reads:
                for vntr_id in vntr_read_ids.keys():
                    if read.id in vntr_read_ids[vntr_id]:
                        reads.append(read)
                        break
        return reads, vntr_read_ids

    def find_repeat_counts_from_pacbio_alignment_file(self, alignment_file):
        unmapped_reads_file = extract_unmapped_reads_to_fasta_file(alignment_file, self.working_dir)
        filtered_reads, vntr_reads_ids = self.get_vntr_filtered_reads_map(unmapped_reads_file, False)

        for vid in self.target_vntr_ids:
            reads = [read for read in filtered_reads if read.id in vntr_reads_ids[vid]]
            copy_numbers, num_of_reads, confidence = self.vntr_finder[vid].find_repeat_count_from_pacbio_alignment_file(alignment_file, reads)
            self.print_genotype(vid, copy_numbers, num_of_reads, confidence)

    def find_repeat_counts_from_pacbio_reads(self, read_file, naive=False):
        filtered_reads, vntr_reads_ids = self.get_vntr_filtered_reads_map(read_file, False)
        for vid in self.target_vntr_ids:
            unmapped_reads = [read for read in filtered_reads if read.id in vntr_reads_ids[vid]]
            copy_numbers = self.vntr_finder[vid].find_repeat_count_from_pacbio_reads(unmapped_reads, naive)
            self.print_genotype(vid, copy_numbers)

    def find_frameshift_from_alignment_file(self, alignment_file):
        for vid in self.target_vntr_ids:
            result = self.vntr_finder[vid].find_frameshift_from_alignment_file(alignment_file, [])
            print(vid)
            print(result)

    def find_repeat_counts_from_alignment_file(self, alignment_file, average_coverage):
        unmapped_reads_file = extract_unmapped_reads_to_fasta_file(alignment_file, self.working_dir)
        filtered_reads, vntr_reads_ids = self.get_vntr_filtered_reads_map(unmapped_reads_file)
        for vid in self.target_vntr_ids:
            unmapped_reads = [read for read in filtered_reads if read.id in vntr_reads_ids[vid]]
            copy_number, num_of_reads, confidence = self.vntr_finder[vid].find_repeat_count_from_alignment_file(alignment_file, unmapped_reads,
                                                                                      average_coverage)
            self.print_genotype(vid, copy_number, num_of_reads, confidence)

    def find_repeat_counts_from_short_reads(self, read_file):
        for vid in self.target_vntr_ids:
            copy_number = self.vntr_finder[vid].find_repeat_count_from_short_reads(read_file)
            self.print_genotype(vid, copy_number)
