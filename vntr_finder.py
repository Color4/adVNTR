from Bio import SeqIO
from multiprocessing import Semaphore, Lock, Process
import numpy
import os
import pysam
from random import random
from uuid import uuid4

from blast_wrapper import get_blast_matched_ids, make_blast_database
from coverage_bias import CoverageBiasDetector, CoverageCorrector
from hmm_utils import *
from sam_utils import get_related_reads_and_read_count_in_samfile, extract_unmapped_reads_to_fasta_file
import settings


class VNTRFinder:
    """Find the VNTR structure of a reference VNTR in NGS data of the donor."""

    def __init__(self, reference_vntr):
        self.reference_vntr = reference_vntr
        self.min_repeat_bp_to_add_read = 2
        if len(self.reference_vntr.pattern) < 30:
            self.min_repeat_bp_to_add_read = 2
        self.min_repeat_bp_to_count_repeats = 2

    def build_vntr_matcher_hmm(self, copies, read_length=150):
        patterns = self.reference_vntr.get_repeat_segments() * 100
        flanking_region_size = read_length - 10
        left_flanking_region = self.reference_vntr.left_flanking_region[-flanking_region_size:]
        right_flanking_region = self.reference_vntr.right_flanking_region[:flanking_region_size]

        vntr_matcher = get_suffix_matcher_hmm(left_flanking_region)
        right_flanking_matcher = get_prefix_matcher_hmm(right_flanking_region)
        repeats_matcher = get_variable_number_of_repeats_matcher_hmm(patterns, copies)
        vntr_matcher.concatenate(repeats_matcher)
        vntr_matcher.concatenate(right_flanking_matcher)
        vntr_matcher.bake(merge=None)
        return vntr_matcher

    def get_vntr_matcher_hmm(self, read_length):
        """Try to load trained HMM for this VNTR
        If there was no trained HMM, it will build one and store it for later usage
        """
        copies = int(round(float(read_length) / len(self.reference_vntr.pattern) + 0.5))

        base_name = str(self.reference_vntr.id) + '_' + str(read_length) + '.json'
        stored_hmm_file = settings.TRAINED_HMMS_DIR + base_name
        if settings.USE_TRAINED_HMMS and os.path.isfile(stored_hmm_file):
            model = Model()
            model = model.from_json(stored_hmm_file)
            return model

        vntr_matcher = self.build_vntr_matcher_hmm(copies, read_length)

        json_str = vntr_matcher.to_json()
        with open(stored_hmm_file, 'w') as outfile:
            outfile.write(json_str)
        return vntr_matcher

    def filter_reads_with_keyword_matching(self, working_directory, short_read_file):
        word_size = int(len(self.reference_vntr.pattern)/3)
        if word_size > 11:
            word_size = 11
        word_size = str(word_size)
        blast_ids = set([])

        db_name = working_directory[:-1]
        blast_db_name = working_directory + db_name
        if not os.path.exists(blast_db_name + '.nal'):
            make_blast_database(short_read_file, blast_db_name)
        search_id = str(uuid4()) + str(self.reference_vntr.id)
        for repeat_segment in self.reference_vntr.get_repeat_segments():
            blast_ids |= get_blast_matched_ids(repeat_segment, blast_db_name, max_seq='50000',
                                               evalue=10, word_size=word_size, search_id=search_id)

        print('blast selected ', len(blast_ids), ' reads')
        if len(blast_ids) == len(self.reference_vntr.get_repeat_segments()) * 50 * 1000:
            with open('errors.txt', 'a') as out:
                out.write('maximum number of read selected in filtering for pattern %s\n' % self.reference_vntr.id)
        return blast_ids

    @staticmethod
    def add_hmm_score_to_list(semaphore, list_lock, hmm, read, result_scores):
        logp, vpath = hmm.viterbi(str(read.seq))
        list_lock.acquire()
        result_scores.append(logp)
        list_lock.release()
        semaphore.release()

    def calculate_min_score_to_select_a_read(self, hmm, alignment_file):
        """Calculate the score distribution of false positive reads
        and return score to select the 0.0001 percentile of the distribution
        """
        semaphore = Semaphore(settings.CORES)
        list_lock = Lock()
        process_list = []
        false_scores = []
        read_mode = 'r' if alignment_file.endswith('sam') else 'rb'
        samfile = pysam.AlignmentFile(alignment_file, read_mode)
        for read in samfile.fetch():
            if read.is_unmapped:
                continue
            if random() > settings.SCORE_FINDING_READS_FRACTION:
                continue
            read_start = read.reference_start
            read_end = read.reference_end if read.reference_end else read_start + len(read.seq)
            vntr_start = self.reference_vntr.start_point
            vntr_end = vntr_start + self.reference_vntr.get_length()
            reference_name = read.reference_name
            if not reference_name.startswith('chr'):
                reference_name = 'chr' + reference_name
            if reference_name == self.reference_vntr.chromosome and (
                        vntr_start <= read_start < vntr_end or vntr_start < read_end <= vntr_end):
                continue
            if read.seq.count('N') > 0:
                continue
            semaphore.acquire()
            p = Process(target=VNTRFinder.add_hmm_score_to_list, args=(semaphore, list_lock, hmm, read, false_scores))
            process_list.append(p)
            p.start()
        for p in process_list:
            p.join()

        score = numpy.percentile(false_scores, 100 - 0.0001)
        return score

    def get_min_score_to_select_a_read(self, hmm, alignment_file, read_length):
        """Try to load the minimum score for this VNTR

        If the score is not stored, it will compute the score and write it for this VNTR in precomputed data.
        """
        base_name = str(self.reference_vntr.id) + '_' + str(read_length) + '.scores'
        stored_scores_file = settings.TRAINED_HMMS_DIR + base_name
        if settings.USE_TRAINED_HMMS and os.path.isfile(stored_scores_file):
            with open(stored_scores_file, 'r') as infile:
                frac_score = [(line.split()[0], line.split()[1]) for line in infile.readlines() if line.strip() != '']
                fraction_score_map = {float(reads_fraction): float(score) for reads_fraction, score in frac_score}
            if settings.SCORE_FINDING_READS_FRACTION in fraction_score_map.keys():
                return fraction_score_map[settings.SCORE_FINDING_READS_FRACTION]

        print('Minimum score is not precomputed for vntr id: %s' % self.reference_vntr.id)
        score = self.calculate_min_score_to_select_a_read(hmm, alignment_file)
        print('computed score: ', score)
        with open(stored_scores_file, 'a') as outfile:
            outfile.write('%s %s\n' % (settings.SCORE_FINDING_READS_FRACTION, score))

        return score

    def process_unmapped_read(self, sema, result_lock, read_segment, hmm, filtered_read_ids, min_score_to_count_read,
                              vntr_bp_in_unmapped_reads, selected_reads):
        if read_segment.id in filtered_read_ids and read_segment.seq.count('N') <= 0:
            logp, vpath = hmm.viterbi(str(read_segment.seq))
            rev_logp, rev_vpath = hmm.viterbi(str(read_segment.seq.reverse_complement()))
            if logp < rev_logp:
                logp = rev_logp
                vpath = rev_vpath
            repeat_bps = get_number_of_repeat_bp_matches_in_vpath(vpath)
            if logp > min_score_to_count_read:
                result_lock.acquire()
                if repeat_bps > self.min_repeat_bp_to_count_repeats:
                    vntr_bp_in_unmapped_reads += repeat_bps
                if repeat_bps > self.min_repeat_bp_to_add_read:
                    selected_reads.append(read_segment.id)
                result_lock.release()
        sema.release()

    def find_repeat_count_from_alignment_file(self, alignment_file, working_directory='./'):

        if working_directory == './':
            working_directory = os.path.dirname(alignment_file) + '/'

        unmapped_read_file = extract_unmapped_reads_to_fasta_file(alignment_file, working_directory)
        print('unmapped reads extracted')

        filtered_read_ids = self.filter_reads_with_keyword_matching(working_directory, unmapped_read_file)
        print('unmapped reads filtered')

        hmm = None
        min_score_to_count_read = None
        selected_reads = []
        vntr_bp_in_unmapped_reads = 0

        number_of_reads = 0
        read_length = 150

        semaphore = Semaphore(settings.CORES)
        result_lock = Lock()
        process_list = []

        unmapped_reads = SeqIO.parse(unmapped_read_file, 'fasta')
        for read_segment in unmapped_reads:
            if number_of_reads == 0:
                read_length = len(str(read_segment.seq))
            number_of_reads += 1
            if not hmm:
                hmm = self.get_vntr_matcher_hmm(read_length=read_length)
                min_score_to_count_read = self.get_min_score_to_select_a_read(hmm, alignment_file, read_length)

            semaphore.acquire()
            p = Process(target=self.process_unmapped_read, args=(semaphore, result_lock, read_segment, hmm,
                                                                 filtered_read_ids, min_score_to_count_read,
                                                                 vntr_bp_in_unmapped_reads, selected_reads))
            process_list.append(p)
            p.start()
        for p in process_list:
            p.join()

        print('vntr base pairs in unmapped reads:', vntr_bp_in_unmapped_reads)

        vntr_bp_in_mapped_reads = 0
        vntr_start = self.reference_vntr.start_point
        vntr_end = self.reference_vntr.start_point + self.reference_vntr.get_length()
        region_start = vntr_start - settings.MAX_INSERT_SIZE
        region_end = vntr_end + settings.MAX_INSERT_SIZE
        chromosome = self.reference_vntr.chromosome[3:]
        read_mode = 'r' if alignment_file.endswith('sam') else 'rb'
        samfile = pysam.AlignmentFile(alignment_file, read_mode)
        for read in samfile.fetch(chromosome, region_start, region_end):
            if read.is_unmapped:
                continue
            read_end = read.reference_end if read.reference_end else read.reference_start + len(read.seq)
            if vntr_start <= read.reference_start < vntr_end or vntr_start < read_end <= vntr_end:
                end = min(read_end, vntr_end)
                start = max(read.reference_start, vntr_start)
                vntr_bp_in_mapped_reads += end - start
        print('vntr base pairs in mapped reads:', vntr_bp_in_mapped_reads)

        total_counted_vntr_bp = vntr_bp_in_unmapped_reads + vntr_bp_in_mapped_reads
        pattern_occurrences = total_counted_vntr_bp / float(len(self.reference_vntr.pattern))
        bias_detector = CoverageBiasDetector(alignment_file, self.reference_vntr.chromosome, 'GRCh37')
        coverage_bias_corrector = CoverageCorrector(bias_detector.get_gc_content_coverage_map())

        observed_copy_number = pattern_occurrences / coverage_bias_corrector.get_sequencing_mean_coverage()
        scaled_copy_number = coverage_bias_corrector.get_scaled_coverage(self.reference_vntr, observed_copy_number)
        print('observed copy number and scaled copy number: ', observed_copy_number, scaled_copy_number)
        return scaled_copy_number

    def find_repeat_count_from_short_reads(self, short_read_files, working_directory='./'):
        """
        Map short read sequencing data to human reference genome (hg19) and call find_repeat_count_from_alignment_file
        :param short_read_files: short read sequencing data
        :param working_directory: directory for generating the outputs
        """
        if working_directory == './':
            working_directory = os.path.dirname(short_read_files[0]) + '/'
        alignment_file = ''
        return self.find_repeat_count_from_alignment_file(alignment_file, working_directory)

    def find_accuracy(self, samfile='original_reads/paired_dat.sam'):
        """Find sensitivity and false positive reads for a set of simulated data
        """
        reference_end_pos = self.reference_vntr.start_point + self.reference_vntr.get_length()
        related_reads, read_count = get_related_reads_and_read_count_in_samfile(self.reference_vntr.pattern,
                                                                                self.reference_vntr.start_point,
                                                                                read_file=samfile,
                                                                                pattern_end=reference_end_pos)
        # TODO
        selected_reads = []
        occurrences = 0
        avg_coverage = 1
        true_positives = [read for read in selected_reads if read in related_reads]
        false_positives = [read for read in selected_reads if read not in true_positives]
        false_negatives = [read for read in related_reads if read not in selected_reads]
        # print('TP:', len(true_positives), 'FP:', len(false_positives), 'selected:', len(selected_reads))
        # print('FN:', len(false_negatives))
        sensitivity = float(len(true_positives)) / len(related_reads) if len(related_reads) > 0 else 0
        if sensitivity > 0.9:
            print(sensitivity, len(false_positives))
        if 1 > sensitivity > 0.9 and len(false_negatives) > 0 and len(false_positives) > 0:
            print('sensitivity ', sensitivity, ' FN:', false_negatives[0], ' FP:', false_positives[0])
        with open('FP_and_sensitivity_HMM_read_scoring_method.txt', 'a') as outfile:
            outfile.write('%s\t%s\t%s\t%s\t%s\n' % (
                len(false_positives), sensitivity, self.reference_vntr.id, len(self.reference_vntr.pattern),
                len(true_positives)))
        error = abs(len(self.reference_vntr.get_repeat_segments()) - occurrences / avg_coverage)
        print(error)