import logging
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO

from advntr.distance import hamming
from advntr.hierarchical_clustering import hierarchical_clustering


class PacBioHaplotyper:
    """Find error corrected VNTR section in two haplotypes from PacBio reads and genotype the number of repeats"""

    def __init__(self, reads):
        logging.debug('Number of reads for finding haplotypes: %s' % len(reads))
        self.reads = [read.upper() for read in reads]

    def get_error_corrected_haplotypes(self, number_of_clusters=2):
        if len(self.reads) < 2:
            logging.info('There is only one spanning read. Returning it as the only available haplotype')
            return self.reads
        haplotypes = []
        clusters = self.get_read_clusters(number_of_clusters)
        for i, cluster in enumerate(clusters):
            logging.debug('Cluster[%s] size: %s' % (i, len(cluster)))
        smaller_cluster_size = min([len(cluster) for cluster in clusters])
        larger_cluster_size = max([len(cluster) for cluster in clusters])
        homozygous = larger_cluster_size >= 7 * smaller_cluster_size
        for cluster in clusters:
            if len(cluster) == smaller_cluster_size and homozygous:
                logging.info('Both haplotypes are similar (homozygous copy number)')
                continue
            if len(cluster) < 2:
                logging.debug('Cluster has only one sequence. Skipping multiple alignment')
                haplotypes.append(cluster[0])
                continue
            muscle_cline = MuscleCommandline('muscle', clwstrict=True)
            data = '\n'.join(['>%s\n' % str(i) + cluster[i] for i in range(len(cluster))])
            stdout, stderr = muscle_cline(stdin=data)
            alignment = AlignIO.read(StringIO(stdout), "clustal")
            aligned_reads = [str(aligned.seq) for aligned in alignment]
            seq = self.get_consensus_sequence_from_multiple_alignment(aligned_reads)
            haplotypes.append(seq)
        if len(haplotypes) < 2:
            haplotypes.append(haplotypes[0])
        return haplotypes

    @staticmethod
    def get_consensus_sequence_from_multiple_alignment(aligned_reads):
        """
        It finds the most frequent element in each column using counting sort,
        then it add the most frequent element to the result if it is not a gap element.
        If the frequency of the gap and one base pair is equal, the base pair will be selected.
        """
        seq = ''
        for i in range(len(aligned_reads[0])):
            bins = {}
            for row in aligned_reads:
                if row[i] in bins.keys():
                    bins[row[i]] += 1
                else:
                    bins[row[i]] = 0
            sorted_frequencies = sorted(bins.items(), key=lambda x: x[1])
            most_frequent_element = sorted_frequencies[-1][0]
            if most_frequent_element != '-':
                seq += most_frequent_element

        return seq

    def get_read_clusters(self, number_of_clusters=2):
        """Cluster reads to two group based on informative base pairs to separate the reads of each haplotype"""
        muscle_cline = MuscleCommandline('muscle', clwstrict=True)
        data = '\n'.join(['>%s\n' % str(i) + self.reads[i] for i in range(len(self.reads))])
        stdout, stderr = muscle_cline(stdin=data)
        alignment = AlignIO.read(StringIO(stdout), "clustal")
        aligned_reads = [str(aligned.seq) for aligned in alignment]
        aligned_read_ids = [str(aligned.id) for aligned in alignment]

        seqs = self.get_informative_columns(aligned_reads)
        distance_matrix = []
        for i in range(len(seqs)):
            distance_matrix.append([])
            for seq in seqs:
                distance_matrix[i].append(hamming(seq, seqs[i]))
        clusters = hierarchical_clustering(number_of_clusters, distance_matrix)

        result = [[self.reads[int(aligned_read_ids[i])] for i in cluster] for cluster in clusters]
        return result

    @staticmethod
    def get_informative_columns(aligned_reads):
        result = ['' for _ in aligned_reads]
        number_of_columns = len(aligned_reads[0])-1
        for col in range(number_of_columns):
            bins = {}
            for row in aligned_reads:
                if row[col] in bins.keys():
                    bins[row[col]] += 1
                else:
                    bins[row[col]] = 0
            sorted_frequencies = sorted(bins.items(), key=lambda x: x[1])
            highest_frequency = sorted_frequencies[-1][1]
            if highest_frequency <= len(aligned_reads) * 0.7:
                for i in range(len(aligned_reads)):
                    result[i] += aligned_reads[i][col]
        return result
