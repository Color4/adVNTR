from subprocess import call
import shlex
import os
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import Seq, SeqRecord, SeqIO

from advntr import settings
from advntr.profiler import time_usage


def make_blast_database(fasta_file, db_name):
    make_db_args = '-in %s -parse_seqids -max_file_sz 2GB -dbtype nucl -out %s -title "%s"' % (fasta_file, db_name, db_name)
    make_db_args = shlex.split(make_db_args)
    call(['makeblastdb'] + make_db_args)

    if os.path.exists(db_name + '.nal') or os.path.exists(db_name + '.nsq'):
        empty_db = False
    else:
        empty_db = True
    return empty_db


def make_blast_database_of_multiple_files(fasta_files, db_name):
    for file_name in fasta_files:
        make_db_args = '-in %s -parse_seqids -max_file_size 8GB -dbtype nucl -out %s' % (file_name, db_name)
        make_db_args = shlex.split(make_db_args)
        call(['makeblastdb'] + make_db_args)

    file_names = ' '.join(fasta_files)
    merge_db_args = '-dblist "%s" -dbtype nucl -out %s -title "%s"' % (file_names, db_name, db_name)
    merge_db_args = shlex.split(merge_db_args)

    call(['blastdb_aliastool'] + merge_db_args)


def run_blast_search(query_file, db, result_file, num_threads, word_size, max_seq, evalue, task, identity_cutoff):
    blastn_cline = NcbiblastnCommandline(query=query_file, db=db, outfmt='"6 sseqid"', dust='no', out=result_file,
                                         num_threads=num_threads, word_size=word_size, max_target_seqs=max_seq,
                                         evalue=evalue, task=task, perc_identity=identity_cutoff)
    blastn_cline()
    with open(result_file) as result_input:
        ids = result_input.readlines()
        matched_ids = [seq_id.strip() for seq_id in ids]

    return matched_ids


@time_usage
def get_blast_matched_ids(query, blast_db_name, word_size='5', max_seq='6000', evalue=10.0, search_id='', threads=None,
                          identity_cutoff='0'):
    if not os.path.exists(settings.BLAST_TMP_DIR):
        os.makedirs(settings.BLAST_TMP_DIR)
    query_file = settings.BLAST_TMP_DIR + search_id + '_query.fasta'
    result_file = settings.BLAST_TMP_DIR + search_id + '_blast_result.txt'
    with open(query_file, "w") as output_handle:
        my_rec = SeqRecord.SeqRecord(seq=Seq.Seq(query), id='query', description='')
        SeqIO.write([my_rec], output_handle, 'fasta')

    if len(query) <= 15:
        task = 'blastn-short'
    else:
        task = 'blastn'

    if not threads:
        threads = settings.CORES

    matched_ids = run_blast_search(query_file, blast_db_name, result_file, threads, word_size, max_seq, evalue, task,
                                   identity_cutoff)

    os.remove(result_file) if os.path.exists(result_file) else None
    os.remove(query_file) if os.path.exists(query_file) else None

    return matched_ids
