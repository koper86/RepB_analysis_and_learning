#!/usr/bin/env python

import os
import re
import csv

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO, Entrez


def retrieve_blast_first_hit(seq_record):
    # Tworzymy plik, zeby za kazdym razem nie strzelac do blasta
    blast_filename = seq_record.id + '_blast.xml'
    if not os.path.isfile(blast_filename):
        result_handle = NCBIWWW.qblast('tblastn', 'nr', seq_record.seq)
        with open(blast_filename, 'w') as out_handle:
            out_handle.write(result_handle.read())
            out_handle.close()
            result_handle.close()
            print('Saved ' + blast_filename)
    # Parsujemy wynik wyszukiwania w blascie na obiekt blast record
    result_handle = open(blast_filename)
    blast_records = NCBIXML.parse(result_handle)
    blast_record = next(blast_records)
    # wyciagniecie pierwszego przyrownania z obiektu blast record
    try:
        first_alignment = blast_record.alignments[0]
    except IndexError:
        pass
    return first_alignment


def designate_region_genome_coordinates(subject_srt_coordinate, subject_nd_coordinate, distance):
    # wyznczenie koordynatow do wyciagniecia dlugiego kawalka, w zaleznosci od tego na ktorej nici jest hit
    if subject_srt_coordinate < subject_nd_coordinate:
        region_srt_coordinate = subject_srt_coordinate - distance
        region_nd_coordinate = subject_nd_coordinate + distance
    else:
        region_srt_coordinate = subject_nd_coordinate - distance
        region_nd_coordinate = subject_srt_coordinate + distance
    if region_srt_coordinate < 0:
        region_srt_coordinate = 0
    return region_srt_coordinate, region_nd_coordinate


def retrieve_genbank_rec_by_accession(accession):
    # wyciagniecie odpowiedniego rekordu DNA z genebanku i zapis do pliku
    Entrez.email = 'piotrkoper86@gmail.com'
    genbank_replicon_filename = accession + '.gbk'

    # tworzymy plik, zeby nie ciagnac za kazdym razem z genbanku
    if not os.path.isfile(genbank_replicon_filename):
        net_handle = Entrez.efetch(db='nucleotide', id=hit_accession, rettype='gb', retmode='text')
        out_handle = open(genbank_replicon_filename, "w")
        out_handle.write(net_handle.read())
        out_handle.close()
        net_handle.close()
        print('Saved')
    # sparsowanie rekordu genebanku
    return SeqIO.read(genbank_replicon_filename, 'genbank')


# Wczytujemy plik fasta - przy pomocy parse, zwracany jest iterator
fasta_filename = 'PATRIC_RepBb_global_families_Rhizobium_Agrobacterium_reduced_set.fasta'
fasta_iterator = SeqIO.parse(open(fasta_filename), 'fasta')

for seq_record in fasta_iterator:
    print(seq_record.id)
    print(seq_record.seq)

    blast_first_alignment = retrieve_blast_first_hit(seq_record)

    # wyciagniecie accession  i innych parametrow dla pierwszego przyrowania
    hit_accession = blast_first_alignment.accession
    blast_record_hsps_object = blast_first_alignment.hsps[0]
    subject_start_coordinate = blast_record_hsps_object.sbjct_start
    subject_end_coordinate = blast_record_hsps_object.sbjct_end

    genbank_record = retrieve_genbank_rec_by_accession(hit_accession)
    region_start_coordinate, region_end_coordinate = designate_region_genome_coordinates(subject_start_coordinate,
                                                                                         subject_end_coordinate, 10000)
    # wyciagniecie samej sekwencji i wyciecie regionu po distance
    replicon_sequence = genbank_record.seq
    region_replicon_sequence = replicon_sequence[region_start_coordinate:region_end_coordinate]

    # skompilowanie regexa do parS-ow (wersja chalupnicza)
    p = re.compile('[g,a]tt[a-z]{4}gc[a-z]{4}aa[t,c]', re.IGNORECASE)

    all_motif_occurences = p.findall(str(region_replicon_sequence))

    print(all_motif_occurences)
    with open('Result.csv', 'a', newline='') as csvfile:
        rowwriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        rowwriter.writerow(
            [seq_record.id]
            + [str(seq_record.seq)]
            + [hit_accession + '.gbk' + ' ' + str(subject_start_coordinate) + '..' + str(subject_end_coordinate)]
            + ['  '.join(all_motif_occurences)])
