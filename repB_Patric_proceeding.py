from Bio import SeqIO
from Bio.Seq import Seq
from Bio import motifs
import csv

def hamming_distance(s1,s2):
    result=0
    if len(s1)!=len(s2):
        print("String are not equal")
    else:
        for x,(i,j) in enumerate(zip(s1,s2)):
            if i!=j:
                result+=1
    return result

tsv_file_path = '/home/pkoper/PycharmProjects/RepB_analysis_and_learning/patric_test/RepB_all_patric.tbl'
csv_result_file_path = '/home/pkoper/PycharmProjects/RepB_analysis_and_learning/patric_test/pssm_parS_search_result.csv'


parS_instances = [Seq('GTGGTCAGCTGACCAC'),
                  Seq('GTTGGGAGCTCCCAAC'),
                  Seq('GTTGGGAGCTCCCAAC'),
                  Seq('GTTGGGAGCTCCCAAC'),
                  Seq('GTTGGGAGCTCCCAAC'),
                  Seq('GTTAGCCGCGGCAAAC'),
                  Seq('GTTGTCAGCTGACAAC'),
                  Seq('GTTGTCAGCTGACAAC'),
                  Seq('GTTAGCCGCGGCTAAC'),
                  Seq('GTTAGCCGCGGCTAAC'),
                  Seq('GTTGACTGCAGTCAAC'),
                  Seq('GTTCTCAGCTGAGAAC'),
                  Seq('GTTCTCAGCTGAGAAC'),
                  Seq('GTTGTCCGCGGACAAC'),
                  Seq('GTTGTCCGCGGACAAC'),
                  Seq('GTTGGCCGCGGCCAAC'),
                  Seq('GTTGGCCGCGGCCAAC'),
                  Seq('GTTGACACCGGTCAAC'),
                  Seq('GTTCCCCGCGGGGAAC'),
                  Seq('GTTGTCCGCGGACAAC'),
                  Seq('GTTGTCCGCGGACAAC'),
                  Seq('CATAACAGCTGTTAAG'),
                  Seq('GTTGGGAGCTCCCAAC'),
                  Seq('GTTCTTTGCAAAGAAC'),
                  Seq('GTTGGGAGCTCCCAAC'),
                  Seq('GTTTACTGCGGTAAAC'),
                  Seq('GTTGTCAGCTGACAAT'),
                  Seq('GTTGTCAGCTGACAAT'),
                  Seq('GTGGTCAGCTGACCAC'),
                  Seq('GTTGGCCGCGGCCAAC'),
                  Seq('GTTGGCCGCGGCCAAC'),
                  Seq('GTTGGTCGCGACCAAC'),
                  Seq('GTTGCCAGCCGGAAAC'),
                  Seq('ATTGTCAGCTGACAAT'),
                  Seq('ATTGTCAGCTGACAAT'),
                  Seq('GTTGGCCGCGGCCAAC'),
                  Seq('GTTGACAGCTGTCAAC'),
                  Seq('GTTCACAGCTGTGAAC'),
                  Seq('GTTCACAGCTGTGAAC'),
                  Seq('GTTCACAGCTGTGAAC'),
                  Seq('GTTTCCGGCCGGAAAC'),
                  Seq('GTTGACAGCTGTCAAC'),
                  Seq('GTTGACAGCTGTCAAC'),
                  Seq('CTTGAACGCGTTCAAG'),
                  Seq('CTTGAACGCGTTCAAG'),
                  Seq('CTTGAACGCGTTCAAG'),
                  Seq('GTCCACCGCGGTGGAC'),
                  Seq('GTCCACCGCGGTGGAC'),
                  Seq('GTTCCCCGCGGGGAAC'),
                  Seq('GTTCCCCGCGGGGAAC'),
                  Seq('GTTCCCCGCGGGGAAC'),
                  Seq('GTTCCCCGCGGGGTAC'),
                  Seq('GTTCCCCGCGGGGAAC'),
                  Seq('GTTCCCCGCGGGGAAT'),
                  Seq('GTTGTCCGCGGACAAC'),
                  Seq('GTTGTCCGCGGACAAC'),
                  Seq('ATTGTCAGCTGACAAT'),
                  Seq('ATTGTCAGCTGACAAT')]
# buduję motyw
m = motifs.create(parS_instances)
# Position-Weight Matrix
pwm = m.counts.normalize(pseudocounts=0.5)

# Position-Weight Matrix (odrwócone powtórzenie)
rpwm = pwm.reverse_complement()

# zamiana na Position-Specific Scoring Matrices
pssm = pwm.log_odds()
rpssm = rpwm.log_odds()

# otwieram plik z danymi z Patrica
tsv_file = open(tsv_file_path)
read_tsv = csv.reader(tsv_file, delimiter='\t')

next(read_tsv)
for row in read_tsv:
    id = row[0]
    upstream = Seq(row[3])
    cds_downstream = Seq(row[4])
    aa_sequence = Seq(row[5])
    print(id)
    # szukanie po dokladnych motywach
    # for pos, seq in m.instances.search(upstream):
    #     print("%i %s" % (pos, seq))
    #
    # for pos, seq in m.instances.search(cds_downstream):
    #     print("%i %s" % (pos, seq))

    # szukanie w oparciu o PSSM (ustawiłem próg na 10),
    # trzeba pamiętać, że jest to odwrócone potwórzenie i czasami ten sam motyw pojawia się dwa razy tylko z indexem na minusie
    parS_match = set()
    if len(upstream) >= 16:
        for pos, score in pssm.search(upstream, threshold=15.0):
            print("Position %d: score = %5.3f seq: %s" % (pos, score, upstream[pos:pos + len(m)]))
            parS_match.add(str(upstream[pos:pos + len(m)].upper()))

    if len(cds_downstream) >= 16:
        for pos, score in pssm.search(cds_downstream, threshold=15.0):
            print("Position %d: score = %5.3f seq: %s" % (pos, score, cds_downstream[pos:pos + len(m)]))
            parS_match.add(str(cds_downstream[pos:pos + len(m)].upper()))
    iterator = iter(parS_match)
    if len(parS_match) == 1 or (len(parS_match) == 2 and hamming_distance(next(iterator), next(iterator)) < 2):
        with open(csv_result_file_path, mode='a') as csv_file:
            csv_writer = csv.writer(csv_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            csv_writer.writerow([id] + [parS_match.pop()] + [aa_sequence])


