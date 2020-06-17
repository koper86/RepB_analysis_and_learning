import csv
import sys
import Levenshtein as lev

from Bio.Seq import Seq

number_of_mismatches = sys.argv[1]
print(number_of_mismatches)

with open('Result.csv', newline='') as csvfile:
    rowreader = csv.reader(csvfile, delimiter=',', quotechar='|')
    for row in rowreader:
        hipo_parS_sequence = row[3]
        print(hipo_parS_sequence)
        hipo_parS_sequences = hipo_parS_sequence.split(' ')
        flag = False
        inner_sequence = []
        for parS in hipo_parS_sequences:
            if len(parS) != 0:
                # wyciagam srodek z palindromu (N-ki)
                inner_parS = Seq(parS[3:7] + parS[9:13])
                rc_inner_parS = inner_parS.reverse_complement()
                distance = lev.distance(str(inner_parS).lower(), str(rc_inner_parS).lower())
                if distance <= int(number_of_mismatches):
                    flag = True
                    inner_sequence.append(str(inner_parS) + '..' + str(rc_inner_parS))
        if flag:
            with open('Result_filtered_' + number_of_mismatches + '.csv', 'a', newline='') as csvfile:
                rowwriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
                rowwriter.writerow(
                    [row[0]]
                    + [row[1]]
                    + [row[2]]
                    + [row[3]]
                    + [' '.join(inner_sequence)])
