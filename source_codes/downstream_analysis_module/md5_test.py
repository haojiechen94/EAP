# md5_test.py
# 2023-11-20
# Haojie Chen/Zhijie Guo
'''
md5_test.py md5_reference_file md5_check_file md5_out_file
'''

import sys

print(sys.argv[1:3])
md5_reference = dict()
with open(sys.argv[1], 'r') as md5_reference_file:
    for line in md5_reference_file:
        word = line.strip().split('  ')
        md5_reference[word[1]] = word[0]
md5_reference_file.close()

md5_check = dict()
with open(sys.argv[2], 'r') as md5_check_file:
    for line in md5_check_file:
        word = line.strip().split('  ')
        md5_check[word[1]] = word[0]
md5_check_file.close()

total = 0
oks = 0
nos = 0
notexist = 0
with open(sys.argv[3], 'w') as md5_out_file:
    for i in md5_reference:
        total += 1
        if i in md5_check and md5_check[i] == md5_reference[i]:
            oks += 1
            md5_out_file.write('\t'.join([md5_reference[i], i, 'OK']) + '\n')
        elif i in md5_check and md5_check[i] != md5_reference[i]:
            nos += 1
            md5_out_file.write('\t'.join([md5_reference[i], i, 'NO']) + '\n')
        elif i not in md5_check:
            notexist += 1
            md5_out_file.write('\t'.join([md5_reference[i], i, 'Not exist']) + '\n')
    md5_out_file.write('Total:%i; OK:%i; NO:%i; Not exist: %i'%(total, oks, nos, notexist))
md5_out_file.close()