import os
import sys
os.chdir('/Users/christinesun/GitHub/CM124-HaplotypeAssembly/data 3')
# length is how long each haplotype is
# lines is how many read fragments there are
length = lines = 0


def read_input(input):
    data = open(input, 'r')
    matrix = []
    for line in data:
        # to account for the \n in each line
        matrix.append(line.rstrip())
    data.close()
    return matrix


# to remove repetitive read fragments
def removedupe(reads):
    nodupe = []
    for i in reads:
        if i not in nodupe:
            nodupe.append(i)
    return nodupe


# check if read1 is in read2
def checkReads(read1, read2):
    read1 = read1.replace('-', '')
    read2 = read2.replace('-', '')
    check = True
    if len(read1) > len(read2):
        for j in range(0, len(read1)):
            if read1[j] != read2[j]:
                check = False
    else:
        for j in range(0, len(read2)):
            if read1[j] != read2[j]:
                check = False
    return check


# removes dashes before and after each string in the subset
def removeDash(set):
    newSet = [i for i in set if i != '-']
    return newSet


# return largest read in that subset
def largest(set):
    large = ''
    for line in set:
        if len(line) >= len(large):
            large = line
    return large


# take subset of read_matrix of all fragments that start at that index
# returns subset with no dashes in the lines
def takeSubset(matrix, index):
    result = []
    for line in matrix:
        if index == 0:
            if line[index] != '-':
                result.append(line)
        else:
            if line[index] != '-' and line[index-1] == '-':
                result.append(line)
        if line[index] == '-' and len(result) != 0:
            # for now don't remove dashes since you want to create new matrix
            # res = removeDash(result)
            return result
    # res = removeDash(result)
    return result


# filter out reads that don't fit into haplotype 1
def filter(matrix, ref):
    temp = []
    for line in matrix:
        for i in range(0, len(line)):
            if line[i] != ref[i]:
                temp.append(line)
                break
    return temp


# write function that compares overlap of read to hap1 or hap2
# returns true if it overlaps, false if it doesn't
def compareReads(read, hap, index):
    result = True
    for j in range(index, len(hap)):
        # if finish iterating through read
        if j >= len(read) or read[j] == '-':
            break
            # if part of read disagrees with hap
        if read[j] != hap[j] and read[j] != '-':
            result = False
            break
    return result


def newFiltered(reads, len):
    new = []
    for i in range(0, len):
        readset1 = takeSubset(reads, i)
        read1 = largest(readset1)
        readset2 = filter(readset1, read1)
        read2 = largest(readset2)
        new.append(read1)
        new.append(read2)
    return new


# assuming part of read1 overlaps with hap, returns part of read1 that's different from hap
def findDiff(read, hap):
    result = ''
    count = 0
    start = findStart(read)
    for j in range(start, len(hap)):
        if read[j] == hap[j]:
            count += 1
    noDashes = removeDash(read)
    for k in range(start + count, len(noDashes)+start):
        result += read[k]
    # print result
    return result


# finds index where read fragment starts
def findStart(read):
    result = 0
    while read[result] == '-':
        result += 1
    return result


# test if 2nd read is completely overlapped by 1st read
def overlap(first, second):
    result = True
    for i in range(0, len(first)):
        if i < len(second) and first[i] != second[i]:
            result = False
    return result

read_matrix = read_input('easy_low_error_2_chromosomes_training_reads.txt')
read_matrix = removedupe(read_matrix)
# write new file with no duplicates for testing purposes
test = open('test_noDupes.txt', 'w')
for k in range(0, len(read_matrix)):
    test.write(read_matrix[k])
    test.write('\n')
test.close()

print 'new_matrix '
for c in read_matrix:
    print c

# set size of read_matrix to new size
lines = len(read_matrix)
# check how long the haplotype should be later
length = len(read_matrix[0])

# make new matrix with longest reads from each subset
new_matrix = newFiltered(read_matrix, length)

# get rid of empty lines in new_matrix
new_matrix = [i for i in new_matrix if i != '']
# write new file with new_matrix for testing purposes
test = open('test_subset.txt', 'w')
for k in range(0, len(new_matrix)):
    test.write(new_matrix[k])
    test.write('\n')
test.close()


# make smaller subset with more filtered reads
# for i in range(0, len(new_matrix)-1):
#     line1 = new_matrix[i]
#     for k in range(1, len(new_matrix)):
#         line2 = new_matrix[k]
#         d = difflib.SequenceMatcher(None, line1, line2)
#         match = max(d.get_matching_blocks(), key=lambda x: x[2])
#
#         i, j, k = match
#         # if k is zero, then no overlap
#         if k != 0:
#             diff = d.b[j+k-1:len(line2)]
shorter_matrix = []
end = False
hap1 = removeDash(new_matrix[0])
hap2 = removeDash(new_matrix[1])


# create new matrix with less overlapped reads
copy = []
i = 0
while i < len(new_matrix)-1:
    line1 = new_matrix[i]
    line2 = new_matrix[i+1]
    if overlap(removeDash(line1), removeDash(line2)):
        copy.append(line1)
        # i += 1
    # elif overlap(removeDash(line1), removeDash(line1)):
    #    copy.append(line1)
    else:
        copy.append(line1)
        copy.append(line2)
    i += 2


# list comprehensions
# somelist[:] = [x for x in somelist if not determine(x)]

for i in range(2, len(new_matrix)):
    read1 = new_matrix[i]
    # print read1
    ind = findStart(read1)
    compare1 = compareReads(read1, hap1, ind)
    compare2 = compareReads(read1, hap2, ind)
    if compare1:
        hap1 += findDiff(read1, hap1)
        # print 'hap1 ', hap1
    elif compare2:
        hap2 += findDiff(read1, hap2)
        # print 'hap2 ', hap2

# for testing purposes
print 'hap1'
for i in range(0, len(hap1)):
    sys.stdout.write(hap1[i])

print '\nhap2'
for i in range(0, len(hap2)):
    sys.stdout.write(hap2[i])