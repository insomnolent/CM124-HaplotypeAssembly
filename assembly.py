import os
os.chdir('/Users/christinesun/PycharmProjects/Final Project 2/data')
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
    newSet = []
    for line in set:
        newSet.append(line.replace('-', ''))
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
            res = removeDash(result)
            return res
    res = removeDash(result)
    return res


# filter out reads that don't fit into haplotype 1
def filter(matrix, ref):
    temp = []
    for line in matrix:
        for i in range(0, len(line)):
            if line[i] != ref[i]:
                temp.append(line)
                break
    return temp

# compare other reads to haplotype 1
def compareReads(set, ref):
    temp = ''
    for line in set:
        for i in range(0, len(line)):
            if line[i] != ref[i]:
                temp = line
    return temp

read_matrix = read_input('small_no_error_training_reads.txt')
read_matrix = removedupe(read_matrix)
# set size of read_matrix to new size
lines = len(read_matrix)
# check how long the haplotype should be later
length = len(read_matrix[0])

# set beginning of hap1 and hap2
subset1 = takeSubset(read_matrix, 0)
print subset1
hap1 = largest(subset1)
print 'hap1', hap1
subset2 = filter(subset1, hap1)
print subset2
hap2 = largest(subset2)
print 'hap2', hap2

print 'next index [1]'

for i in range(1, 2):
    readset1 = takeSubset(read_matrix, i)
    print readset1
    read1 = largest(readset1)
    print 'read1', read1

    readset2 = filter(readset1, read1)
    print readset2
    read2 = largest(readset2)
    print 'read2', read2

    # if read1 overlaps with haplotype1, then read2 overlaps with haplotype2
    # doesn't work right now but will fix later
    if read1[0] == hap1[i]:
        count = i
        for j in range(0, len(read1)):
            if read1[j] != hap1[count]:
                print 'no overlap', read1[j]
                hap1 += read1[j]
                print 'new hap1', hap1
        for k in range(0, len(read2)):
            if read2[k] != hap2[count]:
                hap2 += read2[k]
    # if read1 doesn't overlap with hap1 (and so overlaps with hap2)
    else:
        print 'if read1 doesnt overlap'
        count = i
        for j in range(0, len(read1)):
            if read1[j] != hap2[count]:
                print 'no overlap', read1[j]
                hap2 += read1[j]
                print 'new hap1', hap1
        for k in range(0, len(read2)):
            if read2[k] != hap1[count]:
                hap2 += read2[k]

    print 'for testing'
    print hap1
    print hap2

trackLine = trackIndex = 0
hap1 = hap2 = ''

# take a subset (all reads that start at index 0) of the whole matrix
# determine hap1 and hap2 in that matrix
# check with next subset (so all reads that start at next index)
# determine largest hap1 and hap2 in that set
# figure out which overlaps with which haplotype
# continue building each haplotype
# continue reading in line by line

