import os
os.chdir('/Users/christinesun/GitHub/CM124-HaplotypeAssembly/data')
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
            # for now don't remove dashes since want to create new matrix
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
def compareReads(read, hap):
    result = True
    for j in range(0, len(hap)-1):
        if j >= len(read):
            break
        if read[j] != hap[j+1] and j < len(hap):
            result = False
    return result


# assuming part of read1 overlaps with hap, returns part of read1 that's different from hap
def findDiff(read, hap):
    result = ''
    count = 0
    for j in range(0, len(hap)-1):
        if read[count] == hap[j+1]:
            count += 1
    for k in range(count+1, len(read)):
        result += read[k]
    return result

# compare other reads to haplotype 1
# def compareReads(set, ref):
#     temp = ''
#     for line in set:
#         for i in range(0, len(line)):
#             if line[i] != ref[i]:
#                 temp = line
#     return temp

read_matrix = read_input('small_no_error_training_reads.txt')
read_matrix = removedupe(read_matrix)
# write new file with no duplicates for testing purposes
test = open('test_noDupes.txt', 'w')
for k in range(0, len(read_matrix)):
    test.write(read_matrix[k])
    test.write('\n')
test.close()
# set size of read_matrix to new size
lines = len(read_matrix)
# check how long the haplotype should be later
length = len(read_matrix[0])

# make new matrix with longest reads from each subset
new_matrix = []

# set beginning of hap1 and hap2
subset1 = takeSubset(read_matrix, 0)
hap1 = largest(subset1)
print 'hap1', hap1
subset2 = filter(subset1, hap1)
hap2 = largest(subset2)
print 'hap2', hap2
new_matrix.append(hap1)
new_matrix.append(hap2)

for i in range(1, length):
    readset1 = takeSubset(read_matrix, i)
    read1 = largest(readset1)
    print 'read1', read1
    readset2 = filter(readset1, read1)
    read2 = largest(readset2)
    print 'read2', read2

    new_matrix.append(read1)
    new_matrix.append(read2)

    # if read1 overlaps with hap1
    # currently doesn't work for every single case
    # if compareReads(read1, hap1):
    #     hap1 += findDiff(read1, hap1)
    #     hap2 += findDiff(read2, hap2)
    # else:
    #     hap1 += findDiff(read2, hap1)
    #     hap2 += findDiff(read1, hap2)

    # print 'for testing'
    # print 'current hap1', hap1
    # print 'current hap2', hap2

# get rid of empty lines in new_matrix
new_matrix = [i for i in new_matrix if i != '']

test = open('test_subset.txt', 'w')
for k in range(0, len(new_matrix)):
    test.write(new_matrix[k])
    test.write('\n')
test.close()

print new_matrix

# take a subset (all reads that start at index 0) of the whole matrix
# determine hap1 and hap2 in that matrix
# check with next subset (so all reads that start at next index)
# determine largest hap1 and hap2 in that set
# figure out which overlaps with which haplotype
# continue building each haplotype
# continue reading in line by line

