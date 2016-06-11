import os
import sys
from collections import defaultdict

os.chdir('/Users/christinesun/GitHub/CM124-HaplotypeAssembly/Final_2')
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


# removes dashes before and after each string in the subset
def removeDash(set):
    newSet = [i for i in set if i != '-']
    return newSet


# return largest read in that subset
# for data sets where read lengths aren't all the same
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


# filter out reads that don't fit into haplotype
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


# finds reads from subset that are similar (only 1 error) to read
def findSimilar(subset):
    # check each col of each line of subset and see how many cols match
    # for some reason import Counter doesn't work
    str = ""
    freq = []
    for i in range(0, len(subset)):
        str += subset[i]
        str += " "
    return freq


def gettopreads(subset):
    return top2(get_frequency(subset))


def get_frequency(subset):
    freq = defaultdict(int)
    for c in subset:
        freq[c] += 1
    return freq


def top2(frequency):
    array = [(k,v) for k,v in frequency.items()]
    sorted_array = sorted(array, key=lambda x: -x[1])
    if len(sorted_array) == 1:
        return [k for k, v in sorted_array[0:1]]
    elif len(sorted_array) == 2:
        return [k for k, v in sorted_array[0:2]]
    elif len(sorted_array) == 3:
        return [k for k, v in sorted_array[0:3]]
    elif len(sorted_array) == 4:
        return [k for k, v in sorted_array[0:4]]
    elif len(sorted_array) == 5:
        return [k for k, v in sorted_array[0:5]]
    elif len(sorted_array) == 6:
        return [k for k, v in sorted_array[0:6]]
    else:
        return [k for k, v in sorted_array[0:7]]


def checkInverse(read1, read2):
    track = 0
    check = True
    for c in read1:
        if c == "1":
            if read2[track] != "0":
                check = False
                break
            track += 1
        elif c == "0":
            if read2[track] != "1":
                check = False
                break
            track += 1
    return check

# for reversing a string
# sum([abs(map(int, '10001')[i] - map(int, '0110')[i]) for i in range(len('1001'))])

# change newFiltered to account for some errors
# in both hap1 and hap2 reads
# maybe compare it column by column
def newFiltered(reads, length):
    new = []
    for i in range(0, length):
        readset1 = takeSubset(reads, i)
        topreads = gettopreads(readset1)
        if len(topreads) == 1:
            new.append(topreads[0])
        elif len(topreads) > 1:
            read1 = removeDash(topreads[0])
            read2 = removeDash(topreads[1])

            if checkInverse(read1, read2):
                new.append(topreads[0])
                new.append(topreads[1])
            else:
                # if the most frequent reads aren't inverses of each other
                for j in range(0, len(topreads)-1):
                    read1 = topreads[j]
                    for k in range(j+1, len(topreads)):
                        read2 = topreads[k]
                        if checkInverse(removeDash(read1), removeDash(read2)):
                            new.append(read1)
                            new.append(read2)
                            break
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


# returns inverse of a read
def inverse(read):
    inv = ""
    for c in read:
        if c == "1":
            inv += "0"
        else:
            inv += "1"
    return inv


read_matrix = read_input('./medium/example/input/reads_low_error.txt')
# read_matrix = read_input('./easy_high_error_2_chromosomes_training_reads.txt')

# set size of read_matrix to new size
lines = len(read_matrix)
# check how long the haplotype should be later
length = len(read_matrix[0])

# make new matrix with longest reads from each subset
new_matrix = newFiltered(read_matrix, length)

print 'new_matrix '
for c in new_matrix[0:500]:
    print c

# get rid of empty lines in new_matrix
new_matrix = [i for i in new_matrix if i != '']

# write new file with new_matrix for testing purposes
# test = open('test_subset.txt', 'w')
# for k in range(0, len(new_matrix)):
#     test.write(new_matrix[k])
#     test.write('\n')
# test.close()

end = False
hap1 = removeDash(new_matrix[0])
# if only one read for the first index of new matrix
if new_matrix[1][0] == '-':
    hap2 = inverse(hap1)
else:
    hap2 = removeDash(new_matrix[1])


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
for i in range(0, len(hap1)):
    sys.stdout.write(hap1[i])

print "\n"
for i in range(0, len(hap2)):
    sys.stdout.write(hap2[i])