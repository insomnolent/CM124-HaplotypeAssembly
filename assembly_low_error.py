import os
import sys
from collections import defaultdict

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
    else:
        return [k for k, v in sorted_array[0:2]]


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


def percentOverlap(hap1, hap2):
    percent = 0
    track = 0
    for c in hap1:
        if c == hap2[track]:
            percent += 1
        track += 1
    return float(percent) / float(len(hap1)) * 100

# read_matrix = read_input('./medium/example/input/reads_low_error.txt')
read_matrix = read_input('./large_low_error_2_chromosomes_training_reads.txt')

# set size of read_matrix to new size
lines = len(read_matrix)
# check how long the haplotype should be later
length = len(read_matrix[0])

# make new matrix with longest reads from each subset
new_matrix = newFiltered(read_matrix, length)

# print 'new_matrix '
# for c in new_matrix[0:500]:
#     print c

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

print "\n"
# Final_2, low error
# answer1 = "0101000111001010110011110111100001111011010101110001011101100000001111000000110111110000010010101000000111010011110010110010010100100000011010010111000111111000110011000001011101011011010000100011111100010000010100010011011011111011010010110010110101111011101100111000101110001100100100100000101110000110011111011010110000110001010110111101001010111101101100001001011000111011100010111000000100001001101000011110110110100111000111011001101010001101000111110111000110000101111001110001010000101011100111110001101101100111011111111010110011011011100011010110000000000010111000011000011110000010011000001001011100000110010100001011110100100101000110001011100110111010100000100010011001001001110010110111010000111101111001111010001110111100001101011011001111100001100110101000000011010101110001100110111111000001001100110011101000110101101100100001110011101011110011010101110001110101001111011100100010010000010001001001100001011110100010010101111001110101111010010000001010101110111011111100001110011110"
# answer2 = "1010111000110101001100001000011110000100101010001110100010011111110000111111001000001111101101010111111000101100001101001101101011011111100101101000111000000111001100111110100010100100101111011100000011101111101011101100100100000100101101001101001010000100010011000111010001110011011011011111010001111001100000100101001111001110101001000010110101000010010011110110100111000100011101000111111011110110010111100001001001011000111000100110010101110010111000001000111001111010000110001110101111010100011000001110010010011000100000000101001100100100011100101001111111111101000111100111100001111101100111110110100011111001101011110100001011011010111001110100011001000101011111011101100110110110001101001000101111000010000110000101110001000011110010100100110000011110011001010111111100101010001110011001000000111110110011001100010111001010010011011110001100010100001100101010001110001010110000100011011101101111101110110110011110100001011101101010000110001010000101101111110101010001000100000011110001100001"
# all 100%

# data 3, easy low error
# answer1 = "1011100111111101110100001000100101011011111011010111110000000100101000101100000110000001001101010010"
# answer2 = "0100011000000010001011110111011010100100000100101000001111111011010111010011111001111110110010101101"
# all 100%

# data 3, large low error
answer1 = "1100000011100000101001111110111110010001000101111100100101001100110011001010100011001100000001111111010111010000111001100110100111001110001001110100110010011010110011010011101001101101010110001001001101100010000001000011011111101111010010111111101010101000010110010111110001101001111110100111111001000110100100110100000100011101111111100000100100011000101010011111111111011010011001010111000011111100011111010100000111110111110000000000001011001100011001111111011111000001111101011001010110100111110001010100000001000110011101110101001010001100111110100010010000111110100110111001010010110100100010111011100101111100100100110010001100010010100101100000001110101011110010010011111000001101010110101001001011110011110011011000111010010110100011010110110011101100000111100111010000110100000111110111011001110000001111100100110100000100100000001010001000010010001001011100011011011001001110110110010101101010110110000011111000010011010011001011000010010100011101111110100100100110011001101100100101110111"
answer2 = "0011111100011111010110000001000001101110111010000011011010110011001100110101011100110011111110000000101000101111000110011001011000110001110110001011001101100101001100101100010110010010101001110110110010011101111110111100100000010000101101000000010101010111101001101000001110010110000001011000000110111001011011001011111011100010000000011111011011100111010101100000000000100101100110101000111100000011100000101011111000001000001111111111110100110011100110000000100000111110000010100110101001011000001110101011111110111001100010001010110101110011000001011101101111000001011001000110101101001011011101000100011010000011011011001101110011101101011010011111110001010100001101101100000111110010101001010110110100001100001100100111000101101001011100101001001100010011111000011000101111001011111000001000100110001111110000011011001011111011011111110101110111101101110110100011100100100110110001001001101010010101001001111100000111101100101100110100111101101011100010000001011011011001100110010011011010001000"
# all 100%

# data, reads of different lengths


print percentOverlap(hap1, answer1)
print percentOverlap(hap2, answer2)