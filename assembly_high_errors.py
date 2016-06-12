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
    print "index", index
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
    recentread1 = ""
    recentread2 = ""
    for i in range(0, length):
        readset1 = takeSubset(reads, i)
        topreads = gettopreads(readset1)
        # if beginning of read matrix then set recentreads to haplotype beginnings
        if i == 0:
            recentread1 = topreads[0]
            recentread2 = topreads[1]

        if len(topreads) == 1:
            if overlap(read1, recentread1, i) or overlap(read1, recentread2, i):
                new.append(topreads[0])
                recentread1 = topreads[0]
        elif len(topreads) > 1:
            read1 = topreads[0]
            read2 = topreads[1]

            if checkInverse(removeDash(read1), removeDash(read2)):
                #if (overlap(recentread1, read1, i) and overlap(recentread2, read2, i)) or (overlap(recentread1, read2, i) and overlap(recentread2, read1, i)):
                print read1
                print read2
                recentread1 = topreads[0]
                new.append(topreads[0])
                recentread2 = topreads[1]
                new.append(topreads[1])
            else:
                # if the most frequent reads aren't inverses of each other
                for j in range(0, len(topreads)-1):
                    read1 = topreads[j]
                    # if overlap(recentread1, read1, i) or overlap(recentread2, read1, i):
                    #     recentread1 = read1
                    #     new.append(read1)
                    #


                    for k in range(j+1, len(topreads)):
                        read2 = topreads[k]
                        if checkInverse(removeDash(read1), removeDash(read2)):
                        #     print "does it reach here"
                        #     print recentread1
                        #     print 'read1', read1
                        #     print recentread2
                        #     print 'read2', read2

                        # if (overlap(recentread1, read1, i) and overlap(recentread2, read2, i)) or (overlap(recentread1, read2, i) and overlap(recentread2, read1, i)):
                            recentread1 = read1
                            new.append(read1)
                            recentread2 = read2
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


# test if a read at that index is overlapped by the read above it
def overlap(first, second, index):
    result = True
    end = index + len(removeDash(first))-1
    for i in range(index, end):
        if first[i] != second[i]:
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


# calculates the percent similarity between haplotypes
def percentOverlap(hap1, hap2):
    percent = 0
    track = 0
    for c in hap1:
        if c == hap2[track]:
            percent += 1
        track += 1
    return float(percent)/float(len(hap2)) * 100


def buildmatrix(reads, hap):
    current = hap
    new_matrix = []
    new_matrix.append(hap)
    for k in range(1, len(reads)):
        subset = takeSubset(reads, k)
        for line in subset:
            if percentOverlap(current, line) > 80:
                new_matrix.append(line)
                current = line


read_matrix = read_input('./medium/example/input/reads_high_error.txt')
# read_matrix = read_input('./easy_high_error_2_chromosomes_training_reads.txt')

# set size of read_matrix to new size
lines = len(read_matrix)
# check how long the haplotype should be later
length = len(read_matrix[0])

# set haplotypes 1 and 2
firstset = takeSubset(read_matrix, 0)
firsttopreads = gettopreads(firstset)

if len(firsttopreads) == 1:
        hap1 = firsttopreads[0]
        hap2 = inverse(hap1)
elif len(firsttopreads) > 1:
    read1 = firsttopreads[0]
    read2 = firsttopreads[1]

    if checkInverse(removeDash(read1), removeDash(read2)):
        hap1 = firsttopreads[0]
        hap2 = firsttopreads[1]
    else:
        # if the most frequent reads aren't inverses of each other
        for j in range(0, len(firsttopreads) - 1):
            read1 = firsttopreads[j]
            for k in range(j + 1, len(firsttopreads)):
                read2 = firsttopreads[k]
                if checkInverse(removeDash(read1), removeDash(read2)):
                    hap1 = firsttopreads[j]
                    hap2 = firsttopreads[k]
                    break

# build matrix for each haplotype
hap1matrix = buildmatrix(read_matrix, hap1)




# for testing purposes
for i in range(0, len(hap1)):
    sys.stdout.write(hap1[i])

print "\n"
for i in range(0, len(hap2)):
    sys.stdout.write(hap2[i])

print "\n"
# Final_2, high error data - currently 50.4% and 50.6%
answer1 = "0100101000011100011000010000101101010000111111011101111110000100001000000011001100100001001001010111011110001100101001100111101100001010101000000100011101101101010010011011111100000101000101111111010011101010001000001110100011110111000111111110111111001101000110111101011011111100101101101011011100110110110111101011110110101010101111011011000100010010100101111100000111011010010001111101010000011101100101101010001110011101101001100000100010000000011100011011000000011000010111001110111001000000110100111111111110011001000010010011001110111011000000110000110111010001110110001000011101101110001111110010111011001100110001110000111110000111110010000100101110011001110111110010101100011001111101000110110100010111001010111011011101011000011110111111111101001100100011011100101010001001000100001011001110111000101001000100000110100100100110000000101100110001011011010100000000110110111001100100111011010110110100001100010010111110100011011111101001100000101010110100010010000011001110111010101101101001"
answer2 = "1011010111100011100111101111010010101111000000100010000001111011110111111100110011011110110110101000100001110011010110011000010011110101010111111011100010010010101101100100000011111010111010000000101100010101110111110001011100001000111000000001000000110010111001000010100100000011010010010100100011001001001000010100001001010101010000100100111011101101011010000011111000100101101110000010101111100010011010010101110001100010010110011111011101111111100011100100111111100111101000110001000110111111001011000000000001100110111101101100110001000100111111001111001000101110001001110111100010010001110000001101000100110011001110001111000001111000001101111011010001100110001000001101010011100110000010111001001011101000110101000100100010100111100001000000000010110011011100100011010101110110111011110100110001000111010110111011111001011011011001111111010011001110100100101011111111001001000110011011000100101001001011110011101101000001011100100000010110011111010101001011101101111100110001000101010010010110"

# data 3, easy high error - currently 87%
# answer1 = "0101111100111101101111010010001000000001100101111100101011000111000000000100010110011101111000100001"
# answer2 = "1010000011000010010000101101110111111110011010000011010100111000111111111011101001100010000111011110"

print percentOverlap(hap1, answer1)
print percentOverlap(hap2, answer2)
