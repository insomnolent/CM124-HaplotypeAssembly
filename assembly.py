import os
os.chdir('/Users/christinesun/PycharmProjects/Final Project 2/data')
num_lines = sum(1 for line in open('small_no_error_training_reads.txt'))
data = open('small_no_error_training_reads.txt')
matrix = data.read()
data.close()
# to check file
print(matrix)
print num_lines


# read in first line
# check with second line
# if overlap then consider them to be the same haplotype
# if differences then set them to be different haplotypes
# continue reading in line by line

