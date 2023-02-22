import sys

# The first argument passed into this program will be the name of the dataset of the samples
dataset_name = sys.argv[1] + ".Count_Matrix.tsv"

# The rest of the arguments will be the files which want to be pasted together. 
# Because this python script will be called from a bash shell, this will effectively be one argument in the form of a bash array.
file_names = sys.argv[2:]
sample_names = []

# Here we get the pure names of the samples so that they can used as a header for the gene count matrix

for file in file_names:
    file_array = file.split(".")
    file = file_array[0]
    sample_names.append(file)

res = {}
name = {}
test = {}
tab = '\t'
newline = '\n'
f = open(dataset_name, 'w+')

for file_name in file_names:
    for line_number, line in enumerate(open(file_name)):
        res.setdefault(line_number, []).append(line.strip().split(tab)[1])

for line_number, line in enumerate(open(file_names[1])):
    name.setdefault(line_number, []).append(line.strip().split(tab)[0])

f.write("Genes")

for sample in sample_names:
    f.write(tab)
    f.write(sample)
    
f.write(newline)

for line_number in sorted(res):
    f.write(tab.join(name[line_number]))
    
    f.write(tab)

    f.write(tab.join(res[line_number]))

    f.write(newline)

f.close()

