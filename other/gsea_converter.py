import sys

res = {}
name = {}
test = {}
tab = '\t'
newline = '\n'
numberLines = 0
sep = "."
duplicates = []
# takes first file name of input, creates new file that it will write to
f = open(sys.argv[1], 'w+')
with open(sys.argv[2],'r') as contents:
      save = contents.read()
with open(sys.argv[2],'w') as contents:
      contents.write(tab)
with open(sys.argv[2],'a') as contents:
      contents.write(save)
# for every other input file name
for file_name in sys.argv[2:]:
    # generates a list, numbering the lines and the corresponding row
    for line_nr, line in enumerate(open(file_name)):
            # creates lits of all items (after the first) in each column
            res.setdefault(line_nr, []).append(line.strip().split(tab)[4])
            # for each line in the second input file
for line_nr, line in enumerate(open(sys.argv[2])):
    # generates a list that is just the name column
    name.setdefault(line_nr, []).append(line.strip().split(sep)[0].split(tab)[0])
    numberLines = line_nr
print(numberLines)

for line_nr in sorted(res):
    if line_nr < 58780:
        if name[line_nr] == name[line_nr + 1]:
            duplicates.append(line_nr+1)
duplicates.reverse()
print(len(duplicates))

for i in duplicates:
       name.pop(i)
       res.pop(i)

# writes a new lines
for line_nr in sorted(res):
    # writes the name
    f.write(tab.join(name[line_nr]))
    # adds a tab
    f.write(tab)
    f.write(tab.join(res[line_nr]))
    f.write(newline)
f.close()
