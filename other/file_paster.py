import sys

res = {}
name = {}
test = {}
tab = '\t'
newline = '\n'
f = open(sys.argv[1], 'w+')
for file_name in sys.argv[2:]:
    for line_nr, line in enumerate(open(file_name)):
        for i in range(1, 5):
            res.setdefault(line_nr, []).append(line.strip().split(tab)[i])
for line_nr, line in enumerate(open(sys.argv[2])):
    name.setdefault(line_nr, []).append(line.strip().split(tab)[0])

f.write(newline)
for line_nr in sorted(res):
    f.write(tab.join(name[line_nr]))
    f.write(tab)
    f.write(tab.join(res[line_nr]))
    f.write(newline)
f.close()
