import sys

def id_simplify(input_file, output_file):
    with open(input_file,'rt') as f1:
        with open(output_file,'wt') as f2:
            for eachline in f1:
                if eachline[0] == '>':
                    f2.write(eachline.strip().split()[0])
                    f2.write('\n')
                else:
                    f2.write(eachline)

input_file = sys.argv[1]
output_file = sys.argv[2]

id_simplify(input_file, output_file)