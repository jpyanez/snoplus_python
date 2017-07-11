import numpy as np
import ast

def readfile(infile_name):
    infile = open(infile_name)

    data = {}
    skip_symbols = list('{}/*')
    remove_symbols = [" ", '"', ","]
    brackets = ["[","]"]

    for one_line in infile:
        if one_line[0] in skip_symbols:
            continue
        one_line = one_line.replace(" ", "")
        line_contents = one_line.split(':')

        # Not reading commented parts in the same line
        line_contents[1] = line_contents[1].split('/')[0]

        if '"' in line_contents[1]:
            print 'String ', line_contents[0]
            data[line_contents[0]] = \
                                     line_contents[1][:-1].translate(None, ''.join(remove_symbols))
        elif ',' in line_contents[1][:-2]:
            print 'Array ', line_contents[0]

            # This is an array of numbers
            aux = line_contents[1].translate(None, ''.join(brackets))
            aux = aux.replace(',,','')
            data[line_contents[0]] = np.array(ast.literal_eval(aux[:-1]))
        else:
            print 'Number ', line_contents[0]
            # This is a simple number
            data[line_contents[0]] = np.float(line_contents[1][:-2])    

    return data
