import numpy as np
import ast
from copy import deepcopy

def readfile(infile_name):
    infile = open(infile_name)

    data = {}
    skip_symbols = list('{}/* \n')
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

def dbreader(infile_name):
    infile = open(infile_name)

    skip_lines = ['/']
    data = {}
    for i, line in enumerate(infile):
        line = line.replace('\n', '')
        line = line.replace('"', '')
        line = line.replace(' ', '')
    
        if len(line) == 0:
            continue
        if line[0] in skip_lines:
            continue

        if '{' in line[0]:
            this_dict = {}
            print 'reader: Creating a dict in line', i
        elif '}' in line[0]:
            print 'reader: Closing a dict in line', i
            data[this_dict['index']] = deepcopy(this_dict)
        else:
            split_line = line.split(':')
            this_dict[split_line[0]] = split_line[1][:-1]
        
            if split_line[1][0] == '[':
                # This is an array
                val = np.fromstring(this_dict[split_line[0]][1:-1], sep=',')
                this_dict[split_line[0]] = val
            elif split_line[1][0].isdigit():
                # This is a number
                this_dict[split_line[0]] = float(this_dict[split_line[0]])

    infile.close()
    return data
