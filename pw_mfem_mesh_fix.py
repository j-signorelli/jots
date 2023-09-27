#!/usr/bin/env python3
"""
Script to fix PW-exported MFEM mesh files
"""

from optparse import OptionParser

def decrementIndices(line, startingIndex):
    line_arr = line.split()

    for i in range(startingIndex, len(line_arr)):
        line_arr[i] = str(int(line_arr[i]) - 1)

    newline = " ".join(line_arr)
    
    return newline

def main():

    parser = OptionParser()
    parser.add_option("-f", "--file", dest="filename", help="Mesh file to fix", metavar="FILE")

    (options, args) = parser.parse_args()


    reading_elements = False
    reading_boundaries = False

    w = open(str("fixed_" + options.filename), "w")
    f = open(options.filename, "r")

    for i, line in enumerate(f):
        
        if i == 2:
            w.write("# PW Export Fix Script Applied To This File\n")

        if reading_elements:
            # Update vertex indices

            # If have reached blank line, continue:
            if line == "" or line == "\n":
                reading_elements = False
                continue
            
            w.write(decrementIndices(line,2) + "\n")
            
        elif reading_boundaries:
            # If have reached blank line, continue:
            if line == "" or line == "\n":
                reading_boundaries = False
                continue
            
            if "vertices" in line:
                breakpoint()

            w.write(decrementIndices(line,2) + "\n")
        else:

            if "elements" in line:
                reading_elements = True

            if "boundary" in line:
                reading_boundaries = True

            w.write(line)

    w.close()
    f.close()

# Can run only from terminal
if __name__ == '__main__':
    main()
