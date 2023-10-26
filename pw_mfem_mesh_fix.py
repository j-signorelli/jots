#!/usr/bin/env python3
"""
Script to fix PW-exported MFEM mesh files
"""

from optparse import OptionParser
import os

def decrementIndices(line_arr, startingIndex):

    for i in range(startingIndex, len(line_arr)):
        line_arr[i] = str(int(line_arr[i]) - 1)

    newline = " ".join(line_arr)
    
    return newline

def main():

    parser = OptionParser()
    parser.add_option("-f", "--file", dest="filename", help="Mesh file to fix", metavar="FILE")

    (options, args) = parser.parse_args()

    max_bdr_attr = 0
    reading_elements = False
    reading_boundaries = False

    fixed_filename = str("fixed_" + options.filename)
    w = open(fixed_filename, "w")
    f = open(options.filename, "r")

    for i, line in enumerate(f):
        
        if i == 2:
            w.write("# PW Export Fix Script Applied To This File\n")

        if reading_elements or reading_boundaries:
            # Update vertex indices

            # If have reached blank line, continue
            if line == "" or line == "\n":
                reading_elements = False
                reading_boundaries = False
                w.write(line)
                continue
            

            # Skip number of elements/bdrs line
            line_arr = line.split()
            if len(line_arr) == 1:
                w.write(line)
                continue
            
            if reading_boundaries:
                # Get max bdr attribute
                bdr_attr = int(line_arr[0])
                if bdr_attr > max_bdr_attr:
                    max_bdr_attr = bdr_attr

            # Decrement indices
            w.write(decrementIndices(line_arr,2) + "\n")
            
        else:

            if "elements" in line:
                reading_elements = True

            if "boundary" in line:
                reading_boundaries = True

            w.write(line)

    w.close()
    f.close()

    # Update all max integer element attributes to be max_bdr_attr + 1
    # Credit to mfem/mfem#3927
    os.system("sed -i 's/2147483647/" + str(max_bdr_attr+1) + "/g' " + fixed_filename)

# Can run only from terminal
if __name__ == '__main__':
    main()
