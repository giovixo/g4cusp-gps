import os
import glob

filelist = sorted(glob.glob("isotope_list_for_each_volume/*.dat"))

for filein in filelist:
    print(filein)
    f = open(filein, "r")
    lines = f.readlines()
    f.close()
    fo = open(filein, "w")
    for line in set(lines):
        fo.write(line)
    fo.close()
