# RPC 170920
# This will make a pseudo PED file for somalier to run.
# Uses the file naming system to correctly name the sample and family names
#!/usr/bin/env python

import argparse
import os

parser = argparse.ArgumentParser()
files = parser.add_argument('--files', required=True)
files2 = files.split(" ")

for f in files2:
    list = f.split("-")
    sample = list[0]
    family = list[1]

    PED = family + "\t" + sample + "\t0\t0\t1\t2"

    cmd = "echo '" + PED + "' >> project.PED"
    os.system(cmd)
    print(PED)
