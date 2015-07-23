#!/usr/bin/python
# This script automatically creates a list of examples by reading the header in all problem.c files.
import glob

with open("README.rst") as f:
    readme = f.readlines()

keep_lines_after_header = 5
with open("README.rst","w") as f:
    start_delete = -1
    for i in xrange(0,len(readme)):
        if readme[i]=="Examples\n":
            if readme[i+1]=="========\n":
                start_delete = i + keep_lines_after_header
        if start_delete!=-1 and readme[i+2]=="-----------------------\n":
            for problemc in glob.glob("./examples/*/problem.c"):
                will_output = 0
                with open(problemc) as pf:
                    did_output=0
                    empty_lines = 0
                    for line in pf:
                        if line[0:3] == "/**":
                            will_output += 1
                        if line[0:3] == " */":
                            will_output = -1
                        if will_output>1:
                            if will_output == 2:
                                line = "   * **" + line[3:].strip() + "**"
                            will_output = 2
                            if len(line[3:].strip())==0:
                                f.write("\n\n  "+line[3:].strip())
                            else:
                                f.write(line[3:].strip())
                        if will_output==-1:
                            f.write("\n\n  Directory: examples/"+problemc.split("/")[2]+"\n\n")
                            will_output= -2
                            did_output = 1
                        if will_output>0:
                            will_output += 1
                    if did_output==0:
                        print "Warning: Did not find description in "+problemc
            start_delete = -1
        if start_delete==-1 or i<start_delete:
            f.write(readme[i])

import re
with open("src/rebound.c") as f:
    reboundlines = f.readlines()

for i,l in enumerate(reboundlines):
    if "**VERSIONLINE**" in l:
        reboundlines[i] = "const char* reb_version_str = \"2.0.0\";			// **VERSIONLINE** This line gets updated automatically. Do not edit manually."
