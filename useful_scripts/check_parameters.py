import re
import sys
import glob

source_files = [f for f in glob.glob("src/**", recursive=True)]
source_parlist = []
for filename in source_files:
    # skip files to ignore
    if filename[-2:]!=".h" and filename[-4:]!=".cpp":
        continue

    # search every line within file
    with open(filename,'r') as f:
        for line in f:
            if (("lua->" in line) or ("lua." in line)) and ("<" in line):
                parname = line.split("\"")[1]
                if parname not in source_parlist:
                    source_parlist.append(parname)

PARAMETERS_list = []
with open("PARAMETERS",'r') as f:
    for line in f:
        if "=" in line:
            parname = line.split("=")[0].strip()
            if(len(parname)>0 and parname[0].isalpha() and (not "<" in parname) and (not ">" in parname)):
                if(parname not in PARAMETERS_list):
                    PARAMETERS_list.append(parname)

print()
print("######################")
print("# SOURCE PARAMETERS: #")
print("######################")
print()
print(source_parlist)
print()

print("########################")
print("# OBSOLETE PARAMETERS: #")
print("########################")
print()
test_files = [f for f in glob.glob("tests/**", recursive=True)]
for filename in test_files:
    # skip files to ignore
    if filename[-4:]!=".lua":
        continue
    with open(filename,'r') as f:
        for line in f:
            if "=" in line and line.strip()[:2]!="--":
                parname = line.split("=")[0].strip()
                if parname not in source_parlist:
                    print(parname,"\tin",filename)


print()
print("###########################")
print("# UNDOCUMENTED PARAMETERS #")
print("###########################")
print()
for par in source_parlist:
    if par not in PARAMETERS_list:
        print(par)

print()
print("#############################")
print("# DOCUMENTED NON-PARAMETERS #")
print("#############################")
print()
for par in PARAMETERS_list:
    if par not in source_parlist:
        print(par)
