import os
import sys
import fileinput

for line in fileinput.input():
  l = line[:-1]
  print ("efetch -db taxonomy -id " + l + " -format native -mode xml | xtract -pattern TaxaSet -element TaxId,ScientificName,GenbankCommonName,Division")

