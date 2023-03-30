#!/usr/bin/env python3
import sys
import os

samples = "Our_Zou_Kucab_SNVs.txt"
signatures = "SBS96_S10_Signatures.txt"
output = "SBS96_10_Signatures/Decomposition"

from SigProfilerAssignment import Analyzer as Analyze
def main():
  Analyze.decompose_fit(samples,output,signatures=signatures,genome_build="GRCh38",verbose=True)
if __name__ == '__main__':
        main()


samples = "Our_Zou_Kucab_indels.txt"
signatures = "ID83_S4_Signatures.txt"
output = "ID83_4_Signatures//Decomposition"

from SigProfilerAssignment import Analyzer as Analyze
def main():
  Analyze.decompose_fit(samples,output,signatures=signatures,genome_build="GRCh38",verbose=True)
if __name__ == '__main__':
        main()

