import time
import itertools
import numpy as np
from collections import Counter
from functools import reduce
import OSprocess
import Replace
import SNPallele
from Config import Setting

class Parameter:
    '''
    input: microhaplotype allele lists
    output: values of forensic parameters of microhaplotypes
    '''

    def __init__(self, listMHallele):
        self.listMHallele = listMHallele  # listMHallele stores a list of all alleles that contain a MH,
                                        # the first and second being an individual, and so on.

    def MHalleleallele(self):
        '''
        输出微单倍型等位基因类型
        '''
        MHallele = list(Counter(self.listMHallele).keys())
        return MHallele

    def MHalleleAF(self):
        '''
        Output microhaplotype allele frequency
        '''
        MHalleleAF = []
        for i in list(Counter(self.listMHallele).values()):
            i = str(i / lenthMHtmp)
            MHalleleAF.append(i)
        MHalleleAF = list(map(float, MHalleleAF))
        return MHalleleAF

    def MHAEcalculate(self):
        '''
        Calculate the AE value for this microhaplotype
        '''
        MHAF = self.MHalleleAF()
        AEtmp = []
        for i in MHAF:
            i = float(i) ** 2
            AEtmp.append(i)
        AE = [str(1 / sum(AEtmp))]
        return AE

    def MHHetcalculate(self):
        '''
        Calculate the heterozygosity for this microhaplotype
        '''
        i = 0
        hete = 0
        lenthMHallele = len(self.listMHallele)
        while i < lenthMHallele:
            if self.listMHallele[i] != self.listMHallele[i + 1]:
                hete = hete + 1
            i += 2
        heterozygosity = [str(hete * 2 / lenthMHallele)]
        return heterozygosity

    def MHPDcalculate(self):
        '''
        Calculate PD for this microhaplotype
        '''
        i = 0
        genolist = []
        genoallelelist = []
        lenthMHallele = len(self.listMHallele)
        while i < lenthMHallele:
            geno1 = self.listMHallele[i] + self.listMHallele[i + 1]
            if geno1 in genolist:
                genolist.append(geno1)
            elif geno1 not in genolist:
                geno2 = self.listMHallele[i + 1] + self.listMHallele[i]
                if geno2 in genolist:
                    genolist.append(geno2)
                elif geno2 not in genolist:
                    genolist.append(geno1)
            i += 2
        for i in list(Counter(genolist).values()):
            i = (i * 2 / lenthMHallele) ** 2
            genoallelelist.append(i)
        DP = [str(1 - sum(genoallelelist))]
        return DP

    def SNPafcalculate(self, listPOStmp):
        '''
        Using MH's SNP list, output a list containing only SNP allele frequencies to facilitate subsequent calculation
        of minimum and maximum allele frequencies
        '''
        listPOSaftmp = []
        for snp in listPOStmp[:]:
            if len(snp) == 5:
                del snp[0]
                del snp[0]
                del snp[1]
            if len(snp) == 7:
                del snp[0]
                del snp[0]
                del snp[1]
                del snp[2]
            if len(snp) == 9:
                del snp[0]
                del snp[0]
                del snp[1]
                del snp[2]
                del snp[3]
            listPOSaftmp.append(snp)
        listPOSaftmp = Replace.LtF(listPOSaftmp)
        return listPOSaftmp

    def MHalleleINFO(self):
        '''
        Generating a dictionary of microhaplotype alleles and allele frequencies
        '''
        MHallele = self.MHalleleallele()
        MHAF = self.MHalleleAF()  
        MHalleledic = dict(zip(MHallele, MHAF))
        return MHalleledic
