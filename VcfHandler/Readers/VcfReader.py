# -*-coding:Utf-8 -*


"""
Implements a reader of VCF files.
"""


from Generators.SfsGenerator import SfsGenerator
from VariantCallSet.IndividualCallValue import IndividualCallValue, Genotype
from VariantCallSet.VariantCall import VariantCall
from VariantCallSet.VariantCallSet import VariantCallSet

from datetime import datetime

import re



class VcfReader:
    """
    Implements a reader of VCF files.
    """


    # Constructor

    def __init__(self, enable_debug: bool = False, **kwargs):
        """
        Constructor
        """

        self.variant_call_set = VariantCallSet()

        self.sfs_generator = None

        self.enable_debug = enable_debug



    # Methods

    def ReadFile(self, file_name: str, sfs_generator: SfsGenerator = None):
        """
        Reads the VCF filename provided in parameter
        """

        print(datetime.now().strftime("%H:%M:%S") + ": Reading file " + file_name + "...")

        self.sfs_generator = sfs_generator

        file = open(file_name, "r")
        line_cnt = 0
        for line in file:
            self.ReadLine(line)
            line_cnt += 1
            if(line_cnt % 10000 == 0):
                print(datetime.now().strftime("%H:%M:%S") + ": Read " + str(line_cnt) + " lines...")
        file.close()

        print(datetime.now().strftime("%H:%M:%S") + ": Done.")
        


    def ReadLine(self, line: str):
        """
        Reads the provided line from a VCF file
        """

        # this is a variant call
        if(not str.startswith(line, "#")):
            line_parts = re.split(r' +', line)  # elements on a line a separated by a variable number of empty spaces

            if(len(line_parts) < 9):  # variant call must have at least 9 standard fields
                line_parts = re.split(r'\t+', line)  # then we check tab separator
                if(len(line_parts) < 9):  # variant call must have at least 9 standard fields
                    raise Exception("Incorrect format for the varian call " + line)

            # creation of a variant call with the standard fields
            variant_call = VariantCall(line_parts[0], line_parts[1], line_parts[2], line_parts[3], line_parts[4], line_parts[5], line_parts[6], line_parts[7], line_parts[8])

            # adding individual call values
            for i in range (9, len(line_parts)):
                individual_call_info = line_parts[i]
                individual_call_info_parts = individual_call_info.split(":")
                genotype = individual_call_info_parts[0]  # format: 1/0 (1: ref, 0: alt, | or /: phased or non phased)
                genotype_0 = genotype[0]
                genotype_1 = genotype[2]
                if(genotype[1] == "|" or genotype[1] == "/"):
                    is_phased = True if genotype[1] == "|" else False
                else:
                    raise Exception("Invalid genotype format: " + genotype)
                    
                variant_call.AddIndividualCall(IndividualCallValue(genotype_0, genotype_1, is_phased))

            if(self.enable_debug):
                print("I read the line \"" + line + "\"")
                print("The variant call is " + variant_call.ToFullString())

            if(self.sfs_generator != None):
                self.sfs_generator.AddVariantCall(variant_call)