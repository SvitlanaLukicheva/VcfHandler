# -*-coding:Utf-8 -*


"""
Implements a reader of VCF files.
"""


from VariantCallSet.IndividualCallValue import IndividualCallValue
from VariantCallSet.VariantCall import VariantCall
from VariantCallSet.VariantCallSet import VariantCallSet

import re



class VcfReader:
    """
    Implements a reader of VCF files.
    """


    # Constructor

    def __init__(self, **kwargs):
        """
        Constructor
        """

        self.variant_call_set = VariantCallSet()



    # Methods

    def ReadFile(self, file_name: str):
        """
        Reads the VCF filename provided in parameter
        """

        file = open(file_name, "r")

        for line in file:
            self.ReadLine(line)

        file.close()


    def ReadLine(self, line: str):
        """
        Reads the provided line from a VCF file
        """

        # this is a variant call
        if(not str.startswith(line, "#")):
            line_parts = re.split(r' +', line)  # elements on a line a separated by a variable number of empty spaces

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

            print("I read " + variant_call.ToString())