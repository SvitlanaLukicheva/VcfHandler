# -*-coding:Utf-8 -*


"""
Represents a variant call value for an individual of the cohort.
"""


from enum import Enum


class Genotype(Enum):
    REF = 1
    ALT = 2
    UNKNOWN = 3



class IndividualCallValue:
    """
    Represents a variant call value for an individual of the cohort.
    """
    

    # Constructor

    def __init__(self, genotype_0: Genotype, genotype_1: Genotype, is_phased: bool, **kwargs):
        """
        Constructor
        """

        self.genotype_0 = genotype_0
        self.genotype_1 = genotype_1
        self.is_phased = is_phased



    # Getters

    def GetGenotype0(self):
        """
        Returns the value of the first genotype of the individual
        """
        return self.genotype_0


    def GetGenotype1(self):
        """
        Returns the value of the second genotype of the individual
        """
        return self.genotype_1


    def GetIsPhased(self):
        """
        Returns a boolean indicating whether the genotype is phased
        """
        return self.is_phased