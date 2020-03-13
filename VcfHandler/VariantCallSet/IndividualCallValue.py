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

    def __init__(self, genotype_0: str, genotype_1: str, is_phased: bool, **kwargs):
        """
        Constructor
        """

        self.genotype_0 = self.GenotypeFromString(genotype_0)
        self.genotype_1 = self.GenotypeFromString(genotype_1)
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


    def ToString(self):
        """
        Provides the textual representation of the individual call value
        """

        result = "[" + self.GenotypeToString(self.GetGenotype0())
        if(self.GetIsPhased() == True):
            result += "|"
        else:
            result += "/"
        result += self.GenotypeToString(self.GetGenotype1()) + "]"

        return result


    def GenotypeToString(self, genotype: Genotype):
        """
        Converts the provided genotype value to its textual representation
        """

        result = "."

        if(genotype == Genotype.REF):
            result = "0"
        elif(genotype == Genotype.ALT):
            result = "1"

        return result


    def GenotypeFromString(self, genotype_value: str):
        """
        Converts the provided genotype value from its textual representation to the Genotype enum
        """

        result = Genotype.UNKNOWN

        if(genotype_value == "0"):
            result = Genotype.REF
        elif(genotype_value == "1"):
            result = Genotype.ALT

        return result