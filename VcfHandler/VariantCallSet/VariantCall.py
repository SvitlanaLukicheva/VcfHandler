# -*-coding:Utf-8 -*


"""
Represents a variant call.
"""


from VariantCallSet.IndividualCallValue import IndividualCallValue


class VariantCall:
    """
    Represents a variant call.
    """


    # Constructor

    def __init__(self, chrom: str, pos: int, id: str, ref: str, alt: str, qual: str, filter: str, info: str, format: str, **kwargs):
        """
        Constructor
        """

        self.chrom = chrom
        self.pos = pos
        self.id = id
        self.ref = ref
        self.alt = alt
        self.qual = qual
        self.filter = filter
        self.info = info
        self.format = format
        self.individual_calls = []

        self.CheckFields()



    # Getters

    def GetChromosome(self):
        """
        Returns the chromosome name
        """
        return self.chrom


    def GetPosition(self):
        """
        Returns the position of the call
        """
        return self.pos


    def GetId(self):
        """
        Returns the id of the call
        """
        return self.id


    def GetRef(self):
        """
        Returns the reference allele of the call
        """
        return self.ref


    def GetAlt(self):
        """
        Returns the alternative allele of the call
        """
        return self.alt


    def GetQuality(self):
        """
        Returns the quality of the call
        """
        return self.chrom


    def GetFilter(self):
        """
        Returns the filter-related information of the call
        """
        return self.filter


    def GetInfo(self):
        """
        Returns the information field of the call
        """
        return self.info


    def GetFormat(self):
        """
        Returns the format field of the call
        """
        return self.format


    def GetIndividualCalls(self):
        """
        Returns the list of individual calls of the variant call
        """
        return self.individual_calls



    # Methods

    def CheckFields(self):
        """
        Performs validity check on variant call's fields.
        """

        if(not str.startswith(self.GetFormat(), "GT")):
            raise SyntaxError("The format field of the variant call " + self.getChromosome() + " pos. " + self.GetPosition() + " doesn't start with a GT field")



    def AddIndividualCall(self, individual_call: IndividualCallValue):
        """
        Adds the provided individual call value to the list of individual values of the variant call
        """

        self.individual_calls.append(individual_call)


    def ToString(self):
        """
        Returns the textual representation of the variant call
        """

        return "chrom: {" + self.GetChromosome() + "} pos: {" + self.GetPosition() + "} ref: {" + self.GetRef() + "} alt: {" + self.GetAlt() + "} individuals: {" + str(len(self.GetIndividualCalls())) + "}"
