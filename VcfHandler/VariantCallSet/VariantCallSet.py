# -*-coding:Utf-8 -*


"""
Represents a variant call set.
"""


from VariantCallSet.IndividualCallValue import IndividualCallValue
from VariantCallSet.VariantCall import VariantCall


class VariantCallSet:
    """
    Represents a variant call set.
    """


    # Constructor

    def __init__(self, **kwargs):
        """
        Constructor
        """

        # programs used to generate this variant call set
        self.sources = []

        # dictionary of filters applied to this variant call set
        self.filters = {}

        # dictionary of formats used in this variant call set
        self.formats = {}

        # dictionary of infos used in this variant call set
        self.infos = {}

        # list of contigs in this variant call set
        self.contigs = []

        # list of variant calls of this variant call set
        self.variant_calls = []



    # Methods

    def AddVariantCall(self, variant_call: VariantCall):
        """
        Adds the provided variant call to the list of variant call sets
        """
        self.variant_calls.append(variant_call)