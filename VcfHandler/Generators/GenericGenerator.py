# -*-coding:Utf-8 -*


"""
Interface of a generic file generator from VCF data reader.
"""


from abc import ABC, abstractclassmethod
from VariantCallSet.VariantCall import VariantCall


class GenericGenerator(ABC):
    """
    Interface of a generic file generator from VCF data reader.
    """


    # Constructor

    def __init__(self, **kwargs):
        """
        Constructor
        """

        return super().__init__(**kwargs)



    # Methods

    @abstractclassmethod
    def AddVariantCall(self, variant_call: VariantCall):
        """
        Adds the provided vairant call to the file content
        """
        pass


    @abstractclassmethod
    def GenerateOutputfile(self, file_name: str):
        """
        Generated the output file in the format of the reader
        """
        pass