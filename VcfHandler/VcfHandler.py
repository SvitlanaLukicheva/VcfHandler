# -*-coding:Utf-8 -*


from Readers.VcfReader import VcfReader
from Generators.SfsGenerator import SfsGenerator


# Example of code allowing to read a VCF file, a population file, and generated an SFS

vcf_reader = VcfReader()

sfs_generator = SfsGenerator("example-popfile.txt")
vcf_reader.ReadFile("path-to-vcf-file", sfs_generator)
sfs_generator.GenerateOutputfile("example-sfs")