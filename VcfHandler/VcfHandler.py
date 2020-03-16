# -*-coding:Utf-8 -*


from Readers.VcfReader import VcfReader
from Generators.SfsGenerator import SfsGenerator


vcf_reader = VcfReader()

sfs_generator = SfsGenerator("C:\\Users\\svitl\\Desktop\\test_sfs_generator\\popfile_dadi.txt")

vcf_reader.ReadFile("C:\\Users\\svitl\\Desktop\\test_sfs_generator\\example_dadi.vcf", sfs_generator)

sfs_generator.GenerateOutputfile("result.out")