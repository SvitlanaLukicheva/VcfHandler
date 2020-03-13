# -*-coding:Utf-8 -*


from Readers.VcfReader import VcfReader
from Generators.SfsGenerator import SfsGenerator


vcf_reader = VcfReader()

sfs_generator = SfsGenerator("C:\\Users\\svitl\\Desktop\\test_sfs_generator\\population.txt")

vcf_reader.ReadFile("C:\\Users\\svitl\\Desktop\\test_sfs_generator\\test.vcf", sfs_generator)

sfs_generator.GenerateOutputfile("C:\\Users\\svitl\\Desktop\\test_sfs_generator\\result.out")