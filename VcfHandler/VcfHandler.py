# -*-coding:Utf-8 -*


from Readers.VcfReader import VcfReader
from Generators.SfsGenerator import SfsGenerator


vcf_reader = VcfReader()

test_example_file = False
test_full_vcf_with_t10_processed = False
test_full_vcf_without_t10_processed = False
test_vcf_without_t10 = False
test_gi_vcf = False
test_gq_vcf = True

if(test_example_file == True):
    sfs_generator = SfsGenerator("C:\\Users\\svitl\\Desktop\\test_sfs_generator\\example_dadi\\popfile_dadi.txt")
    vcf_reader.ReadFile("C:\\Users\\svitl\\Desktop\\test_sfs_generator\\example_dadi\\example_dadi.vcf", sfs_generator)
    sfs_generator.GenerateOutputfile("example")

if(test_full_vcf_with_t10_processed == True):
    sfs_generator = SfsGenerator("C:\\Users\\svitl\\Desktop\\test_sfs_generator\\popfile_sv_full_with_t10.txt")
    vcf_reader.ReadFile("C:\\Users\\svitl\\Desktop\\test_sfs_generator\\filtered_cohort_no_nocalls.vcf", sfs_generator)
    sfs_generator.GenerateOutputfile("result_full_with_t10")

if(test_full_vcf_without_t10_processed == True):
    sfs_generator = SfsGenerator("C:\\Users\\svitl\\Desktop\\test_sfs_generator\\popfile_sv_full.txt")
    vcf_reader.ReadFile("C:\\Users\\svitl\\Desktop\\test_sfs_generator\\filtered_cohort_no_nocalls.vcf", sfs_generator)
    sfs_generator.GenerateOutputfile("result_full")

if(test_vcf_without_t10 == True):
    sfs_generator = SfsGenerator("C:\\Users\\svitl\\Desktop\\test_sfs_generator\\popfile_sv_no_t10.txt")
    vcf_reader.ReadFile("C:\\Users\\svitl\\Desktop\\test_sfs_generator\\filtered_cohort_no_nocalls_no_t10.vcf", sfs_generator)
    sfs_generator.GenerateOutputfile("result_no_t10")

if(test_gi_vcf == True):
    sfs_generator = SfsGenerator("C:\\Users\\svitl\\Desktop\\test_sfs_generator\\popfile_gi.txt")
    vcf_reader.ReadFile("G:\\introgression\\results\\gatk\\select_variants\\full_no_phasing\\gi\\remove_nocalls\\result\\gi_filtered_cohort_no_nocalls.vcf", sfs_generator)
    sfs_generator.GenerateOutputfile("result_gi")

if(test_gq_vcf == True):
    sfs_generator = SfsGenerator("C:\\Users\\svitl\\Desktop\\test_sfs_generator\\popfile_gq.txt")
    vcf_reader.ReadFile("G:\\introgression\\results\\gatk\\select_variants\\full_no_phasing\\gq\\remove_nocalls\\result\\gq_filtered_cohort_no_nocalls.vcf", sfs_generator)
    sfs_generator.GenerateOutputfile("result_gq")