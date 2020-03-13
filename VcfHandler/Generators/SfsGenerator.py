# -*-coding:Utf-8 -*


"""
Generator of SFS files from variant call sets.
"""


from Generators.GenericGenerator import GenericGenerator
from VariantCallSet.VariantCall import VariantCall
from VariantCallSet.IndividualCallValue import IndividualCallValue, Genotype

import re


class SfsGenerator(GenericGenerator):
    """
    Generator of SFS files from variant call sets.
    """


    # Constructor

    def __init__(self, config_file_name: str, enable_debug: bool = False, **kwargs):
        """
        Constructor
        """

        self.config_file_name = config_file_name

        # names of all the populations in the variant call set
        self.population_names = []

        # indices of individuals of the population 1 in the variant call set
        self.population_1_indices = []

        # indices of individuals of the population 2 in the variant call set
        self.population_2_indices = []

        # spectrum for the population 1, will be initialized after config parsing
        self.population_1_spectrum = []

        # spectrum for the population 2, will be initialized after config parsing
        self.population_2_spectrum = []

        self.enable_debug = enable_debug

        self.ReadConfig()

        self.InitializeSpectra()

        return super().__init__(**kwargs)



    # Methods

    def ReadConfig(self):
        """
        Reads the config file and initializes the generator.
        Config file consists on one line par individual in the variant call set.
        Each line is composed of two fields separated by a tab:
        1. individual id
        2. population id
        If population id is "N/A", the corresponding individual is not taken
        into account for the computation of the SFS file.
        """

        file = open(self.config_file_name, "r")

        index = 0

        for line in file:
            line_parts = re.split(r'\t+', line)
            if(len(line_parts) < 2):
               raise Exception("Wrong format for config line \"" + line + "\"")

            individual_name = line_parts[0]
            population_name = line_parts[1]

            if(population_name == "N/A"):
                pass

            if(not population_name in self.population_names):
                if(len(self.population_names) == 2):
                    raise Exception("Only two populations can be handled for the moment for SFS file generation")
                self.population_names.append(population_name)

            if(self.population_names.index(population_name) == 0):  # this is the first population
                self.population_1_indices.append(index)
            else:
                self.population_2_indices.append(index)

            index += 1

        file.close()

        print("Config file read. Found " + str(len(self.population_names)) + " populations. Pop 1 indices: [" + str(self.population_1_indices).strip('[]') + "], pop 2 indices: [" + str(self.population_2_indices).strip('[]') + "]");


    def InitializeSpectra(self):
        """
        Initializes the spectra of the populations after config reading
        Each spectrum contains the number of elements equal to 
        2 * its population size (since individuals are diploid) + 1 (for nocalls)
        """

        self.population_1_spectrum = [0] * (1 + 2*len(self.population_1_indices))
        self.population_2_spectrum = [0] * (1 + 2*len(self.population_2_indices))


    def AddVariantCall(self, variant_call: VariantCall):
        """
        Adds the provided vairant call to the file content
        """

        pop_1_ref_cnt = 0
        pop_1_alt_cnt = 0
        pop_2_ref_cnt = 0
        pop_2_alt_cnt = 0

        i = 0
        while i < len(variant_call.GetIndividualCalls()):
            individual_call_value = variant_call.GetIndividualCalls()[i]

            if(individual_call_value.GetGenotype0() == Genotype.REF):
                if(i in self.population_1_indices):
                    pop_1_ref_cnt += 1
                elif(i in self.population_2_indices):
                    pop_2_ref_cnt +=1
            elif(individual_call_value.GetGenotype0() == Genotype.ALT):
                if(i in self.population_1_indices):
                    pop_1_alt_cnt += 1
                elif(i in self.population_2_indices):
                    pop_2_alt_cnt +=1

            if(individual_call_value.GetGenotype1() == Genotype.REF):
                if(i in self.population_1_indices):
                    pop_1_ref_cnt += 1
                elif(i in self.population_2_indices):
                    pop_2_ref_cnt +=1
            elif(individual_call_value.GetGenotype1() == Genotype.ALT):
                if(i in self.population_1_indices):
                    pop_1_alt_cnt += 1
                elif(i in self.population_2_indices):
                    pop_2_alt_cnt +=1

            i+=1

        if(self.enable_debug):
            print("INDICES ARE: 1_ref: " + str(pop_1_ref_cnt) + ", 1_alt: " + str(pop_1_alt_cnt) + ", 2_ref: " + str(pop_2_ref_cnt) + ", 2_alt: " + str(pop_2_alt_cnt))

        self.population_1_spectrum[min(pop_1_ref_cnt, pop_1_alt_cnt)] += 1
        self.population_2_spectrum[min(pop_2_ref_cnt, pop_2_alt_cnt)] += 1

        print("Variant call added. Spectrums: pop1: [" + str(self.population_1_spectrum).strip('[]') + "], pop2: [" + str(self.population_2_spectrum).strip('[]') + "]")


    def GenerateOutputfile(self, file_name: str):
        """
        Generated the output file in the format of the reader
        """

        result = ""

        for i in range(0, len(self.population_1_spectrum)):
            for j in range(0, len(self.population_2_spectrum)):
                result += str(self.population_1_spectrum[i] + self.population_2_spectrum[j]) + " "

        print(result)