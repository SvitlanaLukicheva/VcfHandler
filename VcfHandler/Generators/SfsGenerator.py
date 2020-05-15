# -*-coding:Utf-8 -*


"""
Generator of SFS files from variant call sets.
"""


from Generators.GenericGenerator import GenericGenerator
from VariantCallSet.VariantCall import VariantCall
from VariantCallSet.IndividualCallValue import IndividualCallValue, Genotype

from datetime import datetime

import re
import math


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

        # joint spectrum of both populations: a two-dimensional matrix that will be initialized after config parsing
        self.two_populations_spectrum = []

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

        print(datetime.now().strftime("%H:%M:%S") + ": Reading SFS config...")

        file = open(self.config_file_name, "r")

        index = 0

        for line in file:
            line_parts = re.split(r'\t+', line.replace("\n", ""))
            if(len(line_parts) < 2):
               raise Exception("Wrong format for config line \"" + line + "\"")

            individual_name = line_parts[0]
            population_name = line_parts[1]

            if(population_name != "N/A"):
                if(not population_name in self.population_names):
                    if(self.NumberOfPopulations() == 2):
                        raise Exception("Only two populations can be handled for the moment for SFS file generation")
                    self.population_names.append(population_name)

                if(self.population_names.index(population_name) == 0):  # this is the first population
                    self.population_1_indices.append(index)
                else:
                    self.population_2_indices.append(index)

            index += 1

        file.close()

        if(self.enable_debug):
            print("Config file read. Found " + str(len(self.population_names)) + " populations. Pop 1 indices: [" + str(self.population_1_indices).strip('[]') + "], pop 2 indices: [" + str(self.population_2_indices).strip('[]') + "]");

        print(datetime.now().strftime("%H:%M:%S") + ": Done.")



    def NumberOfPopulations(self):
        """
        Returns the number of populations of the SFS file
        """
        return len(self.population_names)



    def InitializeSpectra(self):
        """
        Initializes the spectra of the populations with 0 after config reading
        Each spectrum contains the number of elements equal to 
        2 * its population size (since individuals are diploid) + 1 (for nocalls)
        Dimensions of two populations spectrum correspond to the dimensions of each spectrum
        """

        self.population_1_spectrum = [0] * (1 + 2*len(self.population_1_indices))
        if(self.NumberOfPopulations() > 1):
            self.population_2_spectrum = [0] * (1 + 2*len(self.population_2_indices))
            self.two_populations_spectrum = [0] * len(self.population_1_spectrum)
            for i in range (0, len(self.population_1_spectrum)):
                self.two_populations_spectrum[i] = [0] * len(self.population_2_spectrum)



    def AddVariantCall(self, variant_call: VariantCall):
        """
        Adds the provided vairant call to the file content
        """

        pop_1_alt_cnt = 0
        pop_2_alt_cnt = 0

        i = 0
        while i < len(variant_call.GetIndividualCalls()):
            individual_call_value = variant_call.GetIndividualCalls()[i]

            if(individual_call_value.GetGenotype0() == Genotype.ALT):
                if(i in self.population_1_indices):
                    pop_1_alt_cnt += 1
                elif(self.NumberOfPopulations() > 1 and i in self.population_2_indices):
                    pop_2_alt_cnt +=1

            if(individual_call_value.GetGenotype1() == Genotype.ALT):
                if(i in self.population_1_indices):
                    pop_1_alt_cnt += 1
                elif(self.NumberOfPopulations() > 1 and i in self.population_2_indices):
                    pop_2_alt_cnt +=1

            i+=1

        if(self.enable_debug):
            print("INDICES ARE: 1_alt: " + str(pop_1_alt_cnt) + ", 2_alt: " + str(pop_2_alt_cnt))

        self.population_1_spectrum[pop_1_alt_cnt] += 1
        if(self.NumberOfPopulations() > 1):
            self.population_2_spectrum[pop_2_alt_cnt] += 1
            self.two_populations_spectrum[pop_1_alt_cnt][pop_2_alt_cnt] += 1
        
        if(self.enable_debug):
            print("Variant call added. Spectrums: pop1: [" + str(self.population_1_spectrum).strip('[]') + "], pop2: [" + str(self.population_2_spectrum).strip('[]') + "]")


    def GenerateOutputfile(self, file_name: str):
        """
        Generates the output file in the format of the reader
        In this case, two files are generated:
        - one in dadi format (prefixed with 'dadi_')
        - one in FastSimCoal format (prefixed with ('fsc_')
        """

        self.generateOnePopulationOutputFiles(self.population_1_spectrum, file_name + "_pop1", False)
        self.generateOnePopulationOutputFiles(self.compute1PopFoldedSpectrum(self.population_1_spectrum), file_name + "_pop1", True)

        if(self.NumberOfPopulations() == 2):
            self.generateOnePopulationOutputFiles(self.population_2_spectrum, file_name + "_pop2", False)
            self.generateOnePopulationOutputFiles(self.compute1PopFoldedSpectrum(self.population_2_spectrum), file_name + "_pop2", True)
            
            self.generateTwoPopulationsOutputFiles(self.two_populations_spectrum, file_name, False)
            self.generateTwoPopulationsOutputFiles(self.compute2PopFoldedSpectrum(), file_name, True)



    def compute1PopFoldedSpectrum(self, one_population_spectrum):
        """
        Computes the folder spectrum for the provided one population spectrum.
        """

        print(datetime.now().strftime("%H:%M:%S") + ": Computing folded spectrum for one population...")

        # initialization
        folded_spectrum = [0] * len(one_population_spectrum)

        # computation of the folded spectrum:
        # first, make a sum of the unfolded spectrum with itself reversed for the left part of the spectrum
        # if the spectrum has an uneven number of elements, take the element in the middle
        # the last elements remain set to 0
        for i in range (0, math.floor(len(one_population_spectrum) / 2)):
            folded_spectrum[i] = one_population_spectrum[i] + one_population_spectrum[len(one_population_spectrum) - i - 1]

        if(math.fmod(len(one_population_spectrum), 2) != 0):
            folded_spectrum[math.floor(len(one_population_spectrum) / 2)] = one_population_spectrum[math.floor(len(one_population_spectrum) / 2)]

        print(datetime.now().strftime("%H:%M:%S") + ": Done")

        return folded_spectrum



    def compute2PopFoldedSpectrum(self):
        """
        Computes the folder spectrum for two populations.
        """

        print(datetime.now().strftime("%H:%M:%S") + ": Computing folded spectrum for two populations...")

        # initialization
        folded_spectrum = [0] * len(self.population_1_spectrum)
        for i in range (0, len(self.population_1_spectrum)):
            folded_spectrum[i] = [0] * len(self.population_2_spectrum)

        # computation of the folded spectrum:
        # first, make a sum of the unfolded spectrum with itself reversed
        # then set the elements of the left down corner to 0
        # finally, divide its values on the diagonal by 2
        for i in range (0, len(self.population_1_spectrum)):
            for j in range (0, len(self.population_2_spectrum)):
                if(i + j > (len(self.population_1_spectrum) + len(self.population_2_spectrum)) / 2 - 1):
                    folded_spectrum[i][j] = 0
                else:
                    folded_spectrum[i][j] = self.two_populations_spectrum[i][j] + self.two_populations_spectrum[len(self.population_1_spectrum) - i - 1][len(self.population_2_spectrum) - j - 1]
                if(i + j == (len(self.population_1_spectrum) + len(self.population_2_spectrum)) / 2 - 1 and folded_spectrum[i][j] != 0):
                    folded_spectrum[i][j] /= 2

        print(datetime.now().strftime("%H:%M:%S") + ": Done")

        return folded_spectrum


    def generateOnePopulationOutputFiles(self, spectrum, file_name: str, is_folded: bool):
        """
        Generates output files for one population spectrum, in dadi and FastSimCoal formats
        """

        result_dadi = str(len(spectrum)) + "\n"

        result_fsc = "1 observations\n"
        for i in range(0, len(spectrum)):
            result_fsc += "\td_" + str(i)
        result_fsc += "\n"

        for i in range(0, len(spectrum)):
            result_dadi += str(spectrum[i]) + " "
            result_fsc += str(spectrum[i]) + " "
            
        self.writeResultInFile(result_dadi, file_name, "dadi", is_folded)
        self.writeResultInFile(result_fsc, file_name, "fsc", is_folded)


    def generateTwoPopulationsOutputFiles(self, spectrum, file_name: str, is_folded: bool):
        """
        Generates output files for two populations spectrum, in dadi and FastSimCoal formats
        """

        rows = len(spectrum)
        cols = len(spectrum[0])

        result_dadi = str(rows) + " " + str(cols) + "\n"

        result_fsc = "1 observations\n"
        for j in range(0, cols):
            result_fsc += "\td_" + str(j)
        result_fsc += "\n"

        for i in range(0, rows):
            for j in range(0, cols):
                result_dadi += str(spectrum[i][j]) + " "
                if(j == 0):
                    result_fsc += "d_" + str(i)
                result_fsc += "\t" + str(spectrum[i][j])
                if(j == cols - 1):
                    result_fsc += "\n"

        self.writeResultInFile(result_dadi, file_name, "dadi", is_folded)
        self.writeResultInFile(result_fsc, file_name, "fsc", is_folded)


    def writeResultInFile(self, result: str, file_name: str, file_format: str, is_folded: bool):
        file_format_to_display = file_format
        folded_suffix = ""
        if(is_folded):
            file_format_to_display = file_format_to_display + "_folded"
            folded_suffix = "_folded"

        print(datetime.now().strftime("%H:%M:%S") + ": Generating output file in " + file_format_to_display + " format...")

        dadi_file = open(file_name + "_" + file_format + folded_suffix + ".out", "w")
        dadi_file.write(result)
        dadi_file.close()

        print(datetime.now().strftime("%H:%M:%S") + ": Done.")