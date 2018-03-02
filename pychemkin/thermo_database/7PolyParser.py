import re
from copy import deepcopy
import xml.etree.ElementTree as ET

class Parser_7_coeffs:
    """
    Parser for the 7-degree reaction coefficients from the BURCAT Database.

    Requires following imports:
    import re
    from copy import deepcopy
    import xml.etree.ElementTree as ET

    The default txt file to parse the 7 polynomialreactions from is defaulted in initialization:
    7poly_scrapper_output.txt which is created by Scrapers.py

    INPUTS
    =======
    The url where the txt file resides.

    RETURNS
    ========
    A database table with the columns: molecular wight, specie name, state,
    low tmin, low thigh,high tmin, high thigh, low coeffs, high coeffs
    """
    
    def __init__(self,f='7poly_scrapper_output.txt'):
        self.txt_file = file

    def species_txt_to_dict(self):
        #read the file
        file = open(self.txt_file,'r')
        lines = file.readlines()
        file.close()
        species_info = {}

        #iterate over lines
        for line in lines:
            #spilt the line
            strings = line.split()

            #skip if string is empty array
            if len(strings) <1:
                continue

            #find line with species
            if strings[-1] == '1':
                #get species name
                specie = strings[0]
                specie_state = strings[1]

                #get the low temp min and max
                low_min = strings[-4]
                low_max = 1000.000

                #get the high temp min and max
                high_min = 1000.000
                high_max = strings[-3]

            if strings[-1] == '2':
                #spilt the line by number as the numbers are not broken up by spaces in txt file
                strings = re.findall(r"[+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?", line)

                #get first 5 high coefs
                high_coeffs = []
                high_coeffs.extend(strings[0:-1])


            if strings[-1] == '3':
                #spilt the line by number as the numbers are not broken up by spaces in txt file
                strings = re.findall(r"[+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?", line)

                #first two are high coeefs
                high_coeffs.extend(strings[0:2])

                #remaining are low coeffs
                low_coeffs = []
                low_coeffs.extend(strings[2:-1])

            if strings[-1] == '4':
                #spilt the line by number as the numbers are not broken up by spaces in txt file
                strings = re.findall(r"[+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?", line)

                # get low coefs
                low_coeffs.extend(strings[0:-1])

                #Add to dictionary
                species_info[specie,specie_state]={'low':{},'high':{}}
                species_info[specie,specie_state]['low']['Tmax'] = low_max
                species_info[specie,specie_state]['low']['Tmin'] = low_min
                species_info[specie,specie_state]['low']['coeffs'] = low_coeffs
                species_info[specie,specie_state]['high']['Tmax'] = high_max
                species_info[specie,specie_state]['high']['Tmin'] = high_min
                species_info[specie,specie_state]['high']['coeffs'] = high_coeffs
        return species_info
