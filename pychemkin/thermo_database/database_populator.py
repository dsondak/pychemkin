import re
from copy import deepcopy
import xml.etree.ElementTree as ET
from xml.etree.ElementTree import Element, SubElement, Comment
from xml.etree.ElementTree import XMLID
from xml.etree import ElementTree
from xml.dom import minidom
import sqlite3
import numpy as np
import pandas as pd

class Parser_7_coeffs:
    """
    Parser for the 7-degree reaction coefficients from the BURCAT Database.

    Requires following imports:
    import re
    from copy import deepcopy
    import xml.etree.ElementTree as ET
    from xml.etree.ElementTree import Element, SubElement, Comment
    from xml.etree.ElementTree import XMLID
    from xml.etree import ElementTree
    from xml.dom import minidom

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
        self.txt_file = f

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
                # get species molec weight
                if strings[-3] in ['A','B','C','D','E','F','G']:
                    specie_weight = float(strings[-2][-8:])
                else:

                    if strings[-3] in ['6000.000','5000.000']:
                        specie_weight = float(strings[-2][-8:])
                    elif strings[-3] in ['ROTATIONS','IA=9.4815']:
                        specie_weight = None
                    else:
                        specie_weight = float(strings[-2][-8:])

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
                species_info[specie,specie_state,specie_weight]={'low':{},'high':{}}
                species_info[specie,specie_state,specie_weight]['low']['Tmax'] = low_max
                species_info[specie,specie_state,specie_weight]['low']['Tmin'] = low_min
                species_info[specie,specie_state,specie_weight]['low']['coeffs'] = low_coeffs
                species_info[specie,specie_state,specie_weight]['high']['Tmax'] = high_max
                species_info[specie,specie_state,specie_weight]['high']['Tmin'] = high_min
                species_info[specie,specie_state,specie_weight]['high']['coeffs'] = high_coeffs

        return species_info

    def prettify(self, elem):
        """Return a pretty-printed XML string for the Element.
        """
        rough_string = ElementTree.tostring(elem, 'utf-8')
        reparsed = minidom.parseString(rough_string)
        return reparsed.toprettyxml(indent="  ")

    def species_dict_to_xml(self):

        species_info = self.species_txt_to_dict()

        #create root for xml file
        root = Element('specieData')
        root.set('id','specie_data')
        comment = Comment('Created by Riddhi Shah')
        root.append(comment)

        for i,k in enumerate(species_info):
            #create a specie entry in xml
            specie = Element('specie')
            specie.set('name',str(k[0]))
            specie.set('state',str(k[1]))
            specie.set('Mweight',str(k[2]))

            thermo = SubElement(specie, 'thermo')

            NASA_low = SubElement(thermo, 'NASA')
            NASA_low.set('Tmin',str(species_info[k]['low']['Tmin']))
            NASA_low.set('Tmax',str(species_info[k]['low']['Tmax']))
            floatArray = SubElement(NASA_low, 'floatArray')
            floatArray.set('name','coeffs')
            low_coeffs =", ".join(species_info[k]['low']['coeffs'])
            floatArray.text = low_coeffs

            NASA_HIGH = SubElement(thermo, 'NASA')
            NASA_HIGH.set('Tmin',str(species_info[k]['high']['Tmin']))
            NASA_HIGH.set('Tmax',str(species_info[k]['high']['Tmax']))
            floatArray = SubElement(NASA_HIGH, 'floatArray')
            floatArray.set('name','coeffs')
            high_coeffs =", ".join(species_info[k]['high']['coeffs'])
            floatArray.text = high_coeffs

            #add specie to root
            root.append(specie)

        result = self.prettify(root)
        f = open("7poly.xml","w")
        f.write(result)
        f.close()

        return result

    def create_tables(self):
        pd.set_option('display.width', 500)
        pd.set_option('display.max_columns', 100)
        pd.set_option('display.notebook_repr_html', True)

        db = sqlite3.connect('7Poly.sqlite')
        self.cursor = db.cursor()
        self.cursor.execute("DROP TABLE IF EXISTS LOW")
        self.cursor.execute("DROP TABLE IF EXISTS HIGH")
        self.cursor.execute("PRAGMA foreign_keys=1")

        #Create High and Low tables
        self.cursor.execute('''CREATE TABLE LOW (
                       SPECIES_NAME TEXT NOT NULL,
                       STATE TEXT NOT NULL,
                       MOLEC_WEIGHT TEXT NOT NULL,
                       TLOW TEXT NOT NULL,
                       THIGH TEXT NOT NULL,
                       COEFF_1 TEXT NOT NULL,
                       COEFF_2 TEXT NOT NULL,
                       COEFF_3 TEXT NOT NULL,
                       COEFF_4 TEXT NOT NULL,
                       COEFF_5 TEXT NOT NULL,
                       COEFF_6 TEXT NOT NULL,
                       COEFF_7 TEXT NOT NULL,
                       COEFF_8 TEXT NOT NULL)''')

        # Commit changes to the database
        db.commit()
        self.cursor.execute('''CREATE TABLE HIGH (
                       SPECIES_NAME TEXT NOT NULL,
                       STATE TEXT NOT NULL,
                       MOLEC_WEIGHT TEXT NOT NULL,
                       TLOW TEXT NOT NULL,
                       THIGH TEXT NOT NULL,
                       COEFF_1 TEXT NOT NULL,
                       COEFF_2 TEXT NOT NULL,
                       COEFF_3 TEXT NOT NULL,
                       COEFF_4 TEXT NOT NULL,
                       COEFF_5 TEXT NOT NULL,
                       COEFF_6 TEXT NOT NULL,
                       COEFF_7 TEXT NOT NULL)''')
        db.commit()

    def species_xml_to_db(self):
        #create xml & db
        self.species_dict_to_xml()
        self.create_tables()

        #Get the xml
        tree = ET.parse('7poly.xml')
        root = tree.getroot()

        #get species
        species = root.findall('specie')

        for specie in species:
            name = specie.get('name')
            state = specie.get('state')
            weight = specie.get('Mweight')

            #get low temp high/low and coeffs for each specie
            NASA = specie.find('thermo').findall('NASA')

            #get low info
            low_tmax = NASA[0].get('Tmax')
            low_tmin = NASA[0].get('Tmin')

            #to handle where there are 8 low coeffs
            lows = NASA[0].find('floatArray').text.split()
            if(len(lows) > 7):
                Low_C_1,Low_C_2,Low_C_3,Low_C_4,Low_C_5,Low_C_6,Low_C_7,Low_C_8 = lows[0:8]
                lows_to_insert = (name,state,weight,low_tmin,low_tmax,Low_C_1.strip(','),Low_C_2.strip(','),Low_C_3.strip(','),Low_C_4.strip(','),Low_C_5.strip(','),Low_C_6.strip(','),Low_C_7.strip(','),Low_C_8.strip(','))
            else:
                Low_C_1,Low_C_2,Low_C_3,Low_C_4,Low_C_5,Low_C_6,Low_C_7 = lows[0:7]
                lows_to_insert = (name,state,weight,low_tmin,low_tmax,Low_C_1.strip(','),Low_C_2.strip(','),Low_C_3.strip(','),Low_C_4.strip(','),Low_C_5.strip(','),Low_C_6.strip(','),Low_C_7.strip(','),"")

            #get low info
            high_tmax = NASA[1].get('Tmax')
            high_tmin = NASA[1].get('Tmin')
            High_C_1,High_C_2,High_C_3,High_C_4,High_C_5,High_C_6,High_C_7 = NASA[1].find('floatArray').text.strip(',').split()[0:7]
            high_to_insert = name,state,weight,high_tmin,high_tmax,High_C_1.strip(','),High_C_2.strip(','),High_C_3.strip(','),High_C_4.strip(','),High_C_5.strip(','),High_C_6.strip(','),High_C_7.strip(',')

            #Insert the values for each species into table
            self.cursor.execute('''INSERT INTO LOW
                          (SPECIES_NAME, STATE,MOLEC_WEIGHT, TLOW, THIGH, COEFF_1, COEFF_2,COEFF_3,COEFF_4,COEFF_5,COEFF_6,COEFF_7, COEFF_8)
                          VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)''', lows_to_insert)
            self.cursor.execute('''INSERT INTO HIGH
                          (SPECIES_NAME, STATE,MOLEC_WEIGHT, TLOW, THIGH, COEFF_1, COEFF_2,COEFF_3,COEFF_4,COEFF_5,COEFF_6,COEFF_7)
                          VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)''', high_to_insert)

    def create_sql_db(self):
          self.species_xml_to_db()
