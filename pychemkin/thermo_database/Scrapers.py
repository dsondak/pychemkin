import requests
from bs4 import BeautifulSoup
import pytest
import sys
import numpy as np
import time
import json


class Scraper_7_coeffs:
    """
    Scraper for the 7-degree reaction coefficients from the BURCAT Database.

    Requires following imports:
    from bs4 import BeautifulSoup
    import requests

    The default url to scrape from is defaulted in initialization:
    http://garfield.chem.elte.hu/Burcat/BURCAT.THR

    INPUTS
    =======
    The url where the database file resides.

    RETURNS
    ========
    String with the database table, to be passed to the parser to create XML or SQL object

    EXAMPLES
    =========
    Showing of the main scraper functionality. Returns true if everything runs properly.

    >>> test_scraper = Scraper_7_coeffs(url='http://garfield.chem.elte.hu/Burcat/BURCAT.THR')
    >>> table_of_coeffs = test_scraper.scrape_coeffs()
    >>> table_of_coeffs.__class__.__name__=='str'
    True
    """

    def __init__(self,url='http://garfield.chem.elte.hu/Burcat/BURCAT.THR'):

        self.url = url
        self.result = None
        self.text_body = None
        self.webpage = None
        self.soup = None
        self.molec_weights = None
        self.range_of_weights = None

    def scrape_url(self):
        self.web_page = requests.get(self.url)

    def create_soup_object(self):
        self.soup = BeautifulSoup(self.web_page.content, 'lxml')

    def scrape_coeffs(self):
        self.scrape_url()
        self.create_soup_object()

        my_tag = self.soup.find('egil.jahnsen')

        #Identify the tag with the 7-coeffs polynomial
        for item in my_tag.contents:
            if item.__class__.__name__ == "Tag":
                text_body = item

        # Find the exact element with the body text
        for child in text_body.contents:
            if child.__class__.__name__ == "Tag":
                raw_coeffs_body = child.get_text()

        #Separate the table
        for idx,line in enumerate(raw_coeffs_body.split('\n')):
            if 'THE NUMBER PRECEDING EACH SPECIES IS THE CHEMICAL ABSTRACT' in line:
                self.result = raw_coeffs_body.split('THE NUMBER PRECEDING EACH SPECIES IS THE CHEMICAL ABSTRACT (CAS) IDENTIFICATION.')[1]
        f = open("7poly_scrapper_output.txt","w")
        f.write(self.result)
        f.close()
        return self.result

        '''
    def scrape_species_molecular_weights(self,range_of_weights,query_lag,save_txt=False):


        self.range_of_weights = range_of_weights
        self.query_lag = query_lag
        self. save_txt = save_txt
        """Returns a dictionary of molecular weights with species as key and molecular weight as a value"""
        
        ########################
        # Helper function to get body of the table from NIST website
        def get_species_body(soup):
            """Returns species body based on NIST search query url"""
            def is_the_only_string_within_a_tag(s):
                """Return True if this string is the only child of its parent tag."""
                return (s == s.parent.string)

            molec_weights_table = soup.find_all(string=is_the_only_string_within_a_tag)
            body_text = " ".join(molec_weights_table)

            # Extract species weights
            lower_half = body_text.split("Click on the name to see more data.")[1]

            species_weights = lower_half.split('2017 by the U.S. Secretary of Commerce')[0]
            species_weights_clean = species_weights.replace('\xa0','').replace('\r','').replace("\n",'')
            return species_weights_clean
        ####################
        ### Scrape the web, utilizing function from above
        start_url = "http://webbook.nist.gov/cgi/cbook.cgi?Value="
        end_url = "&VType=MW&Formula=&AllowExtra=on&Units=SI"
        start_val = 0
        end_val = self.range_of_weights
        molec_weights_dict = {}
        cnt = 0
        
        while start_val < end_val:
        
            url = start_url + str(start_val) + '-' + str(start_val+2) + end_url
            #print(url)
            web_page = requests.get(url)
            soup = BeautifulSoup(web_page.content, 'lxml')
            species_weights_clean = get_species_body(soup)
            current_specie = []

            for member in species_weights_clean.split():
                try:
                    float(member)
                    current_weight = float(member)
                    if current_specie:
                        # Make sure we concatenate the strings if more elements before next element
                        current_specie = " ".join(current_specie)
                        molec_weights_dict[current_specie] = current_weight
                        current_specie = []
                except ValueError:
                    current_specie.append(member.upper())

            # Delete after each website scraped
            current_specie = []
            
            cnt+=1
            #add a check to break the while loop to prevent hitting the website with infinite queries
            if cnt > end_val:
                break
            start_val +=2
            time.sleep(self.query_lag)
            
        if self.save_txt:
            
            with open('Molecular_weights.txt', 'w') as file:
                 file.write(json.dumps(molec_weights_dict)) 
        
        self.molec_weights = molec_weights_dict
        return molec_weights_dict
    '''


if __name__ == "__main__":
    import doctest
    from bs4 import BeautifulSoup
    import requests
    import numpy as np
    import sys
    import time
    import json
    doctest.testmod(verbose=True)
