import requests
from bs4 import BeautifulSoup
import pytest
import sys

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

if __name__ == "__main__":
    import doctest
    from bs4 import BeautifulSoup
    import requests
    import numpy as np
    doctest.testmod(verbose=True)
