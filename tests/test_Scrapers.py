
#Import the desired module
from pychemkin.thermo_database.Scrapers import *

class Test_Scrapers():
    """
    Class for testing the parsers basic functionality.
    """
    
    #def __init__(self):
    #    self.success = True
    #    self.fail = False
    
    @staticmethod
    def test_Scraper_7():
        """
        Tests for Parser of 7-coeff NASA Polynomials
        """
        test7 = Scraper_7_coeffs()
        # The following method scrape_coeffs() implicitly calls scrape_url() and create_soup_object()
        table = test7.scrape_coeffs()
        
        #Test that the server response to the parser was 200.
        assert str(test7.web_page.__repr__) == "<bound method Response.__repr__ of <Response [200]>>"
        
        # Test soup object creation
        assert str(test7.soup.__class__.__name__) == "BeautifulSoup"
        
        # Test isolation of table of interest by testing type and length
        
        assert type(table) == str and len(table)>1400000
              