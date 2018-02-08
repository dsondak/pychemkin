
"""Class for parsing SQLite database containing thermodynamic data"""

import os
import sqlite3

from pychemkin.config import DB_DIRECTORY


class SQLParser:
    """Class for parsing thermodynamic data from sqlite database."""
    def __init__ (self, database):
        self.database = os.path.join(DB_DIRECTORY, database)
        if not os.path.isfile(self.database):
            raise ValueError('Database {} does not exist'.format(self.database))
        self.data = None

    def get_coeffs(self, species_name, temp_range):
        """Returns NASA polynomial coefficients from sql database
        from appropriate temperature range.

        INPUTS:
        -------
        species_name : str
            name of chemical species
        temp_range : str
            temperature range for reaction (options: 'low' or 'high')

        RETURNS:
        --------
        nasa_coeffs : numpy.ndarray
            nasa coefficients for species in reaction
        """
        db = sqlite3.connect(DATA_DIRECTORY + '/' + 'NASA_poly_coeffs.sqlite')
        cursor = db.cursor()

        if temp_range not in ['high', 'low']:
            raise ValueError("Temperature range can only be 'high' or 'low'...")

        if temp_range == "high":
            cursor.execute('''SELECT COEFF_1, COEFF_2, COEFF_3, COEFF_4, 
                           COEFF_5, COEFF_6, COEFF_7 FROM HIGH WHERE 
                           SPECIES_NAME = ?''', (species_name,))
        
        else:
            cursor.execute('''SELECT COEFF_1, COEFF_2, COEFF_3, COEFF_4, 
                           COEFF_5,COEFF_6,COEFF_7 FROM LOW 
                           WHERE SPECIES_NAME = ?''', (species_name,))

        nasa_coeffs = numpy.array(cursor.fetchone())
        db.commit()
        db.close()
        return nasa_coeffs

    def get_Tmid(self, species_name):
        """Returns middle temperature (T) value.

        INPUTS:
        -------
        species_name : str
            name of chemical species

        RETURNS:
        --------
        T_mid : float
            mid temperature value (ie separates
            low and high temperature ranges)
        """
        db = sqlite3.connect(DATA_DIRECTORY + '/' + 'NASA_poly_coeffs.sqlite')
        cursor = db.cursor()
        cursor.execute('''SELECT TLOW FROM HIGH WHERE SPECIES_NAME = ?''', (species_name,))
        Tmid = cursor.fetchone()[0]
        db.commit()
        db.close()
        return Tmid

    def get_NASA_poly_coefs(self):
        """Fetches and updates all NASA polynomial coefficients (high and low temperature)
        for all species.

        RETURNS:
        --------
        NASA_poly_coeffs: dict
            Dictionary of the form {species name : dict of NASA polynomial coefficients}
            where the inner dictionary has keys 'high' (with value, array of coeffs in high T),
            'low' (with value, array of coeffs in low T), and 'Tmid' (with value, float of middle T)
        """
        NASA_poly_coefs = {}
        for species_name in self.species:
            coef = {}
            coef['Tmid'] = self.get_Tmid(species_name)
            coef['low'] = self.get_coeffs(species_name, 'low')
            coef['high'] = self.get_coeffs(species_name, 'high')
            NASA_poly_coefs[species_name] = coef
        self.NASA_poly_coefs = NASA_poly_coefs
        return self.NASA_poly_coefs



