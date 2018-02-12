
"""Class for parsing SQLite database containing thermodynamic data"""

import os
import sqlite3
import numpy

from pychemkin.config import DB_DIRECTORY


class SQLParser:
    """Class for parsing thermodynamic data from sqlite database."""
    def __init__ (self, database, species):
        self.database = os.path.join(DB_DIRECTORY, database)
        if not os.path.isfile(self.database):
            raise OSError('Database {} does not exist'.format(self.database))
        self.species = species
        self.thermo_coefficients = {}

    def get_Tmid(self, db_cursor, species_name):
        """Helper function that returns middle temperature (T) value.

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
        db_cursor.execute('''SELECT TLOW FROM HIGH WHERE SPECIES_NAME = ?''', (species_name,))
        Tmid = db_cursor.fetchone()[0]
        return Tmid

    def get_coeffs_from_T_range(self, db_cursor, species_name, temp_range):
        """Helper function that returns 7th order NASA polynomial coefficients from sql database
        from appropriate temperature range.

        TODO: Generalize this?

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
        if temp_range == "high":
            db_cursor.execute('''SELECT COEFF_1, COEFF_2, COEFF_3, COEFF_4, 
                           COEFF_5, COEFF_6, COEFF_7 FROM HIGH WHERE 
                           SPECIES_NAME = ?''', (species_name,))
        else:
            db_cursor.execute('''SELECT COEFF_1, COEFF_2, COEFF_3, COEFF_4, 
                           COEFF_5,COEFF_6,COEFF_7 FROM LOW 
                           WHERE SPECIES_NAME = ?''', (species_name,))

        nasa_coeffs = numpy.array(db_cursor.fetchone())
        return nasa_coeffs

    def get_thermo_coeffs(self):
        """Fetches and updates all NASA polynomial coefficients (high and low temperature)
        for all species.

        RETURNS:
        --------
        thermo_coefficients : dict
            Dictionary of the form {species name : dict of NASA polynomial coefficients}
            where the inner dictionary has keys 'high' (with value, array of coeffs in high T),
            'low' (with value, array of coeffs in low T), and 'Tmid' (with value, float of middle T)
        """
        db = sqlite3.connect(self.database)
        cursor = db.cursor()
        for species_name in self.species:
            coeff = {}
            coeff['Tmid'] = self.get_Tmid(db_cursor=cursor,
                                         species_name=species_name)
            coeff['low'] = self.get_coeffs_from_T_range(db_cursor=cursor,
                                                       species_name=species_name,
                                                       temp_range='low')
            coeff['high'] = self.get_coeffs_from_T_range(db_cursor=cursor,
                                                        species_name=species_name,
                                                        temp_range='high')
            self.thermo_coefficients[species_name] = coeff
        db.commit()
        db.close()
        return self.thermo_coefficients
