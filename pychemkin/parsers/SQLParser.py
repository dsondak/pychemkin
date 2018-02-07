
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

    def get_coeffs(self, species_name, temperature_range):
        """Returns NASA polynomial coefficients of a specie
        within a given temperature range.

        Args:
        =====
        species_name : str, required
            name of specie
        temperature_range : float, required
            temperature of reaction

        Returns:
        ========
        coeffs : list
            list containing NASA polynomial coefficients
        """
        db = sqlite3.connect(self.db_path)
        cursor = db.cursor()
        query = '''SELECT COEFF_1, COEFF_2, COEFF_3, COEFF_4, COEFF_5, COEFF_6, COEFF_7
                    FROM {} 
                    WHERE SPECIES_NAME = "{}"'''.format(temp_range.upper(), species_name)
        coeffs = list(cursor.execute(query).fetchall()[0])
        db.commit()
        db.close()
        return coeffs

    def get_species(self, temp, temp_range):
        db = sqlite3.connect(self.db_path)
        cursor = db.cursor()
        if temp_range == 'low': # temp_range == 'low'
            query = '''SELECT SPECIES_NAME FROM {} WHERE TLOW < {}'''.format(temp_range.upper(), temp)
        else: # temp_range == 'high'
            query = '''SELECT SPECIES_NAME FROM {} WHERE THIGH > {}'''.format(temp_range.upper(), temp)
        species = []
        for s in cursor.execute(query).fetchall():
            species.append(s[0])
        db.commit()
        db.close()
        return species

