
import os

# Thermodynamics database directory
THIS_DIRECTORY = os.path.dirname(os.path.realpath(__file__))
UNITS_DIRECTORY = os.path.abspath(os.path.join(THIS_DIRECTORY, 'units.csv'))
DB_DIRECTORY = os.path.abspath(os.path.join(THIS_DIRECTORY, 'thermo_database'))

# Constant(s)
R = 8.3144598 # Ideal gas constant, in J/mol/K (https://physics.nist.gov/cgi-bin/cuu/Value?r)
