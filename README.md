# pychemkin

A chemical kinetics library written in `Python`.  

Easy to extend and contribute to.

List of capabilities:
1. Elementary reactions
2. Some non-elementary reactions
  - Duplicate
  - Troe Three-Body
  - Lindmann
  - Troe Falloff
3. Irreversible and reversible rxns
4. Database of species properties:
  - Molecular weights
  - NASA polynomial coefficients
5. Time-integration of system of ODEs
  - Convenient interface with SciPy
  - User can add their own integrator
6. Contains an implementation of the energy 
   equation
  - Formulated in terms of the temperature of 
    the mixture
7. Flexible visualizations 
