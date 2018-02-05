# pychemkin

[![Build Status](https://travis-ci.org/dsondak/pychemkin.svg?branch=develop)](https://travis-ci.org/dsondak/pychemkin.svg?branch=develop)
[![Coverage Status](https://coveralls.io/repos/github/dsondak/pychemkin/badge.svg?branch=develop)](https://coveralls.io/github/dsondak/pychemkin?branch=develop)

A chemical kinetics library written in `Python`.  

Easy to extend and contribute to.

## List of capabilities:
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
  - Formulated in terms of the temperature of the mixture
7. Flexible visualizations 
8. Consistency checks on units
9. Reactions specified through .xml input file.
10. Derivatives wrt state.
  - Try to use automatic differentiation packages?

## Milestones
* March 1
  - Elementary reactions
  - Reversible and irreversible reactions
  - Consistency checks on units
  - Input parsers
  - Time integration module
  - Visualizations
    - Time plots of various quantities
    - Graphical viz?
  - First pass at database
    - NASA polynomials for [Burcat's database](http://garfield.chem.elte.hu/Burcat/burcat.html)
    - Molecular weights
* April 1
  - Energy equation
  - Non-elementary reactions
  - Host database on external site?
  - Web API for database?
* May 1
  - Accessing database
  - Contributing to database
  - Requests?
  - Derivative information
