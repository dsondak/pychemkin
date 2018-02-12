
"""Tests for SQLParser"""

import pytest
import os

import sqlite3

from pychemkin.config import DB_DIRECTORY
TEST_DB_PATH = os.path.join(DB_DIRECTORY + "/NASA7_coeffs.sqlite")

from pychemkin.parsers.XMLParser import XMLParser
from pychemkin.parsers.SQLParser import SQLParser


def test_SQLParser_file_not_found():
    xml_filename = "tests/test_xml_files/rxn.xml"
    xml_parser = XMLParser(xml_filename)
    species = xml_parser.get_species() 
    with pytest.raises(OSError):
        sql_parser = SQLParser("no_such_file", species)

