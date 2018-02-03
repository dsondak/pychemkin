
"""Class for PyChemkin-related errors."""


class PyChemKinError(Exception):
    """Encapsulates errors encountered in this package."""
    def __init__(self, method, info=None):
        """Initializes PyChemkin-related error.

        Args:
        -----
        method : str
            name of method that is causing error
        info : str, optional
            further information about particular error
        """
        msg = 'Error encountered in chemkin.py method: {0}.'.format(method)
        if info is not None:
            msg = msg + ' ' + info
        Exception.__init__(self, msg)

        self.method = method
        self.info = info
