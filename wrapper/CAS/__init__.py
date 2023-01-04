
from .. import Maths

# Import all of the Qt classes
from .. import Qt as _Qt
from .. import Base as _Base

# Import all of the classes and functions from the C++ library
from ._CAS import *

# Now define some pure Python functions and classes that are part of
# this library...

wrap = _Base._add_wrap_function(wrap)

#enable ** operator for exbase types
ExBase.__pow__ = pow
Expression.__pow__ = pow


# Define some oft-used Symbols, so that they can be
# accessed as sire.cas.x, or "from sire.cas import r, theta, lam"
x = Symbol("x")
y = Symbol("y")
z = Symbol("z")
r = Symbol("r")
theta = Symbol("theta")
phi = Symbol("phi")
lam = Symbol("lambda")

from typing import List as _List
from typing import Union as _Union


def create_symbols(symbols: _Union[str,_List[str]], *args) -> _List[Symbol]:
    """Create symbols for each of the passed strings

       Args:
        symbols (str or list[str]):
            The list of symbols to create. This can be passed in
            as a single string, list, or several arguments.

       Returns:
        Symbol or list(Symbol):
            The list of created symbols, or a single Symbol if only
            a single name is passed.

       Examples:
            >>> (x, y, z) = create_symbols(["x", "y", "z"])

            >>> (x, y, z) = create_symbols("x", "y", "z")

            >>> (r, theta, phi) = create_symbols("r", "theta", "phi")

            >>> lam = create_symbol("lam")
    """

    if type(symbols) is not list:
        s = [symbols]
    else:
        import copy
        s = copy.copy(symbols)

    for arg in args:
        s.append(arg)

    symbols = []

    for symbol in s:
        symbols.append(Symbol(str(symbol)))

    if len(symbols) == 1:
        return symbols[0]
    else:
        return symbols
