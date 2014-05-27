Reed-Solomon codec over prime field written in pure Python
==========================================================

Copyright (c) 2010 Andrew Brown <brownan@gmail.com> <brownan@cs.duke.edu>
Copyright (c) 2014 Ryan Castellucci <code@ryanc.org>

The original pure Python Reed-Solomon codec (operating over GF(2^8)) was
written by Andrew Brown (@brownan). It was adapted by Ryan Castellucci
(@ryancdotorg) to operate over an arbitrary prime field with debugging help
from Dennis McKinnon (@dennismckinnon).

If you're interested in the original, please see
https://github.com/brownan/Reed-Solomon

The included tests and example code run over base59. Any base should work,
though for one larger than 59 you'll need to supply a custom symbol mapper.

The example program (`rsprime/__main__.py`) can be run as `python rsprime`

Documentation (**WARNING: NOT CURRENT!**)
-------------------------------------
rs.RSCoder(n, k)
     Creates a new Reed-Solomon Encoder/Decoder object configured with
     the given n and k values.
     n is the length of a codeword, must be less than 256
     k is the length of the message, must be less than n
     
     The code will have error correcting power s where 2s = n - k
     
     The typical RSCoder is RSCoder(255, 223)
 
RSCoder Objects

RSCoder.encode(message, poly=False)
    Encode a given string with reed-solomon encoding. Returns a byte
    string with the k message bytes and n-k parity bytes at the end.
    
    If a message is < k bytes long, it is assumed to be padded at the front
    with null bytes.
    
    The sequence returned is always n bytes long.
    
    If poly is not False, returns the encoded Polynomial object instead of
    the polynomial translated back to a string (useful for debugging)
    
RSCoder.decode(r, nostrip=False)
    Given a received string or byte array r, attempts to decode it. If
    it's a valid codeword, or if there are no more than (n-k)/2 errors, the
    message is returned.
    
    A message always has k bytes, if a message contained less it is left
    padded with null bytes. When decoded, these leading null bytes are
    stripped, but that can cause problems if decoding binary data. When
    nostrip is True, messages returned are always k bytes long. This is
    useful to make sure no data is lost when decoding binary data.

RSCoder.verify(code)
    Verifies the code is valid by testing that the code as a polynomial
    code divides g
    returns True/False


Besides the main RSCoder object, two other objects are used in this
implementation. Their use is not specifically tied to the coder.

polynomial.Polynomial(coefficients=(), \**sparse)
    There are three ways to initialize a Polynomial object.
    1) With a list, tuple, or other iterable, creates a polynomial using
    the items as coefficients in order of decreasing power

    2) With keyword arguments such as for example x3=5, sets the
    coefficient of x^3 to be 5

    3) With no arguments, creates an empty polynomial, equivalent to
    Polynomial((0,))

    >>> print Polynomial((5, 0, 0, 0, 0, 0))
    5x^5

    >>> print Polynomial(x32=5, x64=8)
    8x^64 + 5x^32

    >>> print Polynomial(x5=5, x9=4, x0=2) 
    4x^9 + 5x^5 + 2

Polynomial objects export the following standard functions that perform the
expected operations using polynomial arithmetic. Arithmetic of the coefficients
is determined by the type passed in, so integers or GF256int objects could be
used, the Polynomial class is agnostic to the type of the coefficients.

::

    __add__
    __divmod__
    __eq__
    __floordiv__
    __hash__
    __len__
    __mod__
    __mul__
    __ne__
    __neg__
    __sub__
    evaluate(x)
    degree()
        Returns the degree of the polynomial
    get_coefficient(degree)
        Returns the coefficient of the specified term

ff.GF256int(value)
    Instances of this object are elements of the field GF(2^8)
    Instances are integers in the range 0 to 255
    This field is defined using the irreducable polynomial
    x^8 + x^4 + x^3 + x + 1
    and using 3 as the generator for the exponent table and log table.

The GF256int class inherits from int and supports all the usual integer
operations. The following methods are overridden for arithmetic in the finite
field GF(2^8)

::

    __add__
    __div__
    __mul__
    __neg__
    __pow__
    __radd__
    __rdiv__
    __rmul__
    __rsub__
    __sub__
    inverse()
        Multiplicative inverse in GF(2^8)


Examples
--------
>>> import rs
>>> coder = rs.RSCoder(20,13)
>>> c = coder.encode("Hello, world!")
>>> print repr(c)
'Hello, world!\x8d\x13\xf4\xf9C\x10\xe5'
>>>
>>> r = "\0"*3 + c[3:]
>>> print repr(r)
'\x00\x00\x00lo, world!\x8d\x13\xf4\xf9C\x10\xe5'
>>>
>>> coder.decode(r)
'Hello, world!'
from the full range 0-255. Also note that either the data area or the parity
