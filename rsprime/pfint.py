# Adapted from code Copyright (c) 2010 Andrew Brown <brownan@cs.duke.edu, brownan@gmail.com>
# Copyright (c) 2013 Ryan Castellucci <code@ryanc.org>
# See LICENSE.txt for license terms

import math
# adapted from http://stackoverflow.com/a/18833870
def is_prime(n):
    if (n < 2) or (n % 2 == 0 and n > 2) or (not isinstance(n, (int, long))):
        return False
    return all(n % i for i in range(3, int(math.sqrt(n)) + 1, 2))

class PFint(int):
    "Instances of this object are elements of a prime field of order p."
    # Maps integers to PFint instances
    cache = {}
    invtable = {}

    def __new__(cls, p, n=None):
        if not is_prime(p):
            raise ValueError("Specified field order is not a prime number.")
        """
        if a value wasn't specified for 'n', return a subclass of PFint for a
        field of order p
        """
        if n is None:
            name  = 'PF%dint' % p
            bases = (PFint,)
            attrs = {'__new__': lambda cls, n: PFint(p, n)}
            return type(name, bases, attrs)

        # Check cache
        # Caching sacrifices a bit of speed for less memory usage. This way,
        # the maximum number of instances of this class at any time is limited.
        try:
            return PFint.cache[p][n]
        except KeyError:
            if n >= p or n < 0:
                raise ValueError("Field elements of PF(%d) are between 0 and %d Cannot be %s" % (p, p-1, n))

        newval = int.__new__(cls, n)
        newval.__class__ = PFint(p)
        newval.p = p
        if p not in PFint.cache:
            PFint.cache[p] = {}    
            # multiplicitive inverse table, modulo b
            PFint.invtable[p] = map(lambda x: pow(x, p-2, p), range(0, p))
            # zero doesn't have a multiplicitive inverse
            PFint.invtable[p][0] = None

        PFint.cache[p][n] = newval
        return newval

    def __add__(self, other):
        if isinstance(other, PFint):
            if self.p != other.p:
                raise ValueError("Field elements must be from the same field!")
        else:
            other = PFint(self.p, other)
        "Addition in PF(p) is normal addition modulo p"
        return PFint(self.p, (int(self) + int(other)) % self.p)
    __radd__ = __add__

    def __sub__(self, other):
        if isinstance(other, PFint):
            if self.p != other.p:
                raise ValueError("Field elements must be from the same field!")
        else:
            other = PFint(self.p, other)
        "Subtraction in PF(p) is normal subtraction modulo p"
        # Python's modulo operator handles negative numbers. If we didn't, we
        # could just add p to a before subtracting b
        return PFint(self.p, (int(self) - int(other)) % self.p)

    def __rsub__(self, other):
        if isinstance(other, PFint):
            if self.p != other.p:
                raise ValueError("Field elements must be from the same field!")
        else:
            other = PFint(self.p, other)
        # We have to reverse the argument order for rsub
        return PFint(self.p, (int(other) - int(self)) % self.p)

    def __neg__(self):
        return PFint(self.p, (self.p - int(self)) % self.p)
    
    def __mul__(self, other):
        if isinstance(other, PFint):
            if self.p != other.p:
                raise ValueError("Field elements must be from the same field!")
        else:
            other = PFint(self.p, other)
        "Multiplication in PF(p)"
        return PFint(self.p, (int(self) * int(other)) % self.p)
    __rmul__ = __mul__

    def __pow__(self, power):
        if isinstance(power, PFint):
            raise TypeError("Raising a Field element to another Field element is not defined. power must be a regular integer")
        if (power < 0):
            return PFint(self.p, pow(int(self), -power, self.p)).inverse()
        return PFint(self.p, pow(int(self), power, self.p))

    def inverse(self):
        return PFint(self.p, PFint.invtable[self.p][self])

    def __div__(self, other):
        if isinstance(other, PFint):
            if self.p != other.p:
                raise ValueError("Field elements must be from the same field!")
        else:
            other = PFint(self.p, other)
        return self * other.inverse()

    def __rdiv__(self, other):
        if isinstance(other, PFint):
            if self.p != other.p:
                raise ValueError("Field elements must be from the same field!")
        else:
            other = PFint(self.p, other)
        return self.inverse() * other

    def __repr__(self):
        return "%s(%r)" % (self.__class__.__name__, int(self))

    multiply = __mul__

# vim: sw=4 ts=4 et ai si bg=dark
