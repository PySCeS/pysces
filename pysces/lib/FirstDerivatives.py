# Automatic first-order derivatives
#
# Written by Konrad Hinsen <hinsen@cnrs-orleans.fr>
# last revision: 2002-6-13
#

# Adapted slightly to work independently and use numpy - bgoli
# Brett Olivier 20070308

"""This module provides automatic differentiation for functions with
any number of variables. Instances of the class DerivVar represent the
values of a function and its partial derivatives with respect to a
list of variables. All common mathematical operations and functions
are available for these numbers.  There is no restriction on the type
of the numbers fed into the code; it works for real and complex
numbers as well as for any Python type that implements the necessary
operations.

This module is as far as possible compatible with the n-th order
derivatives module Derivatives. If only first-order derivatives
are required, this module is faster than the general one.

Example:

  >>>print sin(DerivVar(2))

  produces the output

  >>>(0.909297426826, [-0.416146836547])

The first number is the value of sin(2); the number in the following
list is the value of the derivative of sin(x) at x=2, i.e. cos(2).

When there is more than one variable, DerivVar must be called with
an integer second argument that specifies the number of the variable.

Example:

  >>>x = DerivVar(7., 0)
  >>>y = DerivVar(42., 1)
  >>>z = DerivVar(pi, 2)
  >>>print (sqrt(pow(x,2)+pow(y,2)+pow(z,2)))

  produces the output

  >>>(42.6950770511, [0.163953328662, 0.98371997197, 0.0735820818365])

The numbers in the list are the partial derivatives with respect
to x, y, and z, respectively.

Note: It doesn't make sense to use DerivVar with different values
for the same variable index in one calculation, but there is
no check for this. I.e.

  >>>print DerivVar(3, 0)+DerivVar(5, 0)

  produces

  >>>(8, [2])

but this result is meaningless.
"""
from __future__ import division, print_function
from __future__ import absolute_import
from __future__ import unicode_literals

# Changed Numeric to Numpy and updated code to reflect this - bgoli

# import Numeric
try:
    import numpy
except Exception as ex:
    print(ex)
    print('numpy import failed trying to import numeric (this is not ideal) ... ')
    import Numeric as numpy

# The following class represents variables with derivatives:

class DerivVar:

    """Variable with derivatives

    Constructor: DerivVar(|value|, |index| = 0)

    Arguments:

    |value| -- the numerical value of the variable

    |index| -- the variable index (an integer), which serves to
               distinguish between variables and as an index for
               the derivative lists. Each explicitly created
               instance of DerivVar must have a unique index.

    Indexing with an integer yields the derivatives of the corresponding
    order.
    """

    def __init__(self, value, index=0, order=1):
        if order > 1:
            raise ValueError('Only first-order derivatives')
        self.value = value
        if order == 0:
            self.deriv = []
        elif type(index) == type([]):
            self.deriv = index
        else:
            self.deriv = index*[0] + [1]

    def __getitem__(self, item):
        if item < 0 or item > 1:
            raise ValueError('Index out of range')
        if item == 0:
            return self.value
        else:
            return self.deriv

    def __repr__(self):
        return repr((self.value, self.deriv))

    def __str__(self):
        return str((self.value, self.deriv))

    def __coerce__(self, other):
        if isDerivVar(other):
            return self, other
        else:
            return self, DerivVar(other, [])

    def __cmp__(self, other):
        return cmp(self.value, other.value)

    def __neg__(self):
        return DerivVar(-self.value,[-a for a in self.deriv])

    def __pos__(self):
        return self

    def __abs__(self): # cf maple signum # derivate of abs
        absvalue = abs(self.value)
        return DerivVar(absvalue, list(map(lambda a, d=self.value/absvalue:
                                      d*a, self.deriv)))
    def __bool__(self):
        return self.value != 0

    def __add__(self, other):
        return DerivVar(self.value + other.value,
                _mapderiv(lambda a,b: a+b, self.deriv, other.deriv))

    __radd__ = __add__

    def __sub__(self, other):
        return DerivVar(self.value - other.value,
            _mapderiv(lambda a,b: a-b, self.deriv, other.deriv))

    def __rsub__(self, other):
        return DerivVar(other.value - self.value,
            _mapderiv(lambda a,b: a-b, other.deriv, self.deriv))

    def __mul__(self, other):
        return DerivVar(self.value*other.value,
            _mapderiv(lambda a,b: a+b,
                  list(map(lambda x,f=other.value:f*x, self.deriv)),
                  list(map(lambda x,f=self.value:f*x, other.deriv))))

    __rmul__ = __mul__

    def __div__(self, other):
        if not other.value:
            raise ZeroDivisionError('DerivVar division')
        inv = 1./other.value
        return DerivVar(self.value*inv,
                _mapderiv(lambda a,b: a-b,
                      list(map(lambda x,f=inv: f*x, self.deriv)),
                      list(map(lambda x,f=self.value*inv*inv: f*x,
                          other.deriv))))

    def __rdiv__(self, other):
        return other/self

    def __pow__(self, other, z=None):
        if z is not None:
            raise TypeError('DerivVar does not support ternary pow()')
        val1 = pow(self.value, other.value-1)
        val = val1*self.value
        deriv1 = list(map(lambda x,f=val1*other.value: f*x, self.deriv))
        if isDerivVar(other) and len(other.deriv) > 0:
            deriv2 = list(map(lambda x, f=val*numpy.log(self.value): f*x,
                             other.deriv))
            return DerivVar(val,_mapderiv(lambda a,b: a+b, deriv1, deriv2))
        else:
            return DerivVar(val,deriv1)

    def __rpow__(self, other):
        return pow(other, self)

    def exp(self):
        v = numpy.exp(self.value)
        return DerivVar(v, list(map(lambda x,f=v: f*x, self.deriv)))

    def log(self):
        v = numpy.log(self.value)
        d = 1./self.value
        return DerivVar(v, list(map(lambda x,f=d: f*x, self.deriv)))

    def log10(self):
        v = numpy.log10(self.value)
        d = 1./(self.value * numpy.log(10))
        return DerivVar(v, list(map(lambda x,f=d: f*x, self.deriv)))

    def sqrt(self):
        v = numpy.sqrt(self.value)
        d = 0.5/v
        return DerivVar(v, list(map(lambda x,f=d: f*x, self.deriv)))

    def sign(self):
        if self.value == 0:
            raise ValueError("can't differentiate sign() at zero")
        return DerivVar(numpy.sign(self.value), 0)

    def sin(self):
        v = numpy.sin(self.value)
        d = numpy.cos(self.value)
        return DerivVar(v, list(map(lambda x,f=d: f*x, self.deriv)))

    def cos(self):
        v = numpy.cos(self.value)
        d = -numpy.sin(self.value)
        return DerivVar(v, list(map(lambda x,f=d: f*x, self.deriv)))

    def tan(self):
        v = numpy.tan(self.value)
        d = 1.+pow(v,2)
        return DerivVar(v, list(map(lambda x,f=d: f*x, self.deriv)))

    def sinh(self):
        v = numpy.sinh(self.value)
        d = numpy.cosh(self.value)
        return DerivVar(v, list(map(lambda x,f=d: f*x, self.deriv)))

    def cosh(self):
        v = numpy.cosh(self.value)
        d = numpy.sinh(self.value)
        return DerivVar(v, list(map(lambda x,f=d: f*x, self.deriv)))

    def tanh(self):
        v = numpy.tanh(self.value)
        d = 1./pow(numpy.cosh(self.value),2)
        return DerivVar(v, list(map(lambda x,f=d: f*x, self.deriv)))

    def arcsin(self):
        v = numpy.arcsin(self.value)
        d = 1./numpy.sqrt(1.-pow(self.value,2))
        return DerivVar(v, list(map(lambda x,f=d: f*x, self.deriv)))

    def arccos(self):
        v = numpy.arccos(self.value)
        d = -1./numpy.sqrt(1.-pow(self.value,2))
        return DerivVar(v, list(map(lambda x,f=d: f*x, self.deriv)))

    def arctan(self):
        v = numpy.arctan(self.value)
        d = 1./(1.+pow(self.value,2))
        return DerivVar(v, list(map(lambda x,f=d: f*x, self.deriv)))

    def arctan2(self, other):
        den = self.value*self.value+other.value*other.value
        s = self.value/den
        o = other.value/den
        return DerivVar(numpy.arctan2(self.value, other.value),
                _mapderiv(lambda a,b: a-b,
                      list(map(lambda x,f=o: f*x, self.deriv)),
                      list(map(lambda x,f=s: f*x, other.deriv))))

# Don't know if this will work - bgoli
    ##  def gamma(self):
        ##  from transcendental import gamma, psi
        ##  v = gamma(self.value)
        ##  d = v*psi(self.value)
        ##  return DerivVar(v, map(lambda x,f=d: f*x, self.deriv))

# Type check

def isDerivVar(x):
    "Returns 1 if |x| is a DerivVar object."
    return hasattr(x,'value') and hasattr(x,'deriv')

# Map a binary function on two first derivative lists

def _mapderiv(func, a, b):
    nvars = max(len(a), len(b))
    a = a + (nvars-len(a))*[0]
    b = b + (nvars-len(b))*[0]
    return list(map(func, a, b))


# Don't use this - bgoli
##  def DerivVector(x, y, z, index=0):

    ##  """Returns a vector whose components are DerivVar objects.

    ##  Arguments:

    ##  |x|, |y|, |z| -- vector components (numbers)

    ##  |index| -- the DerivVar index for the x component. The y and z
               ##  components receive consecutive indices.
    ##  """

    ##  from Scientific.Geometry.Vector import Vector
    ##  if isDerivVar(x) and isDerivVar(y) and isDerivVar(z):
        ##  return Vector(x, y, z)
    ##  else:
        ##  return Vector(DerivVar(x, index),
                      ##  DerivVar(y, index+1),
                      ##  DerivVar(z, index+2))
