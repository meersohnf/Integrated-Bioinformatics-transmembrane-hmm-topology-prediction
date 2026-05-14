import math
__author__ = 'Wombat'

# Version 2.0 18 April 2024.  Conversion to Python 3.x
# This version implements the comparison operators gt and lt, but not those involving equality. Testing for equality
# with LogFloats (as with normal floats) is a bad idea anyway owing to rounding problems and finite precision.
# The test suite is significantly expanded.
#
# The LogFloat class implements a numeric type by overloading a few of the basic math operators.
# Addition, multiplication, and division are all supported.  Subtraction is not supported, nor is the unary - operator
# In this case, the +, *, /, <, and > operators have all been overloaded in such a way that the end user of this class
# can largely use class instances as if they were normal (zero or positively valued) floats in real space.
# Note that implementation of the comparison operators are required so that comparisons may be directly made,
# but also so that Python functions such as max and sort will work as intended.
# However, values associated with class instances are actually represented internally in log space, and the
# operations are carried out entirely in log space. This helps to avoid numerical stability problems that arise from
# underflow errors in several HMM algorithms that involve successively taking the product of small probabilities.
#
# To use this class, the user need only initialize probability variables in dynamic programming matrices using
# LogFloat(value) instead of value.  Note that once a variable has been initialized as a LogFloat, successive
# operations with that variable using non-negative floats or ints as on of the operands will return a value that is
# a LogFloat, and that the integer or float will have been automatically recast into log space for the purpose of
# the operation.
#
# While the initializer accepts two arguments, the user should normally only use the first, corresponding to the
# value to be stored in the LogFloat.  The second is a flag reserved for internal use, which, if True, is indicates
# that the value passed is already log transformed.  Such usage is discouraged, but if the user wishes to create
# instances using values that are already in log rather than in real space, they may use statements of the form
# LogFloat(value, True). Note that this class assumes and uses the natural logarithm.
# The one element of this simulation of a numeric type that is non-transparent to the user is evaluation of
# LogFloat instances in some expressions.  The evaluate_real and evaluate_log methods return an associated real or
# log space value.


class LogFloat:
    """
    >>> A = LogFloat(10)   # tests setting LogFloat with non-zero positive int  -- requires __init__
    >>> print(A)
    10.000000000000002
    >>> B = LogFloat(1)
    >>> print(B)             # tests printing LogFloat -- requires __str__
    1.0
    >>> C = A + B  # tests LogFloat + LogFloat
    >>> C              # this should print it as well, but I think may use the __repr__ method, not __str__ -- confirm
    11.000000000000002
    >>> C = A + 2.1         # tests LogFloat + float -- requires __add__
    >>> C
    12.1
    >>> C = 2.1 + A         # tests float + LogFloat -- requires __radd__
    >>> print(C)
    12.1
    >>> D = A * C           # tests LogFloat * LogFloat -- requires __mul__
    >>> print(D)
    121.00000000000003
    >>> D = A * 2           # tests LogFloat * int -- requires __mul__
    >>> print(D)
    20.000000000000007
    >>> D = 2.1 * A           # tests float * LogFloat -- requires __rmul__
    >>> print(D)
    21.00000000000001
    >>> D = A * 0           # tests LogFloat * zero int -- requires __mul__
    >>> print(D)
    0.0
    >>> D = A / 2           # tests LogFloat /  int -- requires __div__
    >>> print(D)
    5.000000000000001
    >>> D = 2 / A           # tests int / LogFloat -- requires __rdiv__
    >>> print(D)
    0.19999999999999998
    >>> D = 0 / A           # tests int / LogFloat -- requires __rdiv__
    >>> print(D)
    0.0
    >>> D = A / 0           # tests LogFloat /  zero -- requires __div__
    Traceback (most recent call last):
      File "<console>", line 1, in <module>                                      # indented lines ignored for doctest
      File "C:/Users/Wombat/Documents/Python/LogFloat.py", line 193, in __div__  # indented lines ignored for doctest
        raise ZeroDivisionError                                                  # indented lines ignored for doctest
    ZeroDivisionError
    >>> D = A / 2   # tests LogFloat /  int -- requires __div__
    >>> D
    5.000000000000001
    >>> A = LogFloat(0)    # tests setting LogFloat to zero, internally it should be LOGZERO
    >>> A
    0.0
    >>> C = A + 2           # tests adding a LOGZERO LogFloat + int -- requires __add__
    >>> C
    2.0
    >>> C = 2 + A           # tests adding an int + LOGZERO LogFloat  -- requires __radd__
    >>> C
    2.0
    >>> B = LogFloat(0)    #
    >>> C = A + B           # tests adding two LOGZERO LogFloats together -- requires __add__
    >>> C
    0.0
    >>> C = A * B           # tests multiplying two LOGZERO LogFloats -- requires __mul__
    >>> C
    0.0
    >>> A > B               # tests gt comparison of two LOGZERO LogFloats -- requires __gt__
    False
    >>> A < B               # tests lt comparison of two LOGZERO LogFloats -- requires __lt__
    False

    >>> A > -1              # tests gt comparison of LOGZERO LogFloats and negative integer -- requires __gt__
    True
    >>> A < -1              # tests lt comparison of LOGZERO LogFloat and integer -1 -- requires __lt__
    False
    >>> A > -1.0            # tests gt comparison of LOGZERO LogFloats and negative float -- requires __gt__
    True
    >>> A < -1.0            # tests lt comparison of LOGZERO LogFloat and float -1 -- requires __lt__
    False
    >>> A > 0              # tests gt comparison of LOGZERO LogFloats and integer zero -- requires __gt__
    False
    >>> A < 0               # tests lt comparison of LOGZERO LogFloat and integer 0 -- requires __lt__
    False
    >>> A > 0.0              # tests gt comparison of LOGZERO LogFloats and float zero -- requires __gt__
    False
    >>> A < 0.0             # tests lt comparison of LOGZERO LogFloat and float 0 -- requires __lt__
    False
    >>> A < 1               # tests lt comparison of LOGZERO LogFloat and integer 1 -- requires __lt__
    True
    >>> A > 1               # tests lt comparison of LOGZERO LogFloat and integer 1 -- requires __gt__
    False
    >>> A < 1.0            # tests lt comparison of LOGZERO LogFloat and integer 1 -- requires __lt__
    True
    >>> A > 1.0            # tests lt comparison of LOGZERO LogFloat and integer 1 -- requires __gt__
    False

    >>> -1 > A              # tests gt comparison of LOGZERO LogFloats and negative integer -- requires __lt__
    False
    >>> -1 < A              # tests lt comparison of LOGZERO LogFloat and integer -1 -- requires __gt__
    True
    >>> -1.0 > A            # tests gt comparison of LOGZERO LogFloats and negative float -- requires __lt__
    False
    >>> -1.0 < A            # tests lt comparison of LOGZERO LogFloat and float -1 -- requires __gt__
    True
    >>> 0 > A              # tests gt comparison of LOGZERO LogFloats and integer zero -- requires __lt__
    False
    >>> 0 < A               # tests lt comparison of LOGZERO LogFloat and integer 0 -- requires __gt__
    False
    >>> 0.0 > A              # tests gt comparison of LOGZERO LogFloats and float zero -- requires __lt__
    False
    >>> 0.0 < A             # tests lt comparison of LOGZERO LogFloat and float 0 -- requires __gt__
    False
    >>> 1 < A               # tests lt comparison of LOGZERO LogFloat and integer 1 -- requires __gt__
    False
    >>> 1 > A               # tests lt comparison of LOGZERO LogFloat and integer 1 -- requires __gt__
    True
    >>> 1.0 < A            # tests lt comparison of LOGZERO LogFloat and integer 1 -- requires __gt__
    False
    >>> 1.0 > A            # tests lt comparison of LOGZERO LogFloat and integer 1 -- requires __gt__
    True

    >>> A = LogFloat(10)
    >>> A > B               # tests gt comparison of two LOGZERO LogFloats -- requires __gt__
    True
    >>> A < B               # tests lt comparison of two LOGZERO LogFloats -- requires __lt__
    False

    >>> A > -1              # tests gt comparison of LOGZERO LogFloats and negative integer -- requires __gt__
    True
    >>> A < -1              # tests lt comparison of LOGZERO LogFloat and integer -1 -- requires __lt__
    False
    >>> A > -1.0            # tests gt comparison of LOGZERO LogFloats and negative float -- requires __gt__
    True
    >>> A < -1.0            # tests lt comparison of LOGZERO LogFloat and float -1 -- requires __lt__
    False
    >>> A > 0              # tests gt comparison of LOGZERO LogFloats and integer zero -- requires __gt__
    True
    >>> A < 0               # tests lt comparison of LOGZERO LogFloat and integer 0 -- requires __lt__
    False
    >>> A > 0.0              # tests gt comparison of LOGZERO LogFloats and float zero -- requires __gt__
    True
    >>> A < 0.0             # tests lt comparison of LOGZERO LogFloat and float 0 -- requires __lt__
    False
    >>> A < 1               # tests lt comparison of LOGZERO LogFloat and integer 1 -- requires __lt__
    False
    >>> A > 1               # tests lt comparison of LOGZERO LogFloat and integer 1 -- requires __gt__
    True
    >>> A < 1.0            # tests lt comparison of LOGZERO LogFloat and integer 1 -- requires __lt__
    False
    >>> A > 1.0            # tests lt comparison of LOGZERO LogFloat and integer 1 -- requires __gt__
    True

    >>> -1 > A              # tests gt comparison of LOGZERO LogFloats and negative integer -- requires __lt__
    False
    >>> -1 < A              # tests lt comparison of LOGZERO LogFloat and integer -1 -- requires __gt__
    True
    >>> -1.0 > A            # tests gt comparison of LOGZERO LogFloats and negative float -- requires __lt__
    False
    >>> -1.0 < A            # tests lt comparison of LOGZERO LogFloat and float -1 -- requires __gt__
    True
    >>> 0 > A              # tests gt comparison of LOGZERO LogFloats and integer zero -- requires __lt__
    False
    >>> 0 < A               # tests lt comparison of LOGZERO LogFloat and integer 0 -- requires __gt__
    True
    >>> 0.0 > A              # tests gt comparison of LOGZERO LogFloats and float zero -- requires __lt__
    False
    >>> 0.0 < A             # tests lt comparison of LOGZERO LogFloat and float 0 -- requires __gt__
    True
    >>> 1 < A               # tests lt comparison of LOGZERO LogFloat and integer 1 -- requires __gt__
    True
    >>> 1 > A               # tests lt comparison of LOGZERO LogFloat and integer 1 -- requires __gt__
    False
    >>> 1.0 < A            # tests lt comparison of LOGZERO LogFloat and integer 1 -- requires __gt__
    True
    >>> 1.0 > A            # tests lt comparison of LOGZERO LogFloat and integer 1 -- requires __gt__
    False


    >>> A > B               # tests gt comparison of LogFloats -- requires __gt__
    True
    >>> A > 1               # tests gt comparison LogFloats and positive integer -- requires __gt__
    True
    >>> A > 1.0             # tests gt comparison of LogFloats and positive float -- requires __gt__
    True
    >>> A > 0               # tests gt comparison of LogFloat and integer zero -- requires __gt__
    True
    >>> A > 0.0             # tests gt comparison of LogFloats and float zero -- requires __gt__
    True
    >>> A > -1              # tests gt comparison of LogFloats -- requires __gt__
    True
    >>> A > -1.0            # tests gt comparison of LogFloats -- requires __gt__
    True
    >>> A < B               # tests lt comparison of LogFloats -- requires __lt__
    False
    >>> B < A               # tests lt comparison of LogFloats -- requires __lt__
    True
    >>> A < 2.1             # tests lt comparison of LogFloat and float -- requires __lt__
    False
    >>> A > 2.1             # tests gt comparison of LogFloat and float -- requires __gt__
    True
    >>> 2.1 < A             # tests lt comparison of float and LogFloat -- requires __lt__
    True
    >>> 2.1 > A             # tests gt comparison of float and LogFloat -- requires __gt__
    False
    >>> B = LogFloat(3.1459)
    >>> print(A, B, C, D)
    10.000000000000002 3.1459 0.0 5.000000000000001
    >>> E = [A,B,C,D]
    >>> print(E)
    [10.000000000000002, 3.1459, 0.0, 5.000000000000001]
    >>> print(max(E))
    10.000000000000002
    >>> print(sorted(E))
    [0.0, 3.1459, 5.000000000000001, 10.000000000000002]
    >>> D
    5.000000000000001
    >>> D = D + 0   #  Tests adding a zero to a log float - should be ignored
    >>> print(D)
    5.000000000000001
    """

    def __init__(self, value=None, mode_tag=None):

        if mode_tag:
            self.value = value
        else:
            if value:
                if value < 0:
                    raise ValueError('The LogFloat type supports only 0 and positively valued numbers')
                self.value = math.log(value)
            else:
                self.value = None

    def evaluate_real(self):
        if self.value is None:
            return 0
        else:
            return math.exp(self.value)

    def evaluate_log(self): return self.value

    def __gt__(self, other):

        try:
            return self.value > other.value

        except TypeError:  # either self.value or other.value (or both) are LOGZERO

            if self.value is None:
                return False

            return True  # other.value must be None, so everything is bigger than that

        except AttributeError:  # other is not a LogFloat

            if other > 0:
                if self.value is not None:
                    return self.value > math.log(other)
                else:
                    return False

            if self.value is None:
                if other < 0:
                    return True

                return False

            return True

    def __lt__(self, other):

        try:
            return self.value < other.value

        except TypeError:  # either self.value or other.value (or both) are LOGZERO

            if other.value is None:  # and self.value is None:
                return False

            return True  # self.value must be None so everything is bigger than that other than another LogFloat

        except AttributeError:  # other is not a LogFloat

            if other > 0:
                if self.value is not None:
                    return self.value < math.log(other)
                else:
                    return True

            return False

    def __le__(self, other):

        raise ValueError('Less than or equal to operator not supported. Testing for equality is not recommended.')

    def __ge__(self, other):

        raise ValueError('Greater than or equal to operator not supported. Testing for equality is not recommended.')

    # def __eq__(self, other):
    #
    #     raise ValueError('equality operator not supported. Testing for equality is not recommended.')

    def __ne__(self, other):

        raise ValueError('non-equality operator not supported. Testing for equality is not recommended.')

    def __str__(self):

        if self.value is None:
            return '0.0'
            # return "Public sees 0.0 and internally the value is None"
        else:
            return str(math.exp(self.value))
            # return str("Public sees " + str(math.exp(self.value)) + " but internally the value is " + str(self.value))

    def __repr__(self):

        if self.value is None:
            return '0.0'
            # return "Public sees 0 and internally the value is None"
        else:
            return str(math.exp(self.value))
            # return str("Public sees " + str(math.exp(self.value)) + " but internally the value is " + str(self.value))

    def __float__(self):

        if self.value is None:
            return 0.0
            # return "Public sees 0 and internally the value is None"
        else:
            return math.exp(self.value)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __mul__(self, other):

        try:
            if (self.value is None) or (other.value is None):
                return LogFloat(None, True)
            else:

                return LogFloat(self.value + other.value, True)

        except AttributeError:
            # other had no attribute named value... assume therefore that we were passed an int or a float
            if not other:
                return LogFloat(None, True)
            else:
                return LogFloat(self.value + math.log(other), True)

    def __truediv__(self, other):       # just in case division has been imported from future.

        return self.__div__(other)

    def __rtruediv__(self, other):

        return self.__rdiv__(other)

    def __rdiv__(self, other):       # Here we can assume the value on the left hand side is not a LogFloat

        if not other:
            return LogFloat(None, True)
        if self.value is None:
            raise ZeroDivisionError

        return LogFloat(math.log(other) - self.value, True)

    def __div__(self, other):

        try:
            if other.value is None:
                raise ZeroDivisionError
            if self.value is None:
                return LogFloat(None, True)

            return LogFloat(self.value - other.value, True)

        except AttributeError:
            # other had no attribute named value.. assume therefore that we were passed an int or a float
            if not other:
                raise ZeroDivisionError
            if self.value is None:
                return LogFloat(None, True)

            return LogFloat(self.value - math.log(other), True)

    def __radd__(self, other):
        return self.__add__(other)

    def __add__(self, other):

        try:
            if self.value is None or other.value is None:

                if self.value is None:
                    return LogFloat(other.value, True)
                else:
                    return LogFloat(self.value, True)

            if self.value > other.value:

                return LogFloat(self.value + math.log(1 + math.exp(other.value - self.value)), True)
                # See http://bozeman.genome.washington.edu/compbio/mbt599_2006/hmm_scaling_revised.pdf
            else:
                return LogFloat(other.value + math.log(1 + math.exp(self.value-other.value)), True)

        except AttributeError:

            if not other:
                return LogFloat(self.value, True)

            try:
                log_other = math.log(other)
            except ValueError:
                if other == 0:
                    return self
                raise ValueError('LogFloat does not support adding negative values')

            if self.value is None:
                return LogFloat(log_other, True)

            if self.value > log_other:

                return LogFloat(self.value + math.log(1 + math.exp(log_other - self.value)), True)
            else:
                return LogFloat(log_other + math.log(1 + math.exp(self.value - log_other)), True)


if __name__ == "__main__":

    import doctest
    doctest.testmod()
    # print('Doctest should have run')

    # Usage examples
    A = LogFloat(2)
    B = LogFloat(3)

    print(A)
    print(B)

    C = A * B

    print('Multiplying in log space')

    print(C)

    D = C + 2.0

    print('Adding in log space')

    print(D)
