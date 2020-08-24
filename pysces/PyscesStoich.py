"""
PySCeS - Python Simulator for Cellular Systems (http://pysces.sourceforge.net)

Copyright (C) 2004-2020 B.G. Olivier, J.M. Rohwer, J.-H.S Hofmeyr all rights reserved,

Brett G. Olivier (bgoli@users.sourceforge.net)
Triple-J Group for Molecular Cell Physiology
Stellenbosch University, South Africa.

Permission to use, modify, and distribute this software is given under the
terms of the PySceS (BSD style) license. See LICENSE.txt that came with
this distribution for specifics.

NO WARRANTY IS EXPRESSED OR IMPLIED.  USE AT YOUR OWN RISK.
Brett G. Olivier
"""
from __future__ import division, print_function
from __future__ import absolute_import
from __future__ import unicode_literals

from pysces.version import __version__
import time, os
import scipy
import scipy.linalg

# if int(scipy.__version__.split('.')[1]) < 12:
# thanks scipy
if (
    int(scipy.__version__.split('.')[0]) < 1
    and int(scipy.__version__.split('.')[1]) < 12
):
    myfblas = scipy.linalg.fblas
    myflapack = scipy.linalg.flapack
else:
    myfblas = scipy.linalg.blas
    myflapack = scipy.linalg.lapack

import copy

__doc__ = """
          PyscesStoich
          ------------
          PySCeS stoichiometric analysis classes.

          """

##  print 'Stoichiometry ver ' + __version__ + ' runtime: '+ time.strftime("%H:%M:%S")

##stoich_zero_valM = scipy.machar.MachAr().eps*2.0e4 # safer --> about 4e-11 (unofficially - a minpivot size)
##stoich_zero_valM = scipy.machar.MachAr().eps*2.0e4 # safe --> about 4e-12 (unofficially - a minpivot size)
##print '\tStoichiometric precision = ', stoich_zero_valM


class StructMatrix:
    """
    This class is specifically designed to store structural matrix information
    give it an array and row/col index permutations it can generate its own
    row/col labels given the label src.
    """

    array = None
    ridx = None
    cidx = None
    row = None
    col = None
    shape = None

    def __init__(self, array, ridx, cidx, row=None, col=None):
        """
        Instantiate with array and matching row/col index arrays, optional label arrays
        """
        self.array = array
        self.ridx = ridx
        self.cidx = cidx
        self.row = row
        self.col = col
        self.shape = array.shape

    def __call__(self):
        return self.array

    def getRowsByIdx(self, *args):
        """Return the rows referenced by index (1,3,5)"""
        return self.array.take(args, axis=0)

    def getColsByIdx(self, *args):
        """Return the columns referenced by index (1,3,5)"""
        return self.array.take(args, axis=1)

    def setRow(self, src):
        """
        Assuming that the row index array is a permutation (full/subset)
        of a source label array by supplying that source to setRow it
        maps the row labels to ridx and creates self.row (row label list)
        """
        self.row = [src[r] for r in self.ridx]

    def setCol(self, src):
        """
        Assuming that the col index array is a permutation (full/subset)
        of a source label array by supplying that src to setCol
        maps the row labels to cidx and creates self.col (col label list)
        """

        self.col = [src[c] for c in self.cidx]

    def getRowsByName(self, *args):
        """Return the rows referenced by label ('s','x','d')"""
        assert self.row != None, "\nI need row labels"
        try:
            return self.array.take([self.row.index(l) for l in args], axis=0)
        except Exception as ex:
            print(ex)
            print("\nValid row labels are: %s" % self.row)
            return None

    def getColsByName(self, *args):
        """Return the columns referenced by label ('s','x','d')"""
        assert self.col != None, "\nI need column labels"
        try:
            return self.array.take([self.col.index(l) for l in args], axis=1)
        except Exception as ex:
            print(ex)
            print("Valid column labels are: %s" % self.col)
            return None

    def getLabels(self, axis='all'):
        """Return the matrix labels ([rows],[cols]) where axis='row'/'col'/'all'"""
        if axis == 'row':
            return self.row
        elif axis == 'col':
            return self.col
        else:
            return self.row, self.col

    def getIndexes(self, axis='all'):
        """Return the matrix indexes ([rows],[cols]) where axis='row'/'col'/'all'"""
        if axis == 'row':
            return self.ridx
        elif axis == 'col':
            return self.cidx
        else:
            return self.ridx, self.cidx

    def getByIdx(self, row, col):
        assert row in self.ridx, '\n%s is an invalid index' % row
        assert col in self.cidx, '\n%s is an invalid index' % col
        return self.array[row, col]

    def getByName(self, row, col):
        assert row in self.row, '\n%s is an invalid name' % row
        assert col in self.col, '\n%s is an invalid name' % col
        return self.array[self.row.index(row), self.col.index(col)]

    def setByIdx(self, row, col, val):
        assert row in self.ridx, '\n%s is an invalid index' % row
        assert col in self.cidx, '\n%s is an invalid index' % col
        self.array[row, col] = val

    def setByName(self, row, col, val):
        assert row in self.row, '\n%s is an invalid name' % row
        assert col in self.col, '\n%s is an invalid name' % col
        self.array[self.row.index(row), self.col.index(col)] = val


class MathArrayFunc(object):
    """A class of basic array functions some LAPACK based"""

    __doc__ = '''PySCeS array functions - used by Stoich'''

    array_kind = {'i': 0, 'l': 0, 'f': 0, 'd': 0, 'F': 1, 'D': 1}
    array_precision = {'i': 1, 'l': 1, 'f': 0, 'd': 1, 'F': 0, 'D': 1}
    array_type = [['f', 'd'], ['F', 'D']]
    LinAlgError = 'LinearAlgebraError'

    def commonType(self, *arrays):
        """
        commonType(\*arrays)

        Numeric detect and set array precision (will be replaced with new scipy.core compatible code when ready)

        Arguments:

        \*arrays: input arrays

        """
        kind = 0
        #    precision = 0
        #   force higher precision in lite version
        precision = 1
        for a in arrays:
            t = a.dtype.char
            kind = max(kind, self.array_kind[t])
            precision = max(precision, self.array_precision[t])

        return self.array_type[kind][precision]

    def castCopyAndTranspose(self, type, *arrays):
        """
        castCopyAndTranspose(type, \*arrays)

        Cast numeric arrays to required type and transpose

        Arguments:

        type: the required type to cast to
        \*arrays: the arrays to be processed

        """
        cast_arrays = ()
        for a in arrays:
            if a.typecode() == type:
                cast_arrays = cast_arrays + (copy.copy(scipy.transpose(a)),)
            else:
                cast_arrays = cast_arrays + (
                    copy.copy(scipy.transpose(a).astype(type)),
                )
        if len(cast_arrays) == 1:
            return cast_arrays[0]
        else:
            return cast_arrays

    def assertRank2(self, *arrays):
        """
        assertRank2(\*arrays)

        Check that we are using a 2D array

        Arguments:

        \*arrays: input array(s)

        """
        for a in arrays:
            if len(a.shape) != 2:
                raise LinAlgError('Array must be two-dimensional')

    def SwapCol(self, res_a, r1, r2):
        """
        SwapCol(res_a,r1,r2)

        Swap two columns using BLAS swap, arrays can be (or are upcast to) type double (d) or double complex (D).
        Returns the colswapped array

        Arguments:

        res_a: the input array
        r1: the first column to be swapped
        r2: the second column to be swapped

        """
        if (
            self.array_kind[self.commonType(res_a)] == 1
        ):  # brett 20041226 added complex support to swap
            res_a = res_a.astype('D')
            return self.SwapColz(res_a, r1, r2)
        else:
            res_a = res_a.astype('d')
            return self.SwapCold(res_a, r1, r2)

    def SwapRow(self, res_a, r1, r2):
        """
        SwapRow(res_a,r1,r2)

        Swaps two rows using BLAS swap, arrays can be (or are upcast to) type double (d) or double complex (D).
        Returns the rowswapped array.

        Arguments:

        res_a: the input array
        r1: the first row index to be swapped
        r2:  the second row index to be swapped

        """
        if (
            self.array_kind[self.commonType(res_a)] == 1
        ):  # brett 20041226 added complex support to swap
            res_a = res_a.astype('D')
            return self.SwapRowz(res_a, r1, r2)
        else:
            res_a = res_a.astype('d')
            return self.SwapRowd(res_a, r1, r2)

    def SwapElem(self, res_a, r1, r2):
        """
        SwapElem(res_a,r1,r2)

        Swaps two elements in a 1D vector

        Arguments:

        res_a: the input vector
        r1: index 1
        r2: index 2

        """
        res_a[r1], res_a[r2] = res_a[r2], res_a[r1]
        return res_a

    def SwapCold(self, res_a, c1, c2):
        """
        SwapCold(res_a,c1,c2)

        Swaps two double (d) columns in an array using BLAS DSWAP. Returns the colswapped array.

        Arguments:

        res_a: input array
        c1: column index 1
        c2: column index 2

        """
        res_a = res_a.astype('d')
        res_a[:, c1], res_a[:, c2] = myfblas.dswap(res_a[:, c1], res_a[:, c2])
        return res_a

    def SwapRowd(self, res_a, r1, r2):
        """
        SwapRowd(res_a,c1,c2)

        Swaps two double (d) rows in an array using BLAS DSWAP. Returns the rowswapped array.

        Arguments:

        res_a: input array
        c1: row index 1
        c2: row index 2

        """
        res_a = res_a.astype('d')
        res_a[r1, :], res_a[r2, :] = myfblas.dswap(res_a[r1, :], res_a[r2, :])
        return res_a

    def SwapColz(self, res_a, c1, c2):
        """
        SwapColz(res_a,c1,c2)

        Swaps two double complex (D) columns in an array using BLAS ZSWAP. Returns the colswapped array.

        Arguments:

        res_a: input array
        c1: column index 1
        c2: column index 2

        """
        res_a = res_a.astype('D')
        res_a[:, c1], res_a[:, c2] = myfblas.zswap(res_a[:, c1], res_a[:, c2])
        return res_a

    def SwapRowz(self, res_a, r1, r2):
        """
        SwapRowz(res_a,c1,c2)

        Swaps two double complex (D) rows in an array using BLAS ZSWAP. Returns the rowswapped array.

        Arguments:

        res_a: input array
        c1: row index 1
        c2: row index 2

        """
        res_a = res_a.astype('D')
        res_a[r1, :], res_a[r2, :] = myfblas.zswap(res_a[r1, :], res_a[r2, :])
        return res_a

    def MatrixFloatFix(self, mat, val=1.0e-15):
        """
        MatrixFloatFix(mat,val=1.e-15)

        Clean an array removing any floating point artifacts defined as being smaller than a specified value.
        Processes an array inplace

        Arguments:

        mat: the input 2D array
        val [default=1.e-15]: the threshold value (effective zero)

        """
        zero_vals = abs(mat) < val

        for x in range(mat.shape[0]):
            for y in range(mat.shape[1]):
                if zero_vals[x, y]:
                    mat[x, y] = 0.0
                '''
                if abs(mat[x,y]) != 0.0 and abs(mat[x,y]) < val:
                    mat[x,y] = round(mat[x,y])
                    if abs(mat[x,y]) < val:
                        #mat[x,y] = round(mat[x,y])
                        mat[x,y] = 0.0
                '''
        del zero_vals

    def MatrixValueCompare(self, matrix):
        """
        MatrixValueCompare(matrix)

        Finds the largest/smallest abs(value) > 0.0 in a matrix.
        Returns a tuple containing (smallest,largest) values

        Arguments:

        matrix: the input 2D array

        """
        val_B = 0.0
        val_S = 1.0e30
        for x in range(matrix.shape[0]):
            for y in range(matrix.shape[1]):
                if abs(matrix[x, y]) > abs(val_B) and abs(matrix[x, y]) != 0.0:
                    val_B = matrix[x, y]
                if abs(matrix[x, y]) < abs(val_S) and abs(matrix[x, y]) != 0.0:
                    val_S = matrix[x, y]
        return (val_S, val_B)


class Stoich(MathArrayFunc):
    '''PySCeS stoichiometric analysis class: initialized with a stoichiometric matrix N (input)'''

    __stoichdiagmode__ = 0
    __version__ = __version__
    __TimeFormat = "%H:%M:%S"
    USE_QR = False
    info_moiety_conserve = False

    def __init__(self, input):
        """Initialize class variables"""

        self.nmatrix = input
        row, col = self.nmatrix.shape
        self.nmatrix_col = tuple(scipy.array(list(range(col))))
        self.nmatrix_row = tuple(scipy.array(list(range(row))))

        # Create a machine specific instance
        from scipy import MachAr

        mach_spec = MachAr()

        self.stoichiometric_analysis_fp_zero = mach_spec.eps * 2.0e4
        self.stoichiometric_analysis_lu_precision = self.stoichiometric_analysis_fp_zero
        self.stoichiometric_analysis_gj_precision = (
            self.stoichiometric_analysis_lu_precision * 10.0
        )

        self.species = None
        self.reactions = None

    def AnalyseK(self):
        """
        AnalyseK()

        Evaluate the stoichiometric matrix and calculate the nullspace using LU decomposition and backsubstitution .
        Generates the MCA K and Ko arrays and associated row and column vectors

        Arguments:
        None

        """
        print('Calculating K matrix .', end=' ')

        self.info_flux_conserve = 0  # added 20020416 for conservation detection

        if self.__stoichdiagmode__:
            print('\nKMATRIX\nGetUpperMatrix: ' + time.strftime(self.__TimeFormat))

        p, u, row_vector, column_vector = self.GetUpperMatrix(self.nmatrix)
        # p,u,row_vector,column_vector = self.GetUpperMatrixUsingQR(self.nmatrix)
        if self.__stoichdiagmode__:
            print('\nKMATRIX\nScalePivots: ' + time.strftime(self.__TimeFormat))

        unipiv_a = self.ScalePivots(u)
        if self.__stoichdiagmode__:
            print('\nKMATRIX\nBackSubstitution: ' + time.strftime(self.__TimeFormat))

        R_a, row_vector, column_vector = self.BackSubstitution(
            unipiv_a, row_vector, column_vector
        )
        if self.__stoichdiagmode__:
            print('\nKMATRIX\nK_split_R: ' + time.strftime(self.__TimeFormat))

        (
            r_ipart,
            r_fpart,
            row_vector,
            column_vector,
            nullspace,
            r_fcolumns,
            self.info_flux_conserve,
        ) = self.K_split_R(R_a, row_vector, column_vector)

        if self.__stoichdiagmode__:
            print('\nKMATRIX\n' + time.strftime(self.__TimeFormat))

        self.kmatrix = nullspace
        self.kmatrix_row = tuple(row_vector)
        self.kmatrix_col = tuple(column_vector)

        self.kzeromatrix = r_fpart
        self.kzeromatrix_row = tuple(r_fcolumns)
        self.kzeromatrix_col = tuple(column_vector)
        print(' done.')

    def AnalyseL(self):
        """
        AnalyseL()

        Evaluate the stoichiometric matrix and calculate the left nullspace using LU factorization and backsubstitution.
        Generates the MCA L, Lo, Nr and Conservation matrix and associated row and column vectors

        Arguments:
        None

        """
        print('Calculating L matrix .', end=' ')

        a = scipy.transpose(self.nmatrix)
        self.info_moiety_conserve = False  # added 20020416 for conservation detection

        if self.__stoichdiagmode__:
            print('\nLMATRIX\nGetUpperMatrix: ' + time.strftime(self.__TimeFormat))

        if not self.USE_QR:
            p, u, row_vector, column_vector = self.GetUpperMatrix(a)
        else:
            p, u, row_vector, column_vector = self.GetUpperMatrixUsingQR(a)
        if self.__stoichdiagmode__:
            print('\nLMATRIX\nScalePivots: ' + time.strftime(self.__TimeFormat))

        unipiv_a = self.ScalePivots(u)
        if self.__stoichdiagmode__:
            print('\nLMATRIX\nBackSubstitution: ' + time.strftime(self.__TimeFormat))

        R_a, row_vector, column_vector = self.BackSubstitution(
            unipiv_a, row_vector, column_vector
        )
        if self.__stoichdiagmode__:
            print('\nLMATRIX\nK_split_R: ' + time.strftime(self.__TimeFormat))

        (
            r_ipart,
            consmatrix,
            cons_row_vector,
            cons_col_vector,
            lmatrix,
            lmatrix_row_vector,
            lmatrix_col_vector,
            lomatrix,
            lomatrix_row_vector,
            lomatrix_col_vector,
            nrmatrix,
            Nred_vector,
            Nred_vector_col,
            self.info_moiety_conserve,
        ) = self.L_split_R(a, R_a, row_vector, column_vector)

        if self.__stoichdiagmode__:
            print('\nLMATRIX\n' + time.strftime(self.__TimeFormat))

        self.lmatrix = lmatrix
        self.lmatrix_row = tuple(lmatrix_row_vector)
        self.lmatrix_col = tuple(lmatrix_col_vector)

        self.lzeromatrix = lomatrix
        self.lzeromatrix_row = tuple(lomatrix_row_vector)
        self.lzeromatrix_col = tuple(lomatrix_col_vector)

        self.conservation_matrix = consmatrix
        self.conservation_matrix_row = tuple(cons_row_vector)
        self.conservation_matrix_col = tuple(cons_col_vector)

        self.nrmatrix = nrmatrix
        self.nrmatrix_row = tuple(Nred_vector)
        self.nrmatrix_col = tuple(scipy.copy(self.nmatrix_col))
        print(' done.')

    def PivotSort(self, a, row_vector, column_vector):
        """
        PivotSort(a,row_vector,column_vector)

        This is a sorting routine that accepts a matrix and row/colum vectors
        and then sorts them so that: there are no zero rows (by swapping with first
        non-zero row) The abs(largest) pivots are moved onto the diagonal to maintain
        numerical stability. Row and column swaps are recorded in the tracking vectors.

        Arguments:

        a: the input array
        row_vector: row tracking vector
        column_vector: column tracking vector

        """
        t = self.commonType(a)
        row, col = a.shape

        for z in range(0, min(row, col)):
            if abs(a[z, z]) < self.stoichiometric_analysis_lu_precision:
                maxV = self.stoichiometric_analysis_lu_precision
                maxP = (None, None)
                zeroP = (None, None)

                mList = []
                for x in range(z, min(row, col)):
                    mList.append((x, max(abs(a[x, x:]))))
                mVal = 0.0
                mRow = None
                # print mList
                for el in mList:
                    if el[1] > mVal:
                        mVal = el[1]
                        mRow = el[0]
                if self.__stoichdiagmode__:
                    print('\tmVal', mVal, mRow)

                if mRow != None:
                    for y in range(z, col):
                        if (
                            abs(a[mRow, y]) > abs(maxV)
                            and abs(a[mRow, y])
                            > self.stoichiometric_analysis_lu_precision
                        ):
                            maxV = a[mRow, y]
                            maxP = (mRow, y)

                if zeroP != (None, None) and zeroP != (z, z) and mRow != None:
                    if self.__stoichdiagmode__:
                        print('  Swapping: ', (z, z), (zeroP[0], zeroP[1]), 1.0)
                    a = self.SwapRowd(a, z, zeroP[0])
                    a = self.SwapCold(a, z, zeroP[1])
                    row_vector = self.SwapElem(row_vector, z, zeroP[0])
                    column_vector = self.SwapElem(column_vector, z, zeroP[1])
                elif (
                    maxP != (None, None)
                    and maxP != (min(row, col), min(row, col))
                    and mRow != None
                ):
                    if self.__stoichdiagmode__:
                        print('  Swapping: ', (z, z), (maxP[0], maxP[1]), a[z, z], maxV)
                    a = self.SwapRowd(a, z, maxP[0])
                    a = self.SwapCold(a, z, maxP[1])
                    row_vector = self.SwapElem(row_vector, z, maxP[0])
                    column_vector = self.SwapElem(column_vector, z, maxP[1])

        return (a, row_vector, column_vector)

    def PivotSort_initial(self, a, row_vector, column_vector):
        """
        PivotSort_initial(a,row_vector,column_vector)

        This is a sorting routine that accepts a matrix and row/colum vectors
        and then sorts them so that: the abs(largest) pivots are moved onto the diagonal to maintain
        numerical stability i.e. the matrix diagonal is in descending max(abs(value)).
        Row and column swaps are recorded in the tracking vectors.

        Arguments:

        a: the input array
        row_vector: row tracking vector
        column_vector: column tracking vector

        """

        # SAME AS THE ABOVE JUST DOES ALL VALUES NOT ONLY NON_ZERO ONES
        # brett 200500802

        t = self.commonType(a)
        row, col = a.shape

        for z in range(0, min(row, col)):
            if z == 0 or abs(a[z, z]) < abs(a[z - 1, z - 1]):
                maxV = self.stoichiometric_analysis_lu_precision
                maxP = (None, None)
                zeroP = (None, None)
                mList = []
                for x in range(z, min(row, col)):
                    mList.append((x, max(abs(a[x, x:]))))
                mVal = 0.0
                mRow = None
                for el in mList:
                    if el[1] > mVal:
                        mVal = el[1]
                        mRow = el[0]
                if self.__stoichdiagmode__:
                    print('\tmVal', mVal, mRow)

                if mRow != None:
                    for y in range(z, col):
                        if (
                            abs(a[x, y]) > abs(maxV)
                            and abs(a[x, y]) > self.stoichiometric_analysis_lu_precision
                        ):
                            maxV = a[x, y]
                            maxP = (x, y)

                if (
                    maxP != (None, None)
                    and maxP != (min(row, col), min(row, col))
                    and maxP != (z, z)
                ):
                    if self.__stoichdiagmode__:
                        print('  Swapping: ', (z, z), (maxP[0], maxP[1]), a[z, z], maxV)
                    a = self.SwapRowd(a, z, maxP[0])
                    a = self.SwapCold(a, z, maxP[1])
                    row_vector = self.SwapElem(row_vector, z, maxP[0])
                    column_vector = self.SwapElem(column_vector, z, maxP[1])

        return (a, row_vector, column_vector)

    def PLUfactorize(self, a_in):
        """
        PLUfactorize(a_in)

        Performs an LU factorization using LAPACK D/ZGetrf. Now optimized for FLAPACK interface.
        Returns LU - combined factorization, IP - rowswap information and info - Getrf error control.

        Arguments:

        a_in: the matrix to be factorized

        """
        print('.', end=' ')

        self.assertRank2(a_in)
        t = self.commonType(a_in)
        if a_in.dtype.char == 'D':
            # print 'Complex matrix ' + a_in.typecode()
            ##  a = copy.copy(scipy.transpose(a_in))
            # brett 201106 flapack optimize
            a = a_in.copy()
        elif a_in.dtype.char == 'd':
            # print 'Float matrix ' + a_in.typecode()
            ##  a = copy.copy(scipy.transpose(a_in))
            # brett 201106 flapack optimize
            a = a_in.copy()
        else:
            # print 'Other matrix casting to double, was: ' + a_in.typecode()
            ##  a = copy.copy(scipy.transpose(a_in).astype('d'))
            # brett 201106 flapack optimize
            a = a_in.copy().astype('d')
        Using_FLAPACK = 1

        if Using_FLAPACK == 1:
            try:
                if self.array_kind[t] == 1:
                    getrf = myflapack.zgetrf  # scipy flapack 20041226
                else:
                    getrf = myflapack.dgetrf  # scipy flapack 20041226
            except Exception as e:
                print("FLAPACK error", e)
        ##  else:
        ##  try:
        ##  if self.array_kind[t] == 1:
        ##  getrf = scipy.linalg.clapack.zgetrf #scipy clapack (ATLAS) 20030605
        ##  else:
        ##  getrf = scipy.linalg.clapack.dgetrf #scipy clapack (ATLAS) 20030605
        ##  except Exception, e:
        ##  print "CLAPACK error", e

        if Using_FLAPACK == 1:
            ##  results = getrf(scipy.transpose(a)) # brett 20041226
            results = getrf(a)  # brett 201106

        results = list(results)

        if results[2] < 0:
            print('Argument ', results['info'], ' had an illegal value')
            raise LinAlgError
        elif results[2] > 0:
            pass

        if Using_FLAPACK == 1:
            result = results[0]  # brett 20041226
        ##  else:
        ##  result = scipy.transpose(results[0]) # -- this is normal

        badlist = []
        for x in range(result.shape[0]):
            for y in range(result.shape[1]):
                if (
                    abs(result[x, y]) != 0.0
                    and abs(result[x, y]) < self.stoichiometric_analysis_lu_precision
                ):
                    result[x, y] = 0.0  # gets rid of the -0.0 irritation
                    if len(badlist) == 0 and (x, y) == (
                        x,
                        x,
                    ):  # catch 1st float on a pivot
                        badlist.append(x)
        if len(badlist) != 0:
            results[2] = (
                badlist[0] - 1
            )  # 20040423 if float was on a pivot - refactorize

        return (result, results[1], results[2])  # scipy cblas (ATLAS?) 20030506

    """
    def PLUfactorizeOLD(self,a_in):
        '''
        PLUfactorize(a_in) OLD PRE SCIPY 0.10 uses clapack

        Performs an LU factorization using LAPACK D/ZGetrf.
        Returns LU - combined factorization, IP - rowswap information and info - Getrf error control.

        Arguments:

        a_in: the matrix to be factorized

        '''
        print '.',

        self.assertRank2(a_in)
        t = self.commonType(a_in)
        if a_in.dtype.char == 'D':
            #print 'Complex matrix ' + a_in.typecode()
            a = copy.copy(scipy.transpose(a_in))
        elif a_in.dtype.char == 'd':
            #print 'Float matrix ' + a_in.typecode()
            a = copy.copy(scipy.transpose(a_in))
        else:
            #print 'Other matrix casting to double, was: ' + a_in.typecode()
            a = copy.copy(scipy.transpose(a_in).astype('d'))
        Using_FLAPACK = 0
        # brett 20041226 - protecting ourselves against flapack
        try:
            scipy.linalg.clapack.empty_module()
            Using_FLAPACK = 1
            print "Using FLAPACK"
        except:
            print "Using CLAPACK"

        if Using_FLAPACK == 1:
            try:
                if self.array_kind[t] == 1:
                    getrf = myflapack.zgetrf #scipy flapack 20041226
                else:
                    getrf = myflapack.dgetrf #scipy flapack 20041226
            except Exception, e:
                print "FLAPACK error", e
        else:
            try:
                if self.array_kind[t] == 1:
                    getrf = scipy.linalg.clapack.zgetrf #scipy clapack (ATLAS) 20030605
                else:
                    getrf = scipy.linalg.clapack.dgetrf #scipy clapack (ATLAS) 20030605
            except Exception, e:
                print "CLAPACK error", e

        if Using_FLAPACK == 1:
            results = getrf(scipy.transpose(a)) # brett 20041226
        else:
            # This is a $%^&*( ... suddenly with latest cvs f2py getrf only accepts arrays
            # not vectors so this song and dance is necessary to fix this 'behaviour'
            # as idiotic as it seems, we pad the vector with a zero row (no difference numerically)
            # run it through getrf and then strip it afterwards !@#$%^&*() - brett 20050707
            if a.shape[0] == 1:
                tarr = scipy.zeros((a.shape[0]+1,a.shape[1]),'d')
                tarr[0] = a[0]

                results = getrf(tarr) #scipy cblas (ATLAS) 20030506 -- this is normal

                results = list(results)
                results[0] = scipy.array([results[0][0]])
                results[1] = scipy.array([results[1][0]],'i')
                results[2] = results[2]
                results = tuple(results)
            else:
                results = getrf(a) #scipy cblas (ATLAS) 20030506 -- this is normal

        results = list(results)

        if results[2] < 0:
            print('Argument ', results['info'], ' had an illegal value')
            raise LinAlgError
        elif results[2] > 0:
            pass

        if Using_FLAPACK == 1:
            result = results[0] # brett 20041226
        else:
            result = scipy.transpose(results[0]) # -- this is normal


        badlist = []
        for x in range(result.shape[0]):
            for y in range(result.shape[1]):
                if abs(result[x,y]) != 0.0 and abs(result[x,y]) < self.stoichiometric_analysis_lu_precision:
                    result[x,y] = 0.0        # gets rid of the -0.0 irritation
                    if len(badlist) == 0 and (x,y) == (x,x): # catch 1st float on a pivot
                        badlist.append(x)
        if len(badlist) != 0:
            results[2] = badlist[0]-1 # 20040423 if float was on a pivot - refactorize

        return(result,results[1],results[2]) #scipy cblas (ATLAS?) 20030506
        """

    def SplitLU(self, plu, row, col, t=None):
        """
        SplitLU(plu,row,col,t)

        PLU takes the combined LU factorization computed by PLUfactorize and extracts the upper matrix.
        Returns U.

        Arguments:

        plu: LU factorization
        row: row tracking vector
        col: column tracking vector
        t [default=None)]: typecode argument (currently not used)

        """
        print('.', end=' ')

        if self.__stoichdiagmode__:
            print('\n', row, col)
        for cl in range(0, min(row, col)):
            plu[cl + 1 :, cl] = 0.0
        return plu

    def GetUpperMatrix(self, a):
        """
        GetUpperMatrix(a)

        Core analysis algorithm; an input is preconditioned using PivotSort_initial and then cycles of PLUfactorize and
        PivotSort are run until the factorization is completed. During this process the matrix is reordered by
        column swaps which emulates a full pivoting LU factorization. Returns the pivot matrix P, upper factorization U
        as well as the row/col tracking vectors.

        Arguments:

        a: a stoichiometric matrix

        """
        print('.', end=' ')

        t = self.commonType(a)
        row, col = a.shape

        # this is a test brett 20050802
        row_vector = scipy.array((list(range(row))))
        column_vector = scipy.array((list(range(col))))
        a, row_vector, column_vector = self.PivotSort_initial(
            a, row_vector, column_vector
        )

        # 18/09/2000 PLUfactorize included in self.GetUpperMatrix
        a, ip, info = self.PLUfactorize(a)

        # get U from getrf's PLU
        if self.__stoichdiagmode__:
            print(info)
        upper_out = self.SplitLU(a, row, col, t)

        p_out = scipy.identity(row)
        for x in range(0, min(row, col)):
            if x + 1 != ip[x]:
                p_out = self.SwapRowd(p_out, x, (ip[x] - 1))
                row_vector = self.SwapElem(row_vector, x, (ip[x] - 1))

        '''31/08/2000 Added to try to make sure that pivots exist only on the upper left side of the matrix
        max(abs(upper_out[x,:])) used instead of abs(max(upper_out[x,:])) otherwise for rows where
        all elements < 0 zero is maximum, this way the max of positive values is used'''

        upper_out, row_vector, column_vector = self.PivotSort(
            upper_out, row_vector, column_vector
        )

        '''20/09/2000 This bit sorts out the echelon matrix by running (if needed) cycles of
        self.PLUfactorize, self.GetUpperMatrix, and pivsort until only a staircase matrix remains. It uses both the error
        generated by self.PLUfactorize(info) and go_flag to control itself'''

        if info > 0:
            # Here we check to see if the error is not in the last or reduced-last column if it is - exit
            #        go_flag = 'yes' # 2001/04/26 changed for Python21 and future compatibility
            go_flag = 1
            for x in range(0, min(row, col)):
                if abs(upper_out[x, x]) > self.stoichiometric_analysis_lu_precision:
                    pos_holder = x
            if pos_holder + 1 < info and pos_holder + 1 < row:
                if (
                    max(abs(upper_out[pos_holder + 1, :]))
                    < self.stoichiometric_analysis_lu_precision
                ):
                    #                go_flag = 'no' # 2001/04/26 changed for Python21 and future compatibility
                    go_flag = 0
            elif pos_holder + 1 == info and pos_holder + 1 == min(row, col):
                #            go_flag = 'no' # 2001/04/26 changed for Python21 and future compatibility
                go_flag = 0
            #        while info > 0 and go_flag == 'yes': # 2001/04/26 changed for Python21 and future compatibility
            while info > 0 and go_flag == 1:
                # sort
                # upper_out,row_vector,column_vector = pivot_sort2k5(upper_out,row_vector,column_vector)
                upper_out, row_vector, column_vector = self.PivotSort(
                    upper_out, row_vector, column_vector
                )
                # PLUfactorize
                if self.__stoichdiagmode__:
                    print(
                        ' Echelon: ' + time.strftime(self.__TimeFormat), '(', info, ')'
                    )
                upper_out, ip, info = self.PLUfactorize(upper_out)

                # get U from getrf's PLU
                if self.__stoichdiagmode__:
                    print(info)
                upper_out = self.SplitLU(upper_out, row, col, t)

                # realign row vector and permutation matrix if necessary
                for x in range(0, min(row, col)):
                    if x + 1 != ip[x]:
                        p_out = self.SwapRowd(p_out, x, (ip[x] - 1))
                        row_vector = self.SwapElem(row_vector, x, (ip[x] - 1))
                # test if we can exit -- same criteria as before
                for x in range(0, min(row, col)):
                    if abs(upper_out[x, x]) > self.stoichiometric_analysis_lu_precision:
                        pos_holder = x
                if pos_holder + 1 < info and pos_holder + 1 < row:
                    if (
                        max(abs(upper_out[pos_holder + 1, :]))
                        < self.stoichiometric_analysis_lu_precision
                    ):
                        #                    go_flag = 'no' # 2001/04/26 changed for Python21 and future compatibility
                        go_flag = 0
                elif pos_holder + 1 == info and pos_holder + 1 == min(row, col):
                    #                go_flag = 'no' # 2001/04/26 changed for Python21 and future compatibility
                    go_flag = 0

        '''This bit will get rid of any zero rows so that we will hopefully only have to work with a
        reduced matrix after this, for completeness the row_vector will also be sliced'''

        for x in range(0, min(row, col)):
            if abs(upper_out[x, x]) > self.stoichiometric_analysis_lu_precision:
                pos_holder = x

        upper_out_r = scipy.zeros((pos_holder + 1, col)).astype(t)
        upper_out_r = upper_out[: pos_holder + 1, :]
        row_vector_r = row_vector[: pos_holder + 1]

        return (p_out, upper_out_r, row_vector_r, column_vector)

    def GetUpperMatrixUsingQR(self, a):
        """
        GetUpperMatrix(a)

        Core analysis algorithm; an input is preconditioned using PivotSort_initial and then cycles of PLUfactorize and
        PivotSort are run until the factorization is completed. During this process the matrix is reordered by
        column swaps which emulates a full pivoting LU factorization. Returns the pivot matrix P, upper factorization U
        as well as the row/col tracking vectors.

        Arguments:

        a: a stoichiometric matrix

        """
        print('.', end=' ')

        t = self.commonType(a)
        row, col = a.shape

        # this is a test brett 20050802
        row_vector = scipy.array((list(range(row))))
        column_vector = scipy.array((list(range(col))))
        ##  a,row_vector,column_vector = self.PivotSort(a,row_vector,column_vector)
        a, row_vector, column_vector = self.PivotSort_initial(
            a, row_vector, column_vector
        )

        Q, upper_out = scipy.linalg.qr(a)
        self.MatrixFloatFix(
            upper_out, val=self.stoichiometric_analysis_lu_precision * 10.0
        )

        '''This bit will get rid of any zero rows so that we will hopefully only have to work with a
        reduced matrix after this, for completeness the row_vector will also be sliced'''

        for x in range(0, min(row, col)):
            if abs(upper_out[x, x]) > self.stoichiometric_analysis_lu_precision:
                pos_holder = x

        upper_out_r = scipy.zeros((pos_holder + 1, col)).astype(t)
        p_out = scipy.diag(upper_out.shape[1] * [1.0])
        upper_out_r = upper_out[: pos_holder + 1, :]
        row_vector_r = row_vector[: pos_holder + 1]

        return (p_out, upper_out_r, row_vector_r, column_vector)

    def ScalePivots(self, a_one):
        """
        ScalePivots(a_one)

        Given an upper triangular matrix U, this method scales the diagonal (pivot values) to one.

        Arguments:

        a_one: an upper triangular matrix U

        """
        print('.', end=' ')

        t = self.commonType(a_one)
        row, col = a_one.shape

        '''13/09/2000 We now assume that the matrix has the correct shape ie. the pivots are in a
        perfect staircase'''

        for x in range(0, row):
            if (
                abs(a_one[x, x]) == 0.0
            ):  # ths is to prevent NAN's when FixFloat zeros a supersmall pivot
                pass
            elif abs(a_one[x, x]) != 1.0:
                a_one[x, :] = a_one[x, :] / abs(a_one[x, x])
            if a_one[x, x] < 0.0:
                a_one[x, :] = -(a_one[x, :])
        return a_one

    def BackSubstitution(self, res_a, row_vector, column_vector):
        """
        BackSubstitution(res_a,row_vector,column_vector)

        Jordan reduction of a scaled upper triangular matrix. The returned array is now in the form [I R] and can
        be used for nullspace determination. Modified row and column tracking vetors are also returned.

        Arguments:

        res_a: unitary pivot upper triangular matrix
        row_vector: row tracking vector
        column_vector: column tracking vector

        """
        print('.', end=' ')

        t = self.commonType(res_a)
        row, col = res_a.shape

        '''
        13/09/2000 removed the copy thing, and eliminated the divnumb variable. Because we
        are now working with a row reduced matrix upward elimination only works on the pivot cols
        '''

        # old back substitution circa 2000! brett - 20050801
        '''
        bigF = 0
        if max(res_a.shape) > 500:
            bigF = 1
        for y in range (min(row,col)-1,-1,-1):
            if bigF and self.__stoichdiagmode__:
                print '.',
            for x in range (min(row,col)-2,-1,-1):
                z = x + 1
                while z < min(row,col):
                    if abs(res_a[x,y]) > self.stoichiometric_analysis_gj_precision and abs(res_a[z,y]) > self.stoichiometric_analysis_gj_precision:
                        res_a[x,:] = res_a[x,:] - res_a[z,:]*(res_a[x,y]/res_a[z,y])
                        z = z + 1
                        continue
                    else:
                        z = z + 1
        '''
        # new back substitution
        # right looking algorithm that uses array slicing speeds up as it moves down the pivots
        # this algorithm is at least 3X as fast as the old one for a large system - brett 20050805
        bigF = 0
        if max(res_a.shape) > 500:
            bigF = 1
        for x in range(min(row, col)):
            if bigF and self.__stoichdiagmode__:
                print('.', end=' ')
            for y in range(x + 1, min(row, col)):
                if abs(res_a[x, y]) > self.stoichiometric_analysis_lu_precision:
                    res_a[x, y:] = res_a[x, y:] - (res_a[x, y] * res_a[y, y:])

        # $%^&* wtf is this for??? brett - 20050801 (besides cleaning fp wierdness ... nothing ;-)
        self.MatrixFloatFix(res_a, val=self.stoichiometric_analysis_gj_precision)
        res_a, row_vector, column_vector = self.PivotSort(
            res_a, row_vector, column_vector
        )

        return (res_a, row_vector, column_vector)

    def K_split_R(self, R_a, row_vector, column_vector):
        """
        K_split_R(R_a,row_vector,column_vector)

        Using the R factorized form of the stoichiometric matrix we now form the K and Ko matrices. Returns
        the r_ipart,Komatrix,Krow,Kcolumn,Kmatrix,Korow,info

        Arguments:

        R_a: the Gauss-Jordan reduced stoichiometric matrix
        row_vector: row tracking vector
        column_vector: column tracking vector

        """
        print('.', end=' ')

        t = self.commonType(R_a)
        row, col = R_a.shape

        '''14/09/2000 Seeing as we now should have a perfect staircase in a reduced matrix we do
        not have to search for the last pivot it will exist at min(row,col)-1 so the pos_holder
        finding code has been removed (actually moved into self.GetUpperMatrix)'''

        pos_holder = min(row, col) - 1

        '''This bit extracts the identity part from R (future note this could be replaced by an I matrix formed by min(row,col))'''

        r_ipart = scipy.zeros((pos_holder + 1, pos_holder + 1)).astype(t)
        r_ipart = R_a[: pos_holder + 1, : pos_holder + 1]

        row_i, col_i = r_ipart.shape

        '''If there are free variables, then this bit will extract them and form the row/col vectors'''

        empty_rf = 0
        K_switch = 0  # added 20020416 class attribute for conservation detection 0 = none, 1 = exists

        if col - col_i > self.stoichiometric_analysis_fp_zero:
            r_fpart = scipy.zeros((pos_holder + 1, col - (pos_holder + 1))).astype(t)
            r_fpart = R_a[: pos_holder + 1, pos_holder + 1 :]

            row_vector_dependent = column_vector[: pos_holder + 1]
            row_vector_independent = column_vector[pos_holder + 1 :]

            # row_vector = scipy.concatenate((row_vector_independent,row_vector_dependent),1)
            row_vector = scipy.hstack(
                (row_vector_independent, row_vector_dependent)
            )  # numpy 0.10
            column_vector = column_vector[pos_holder + 1 :]
            K_switch = 1
        else:
            print('no flux conservation')
            r_fpart = 'no flux conservation'
            empty_rf = 1
            K_switch = 0

        #    if r_fpart == 'F is empty, R = I': # 2001/04/26 changed for Python21 compatibility
        if empty_rf == 1:
            return (
                r_ipart,
                r_ipart,
                column_vector,
                column_vector,
                r_ipart,
                column_vector,
                K_switch,
            )
        else:
            row, col = r_fpart.shape
            id = scipy.identity(col).astype(t)
            nullspace = scipy.concatenate((id, -r_fpart), 0).astype(t)

        # brett 05/11/2002 changed r_fpart to -r_fpart
        return (
            r_ipart,
            -r_fpart,
            row_vector,
            column_vector,
            nullspace,
            row_vector_dependent,
            K_switch,
        )

    def L_split_R(self, Nfull, R_a, row_vector, column_vector):
        """
        L_split_R(Nfull,R_a,row_vector,column_vector)

        Takes the Gauss-Jordan factorized N^T and extract the L, Lo, conservation (I -Lo) and reduced stoichiometric matrices. Returns: lmatrix_col_vector, lomatrix, lomatrix_row, lomatrix_co, nrmatrix, Nred_vector_row, Nred_vector_col, info

        Arguments:

        Nfull: the original stoichiometric matrix N
        R_a: gauss-jordan factorized form of N^T
        row_vector: row tracking vector
        column_vector: column tracking vector

        """
        print('.', end=' ')

        # print '\nR_a'
        # print `R_a`
        t = self.commonType(R_a)
        row, col = R_a.shape

        '''14/09/2000 Seeing as we now should have a perfect staircase in a reduced matrix we do
        not have to search for the last pivot it will exist at min(row,col)-1 so the pos_holder
        finding code has been removed (actually moved into self.GetUpperMatrix)'''

        pos_holder = min(row, col) - 1

        '''Here we extract the identity matrix from R (future note this could be replaced by an I matrix formed by min(row,col))'''

        r_ipart = scipy.zeros((pos_holder + 1, pos_holder + 1)).astype(t)
        r_ipart = R_a[: pos_holder + 1, : pos_holder + 1]

        row_i, col_i = r_ipart.shape

        #    exit1 = 'no' # 2001/04/26 changed for Python21 and future compatibility
        exit1 = 0
        L_switch = False  # added 20020416 class attribute for conservation detection 0 = none, 1 = exists

        '''If there are free variable then they are extracted and packaged, row/col vectors are formed'''

        if col - col_i > self.stoichiometric_analysis_fp_zero:
            r_fpart = scipy.zeros((pos_holder + 1, col - (pos_holder + 1))).astype(t)
            r_fpart = R_a[: pos_holder + 1, pos_holder + 1 :]

            col_vector_dependent = column_vector[pos_holder + 1 :]
            col_vector_independent = column_vector[: pos_holder + 1]
            # cons_col_vector = scipy.concatenate((col_vector_independent,col_vector_dependent),1)
            cons_col_vector = scipy.hstack(
                (col_vector_independent, col_vector_dependent)
            )  # numpy 0.10
            cons_row_vector = col_vector_dependent
            # lmatrix_row_vector = scipy.concatenate((col_vector_independent,col_vector_dependent),1)
            lmatrix_row_vector = scipy.hstack(
                (col_vector_independent, col_vector_dependent)
            )  # numpy 0.10
            lmatrix_col_vector = col_vector_independent
            lomatrix_row_vector = col_vector_dependent
            lomatrix_col_vector = col_vector_independent
            Nred_vector = col_vector_independent
            L_switch = True
        else:
            lomatrix_row_vector = column_vector[: pos_holder + 1]
            lomatrix_col_vector = column_vector[: pos_holder + 1]
            Nred_vector = column_vector[: pos_holder + 1]
            # print('F is empty, R = I, no conservation matrix, L is Lo')
            r_fpart = r_ipart
            #        exit1 = 'yes' # 2001/04/26 changed for Python21 and future compatibility
            exit1 = 1
            L_switch = False

        #    if exit1 == 'yes': # 2001/04/26 changed for Python21 and future compatibility
        if exit1 == 1:
            r_fpart = scipy.transpose(r_fpart)
            lmatrix = r_fpart
            # 02/10/2000 removed so that the thing returns the Lo matrix as L
            # r_fpart = 'no conservation Lo = L = I'
            # my factorization routines now swap things around for numeric stability, so that Nr might be a row/column swapped
            # this simply synchronizes Nr with its labels - brett 20050805
            Nfull = scipy.transpose(Nfull)
            row, col = Nfull.shape
            Nred = scipy.zeros((row_i, col)).astype(t)
            for x in range(0, row_i):
                Nred[x, :] = Nfull[Nred_vector[x], :]
            # return the right stuff - brett 20050805
            return (
                r_ipart,
                'no conservation',
                'no conservation',
                'no conservation',
                lmatrix,
                lomatrix_row_vector,
                lomatrix_col_vector,
                lmatrix,
                lomatrix_row_vector,
                lomatrix_col_vector,
                Nred,
                Nred_vector,
                row_vector,
                L_switch,
            )
            # return(r_ipart,'no conservation','no conservation','no conservation',lmatrix,lomatrix_row_vector,lomatrix_col_vector,lmatrix,lomatrix_row_vector,lomatrix_col_vector,scipy.transpose(Nfull),Nred_vector,row_vector,L_switch)
        else:
            r_fpart = scipy.transpose(r_fpart)
            row, col = r_fpart.shape

            id = scipy.identity(row).astype(t)
            # consmatrix = scipy.concatenate((-r_fpart,id),1).astype(t)
            consmatrix = scipy.hstack((-r_fpart, id)).astype(t)  # numpy 0.10

            id = scipy.identity(col).astype(t)
            lmatrix = scipy.concatenate((id, r_fpart), 0).astype(t)

            '''This bit creates Nr. The transpose is only necessary if Nfull is already transposed in the input function'''

            Nfull = scipy.transpose(Nfull)
            row, col = Nfull.shape
            Nred = scipy.zeros((row_i, col)).astype(t)
            for x in range(0, row_i):
                Nred[x, :] = Nfull[Nred_vector[x], :]

            return (
                r_ipart,
                consmatrix,
                cons_row_vector,
                cons_col_vector,
                lmatrix,
                lmatrix_row_vector,
                lmatrix_col_vector,
                r_fpart,
                lomatrix_row_vector,
                lomatrix_col_vector,
                Nred,
                Nred_vector,
                row_vector,
                L_switch,
            )

    def SVD_Rank_Check(self, matrix=None, factor=1.0e4, resultback=0):
        """
        SVD_Rank_Check(matrix=None,factor=1.0e4,resultback=0)

        Calculates the dimensions of L/L0/K/K) by way of SVD and compares them to the Guass-Jordan results. Please note that for LARGE ill conditioned matrices the SVD can become numerically unstable when used for nullspace determinations

        Arguments:

        matrix [default=None]: the stoichiometric matrix default is self.Nmatrix
        factor [default=1.0e4]: factor used to calculate the 'zero pivot' mask = mach_eps*factor
        resultback [default=0]: return the SVD results, U, S, vh

        """
        if matrix == None:
            matrix = self.nmatrix

        nrow, ncol = matrix.shape
        TrMat = 0
        if ncol > nrow:
            # 'INFO: SVD is more accurate for tall matrices using (N)T\n'
            matrix = scipy.transpose(matrix)
            TrMat = 1

        u, s, vh = scipy.linalg.svd(matrix)
        maskF = scipy.machar.machar_double.eps * factor

        if TrMat == 0:
            print('SVD zero mask:', maskF)
            print(
                'LU factorization effective zero:',
                self.stoichiometric_analysis_lu_precision,
            )

            rank = len([el for el in (abs(s) > maskF) if el > 0.0])
            ##            null_mask = (abs(s) > maskF)
            ##            vhnull = scipy.compress(null_mask,vh,axis=0)
            ##            unull = scipy.compress(null_mask,u,axis=0)
            ##            assert vhnull.shape[0] == unull.shape[0], 'This should be true or I\'m very confused'
            ##            rank = vhnull.shape[0]

            print(
                '\nNmatrix has ' + repr(nrow) + ' rows and ' + repr(ncol) + ' columns'
            )

            print('\nSVD \"considers\" the rank to be:        ' + repr(rank))
            print(
                'LU (Kmatrix) considers the rank to be: '
                + repr(self.kzeromatrix.shape[0])
            )
            print(
                'LU (Lmatrix) considers the rank to be: '
                + repr(self.lzeromatrix.shape[1])
            )

            print('\nComparing lzeromatrix dimensions')
            print('SVD lzeromatrix.shape =', (u.shape[0] - rank, rank))
            print('LU  lzeromatrix.shape =', self.lzeromatrix.shape)

            print('\nComparing lmatrix dimensions')
            print('SVD lmatrix.shape =', (u.shape[0], rank))
            print('LU  lmatrix.shape =', self.lmatrix.shape)

            print('\nComparing kzeromatrix dimensions')
            print('SVD kzeromatrix.shape =', (rank, vh.shape[0] - rank))
            print('LU  kzeromatrix.shape =', self.kzeromatrix.shape)

            print('\nComparing kmatrix dimensions')
            print('SVD kmatrix.shape =', (vh.shape[0], vh.shape[0] - rank))
            print('LU  kmatrix.shape =', self.kmatrix.shape)
        else:
            print('SVD zero mask:', maskF)
            print(
                'LU factorization effective zero:',
                self.stoichiometric_analysis_lu_precision,
            )

            rank = len([el for el in (abs(s) > maskF) if el > 0.0])
            ##            null_mask = (abs(s) > maskF)
            ##            vhnull = scipy.compress(null_mask,vh,axis=0)
            ##            unull = scipy.compress(null_mask,u,axis=0)
            ##            assert vhnull.shape[0] == unull.shape[0], 'This should be true or I\'m very confused'
            ##            rank = vhnull.shape[0]

            print(
                '\nNmatrix has ' + repr(nrow) + ' rows and ' + repr(ncol) + ' columns'
            )

            print('\nSVD considers the rank to be:          ' + repr(rank))
            print(
                'LU (Kmatrix) considers the rank to be: '
                + repr(self.kzeromatrix.shape[0])
            )
            print(
                'LU (Lmatrix) considers the rank to be: '
                + repr(self.lzeromatrix.shape[1])
            )

            print('\nComparing lzeromatrix dimensions')
            print('SVD lzeromatrix.shape =', (vh.shape[0] - rank, rank))
            print('LU  lzeromatrix.shape =', self.lzeromatrix.shape)

            print('\nComparing lmatrix dimensions')
            print('SVD lmatrix.shape =', (vh.shape[0], rank))
            print('LU  lmatrix.shape =', self.lmatrix.shape)

            print('\nComparing kzeromatrix dimensions')
            print('SVD kzeromatrix.shape =', (rank, u.shape[0] - rank))
            print('LU  kzeromatrix.shape =', self.kzeromatrix.shape)

            print('\nComparing kmatrix dimensions')
            print('SVD kmatrix.shape =', (u.shape[0], u.shape[0] - rank))
            print('LU  kmatrix.shape =', self.kmatrix.shape)

        print('\nMatrix value check\n******************')

        print('\nKzeromatrix:')
        sm, bg = self.MatrixValueCompare(self.kzeromatrix)
        print('Smallest value:  ', abs(sm))
        print('Largest value:   ', abs(bg))
        print('Ratio abs(bg/sm):', abs(bg / sm))

        print('\nLzeromatrix:')
        sm, bg = self.MatrixValueCompare(self.lzeromatrix)
        print('Smallest value:  ', abs(sm))
        print('Largest value:   ', abs(bg))
        print('Ratio abs(bg/sm):', abs(bg / sm))

        print(
            '\nPlease note: I\'ve found that for larger models the rank calculated by SVD in this test can become unstable and give incorrect results.'
        )
        print(' ')

        if resultback:
            return u, s, vh
