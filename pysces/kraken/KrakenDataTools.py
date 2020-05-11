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

import numpy, scipy, scipy.io

class DataTools:
    """
    Data manipulation tools, requires scipy, scipy.io
    """

    def data_filter_func(row):
        print('You should overwrite this (data_filter_func) method with your own or supply a function to data_filter')
        return False

    def DataFilter(self, data, filter_func=None, report=True):
        """Takes in an array and a function which returns a boolean for filter==True"""

        if filter_func == None:
            filter_func = self.data_filter_func

        data_out = numpy.zeros(data.shape,'d')
        rr = 0
        filtered = []
        for r in range(data_out.shape[0]):
            if filter_func(data[r,:]):
                filtered.append('Filtered: %s' % data[r,:])
            else:
                data_out[rr,:] = data[r,:]
                rr += 1
        # 'remove filtered rows'
        data_out = data_out[:rr,:]
        if report:
            for a in filtered:
                print(a)
        return data_out

    def GridSortLR(self, arr):
        """
        Sort irregular unstructured grid data into monotonically rising grid data
        (where possible) Work from left column to right (uses ultra cool NumPy functions :-)
        """
        print('arr.shape:', arr.shape)
        sorted = numpy.rot90(arr)
        print('rotated.shape:', sorted.shape)
        sorted_idx = numpy.lexsort(sorted)
        print('sorted_idx.shape', sorted_idx.shape)
        sorted = sorted.take(sorted_idx, axis=-1)
        sorted = numpy.flipud(sorted)
        arr = numpy.transpose(sorted)
        print('arr.shape:', arr.shape)
        return arr


    def writeGnuPlot_3D(self, arr, fname):
        if arr.shape[0] > 5e6:
            print("warning: writeGnuPlot_3D does not work for large arrays")
        try:
            scipy.io.write_array(fname+'_result_gplt.dat', arr, precision=9, separator=' ')
        except Exception as ex:
            print('writeGnuPlot_3D exception raised')
            print(ex)

    def vtk_scalar_func(self, val):
        """
        Manipulate the scalar (or z-axis) values
        Accepts and returns a double
        """
        return val

    def writeVTK_UnstructuredGrid(self, arr, fname, scalar_func=None):
        """
        Write arr as a VTK unstructured grid to file "fname.vtk"
        arr must either be in:
            x,y,z format (in which case z is used for scalar data) or
            z,y,z,s format (s being point scalar value)
        Scalar_func(val) function to manipulate the scalar (z-axis) values

        Original version http://mayavi.sourceforge.net/mwiki/PointClouds
        Hacked quite a bit by brett 20070221
        """
        assert arr.shape[1] == 3 or arr.shape[1] == 4, '\nneed 3 or 4 columns for this'
        if scalar_func == None:
            scalar_func = self.vtk_scalar_func
        if arr.shape[1] == 4:
            HAVE_SCALARS = 1
        else:
            HAVE_SCALARS = 0
            print('No scalar values supplied Z axis values will be used')

        n=arr.shape[0]
        print("n:",n)
        # write data to vtk polydata file
        # write header
        out = open(fname+'.vtk', 'w')
        h1 = "# vtk DataFile Version 2.0\n"
        h1 += "%s\n" % fname
        h1 += "ASCII\n"
        h1 += "DATASET UNSTRUCTURED_GRID\n"
        h1 += "POINTS " + str(n) + " double\n"
        out.write(h1)
        # write xyz data
        for r in range(n):
            #s = '%15.2f %15.2f %15.2f' % (x[i], y[i], z[i])
            out.write(str(arr[r,0])+" "+str(arr[r,1])+" "+str(arr[r,2])+'\n')

        # write cell data
        out.write("CELLS "+ str(n)+ " "+ str(2*n)+'\n')
        for r in range(n):
                #s = '1 %d \n' % (i)
                out.write("1 "+str(r)+"\n")

        # write cell types
        out.write("CELL_TYPES " + str(n)+'\n')
        for r in range(n):
            out.write("1 \n")

        # write z scalar values
        h2 = '\n' + """POINT_DATA """ + str(n) + "\n"
        h3 = "SCALARS %s double 1\n" % fname
        h3 += "LOOKUP_TABLE default\n"
        out.write(h2 + h3)

        for r in range(n):
            if HAVE_SCALARS:
                sc=(scalar_func(arr[r,3]))
            else:
                sc=(scalar_func(arr[r,2]))
            out.write(str(sc)+ "\n")

        out.write('\n')
        out.close()

class IntersectionAnalysis(DataTools):
    zmin = 0.01
    intersect_tol = 0.04
    coords = None
    vtkcoords = None
    intersection = None
    vtkintersection = None

    def setCoords(self, coords):
        self.coords = coords

    def setVTKCoords(self):
        self.vtkcoords = numpy.zeros((self.coords.shape[0], self.coords.shape[1]+1),'d')
        self.vtkcoords[:,:self.coords.shape[1]] = self.coords.copy()

    def data_filter_func(self, row):
        val = False
        ##  if (numpy.abs(row) < 1.0e-10).any() or (row < 0.0).any():
        if (row < self.zmin).any():
            val = True
        return val

    def CalculateIntersect(self, tol=None):
        """calc and add scalar value for intersection to vtkcoords"""
        if tol != None:
            self.intersect_tol = tol
        intersect = []
        cntr = 0
        for r in range(self.coords.shape[0]):
            if numpy.abs(self.coords[r, 2] - self.coords[r, 3]) < self.intersect_tol:
                intersect.append(self.coords[r,:])
                cntr += 1
        print('\n*****\n%s intersections found.\n*****\n' % cntr)
        self.intersection = numpy.array(intersect)
        self.vtkintersection = self.intersection.copy()
        del intersect

    def LogVTKdata(self):
        print("Logging VTK coords")
        # log for vtk
        for r in range(self.vtkcoords.shape[0]):
            self.vtkcoords[r,:self.coords.shape[1]] = scipy.log10(self.vtkcoords[r,:self.coords.shape[1]])

    def LogVTKIntersectiondata(self):
        print("Logging VTK intersection")
        # log for vtk intersect
        for r in range(self.vtkintersection.shape[0]):
            self.vtkintersection[r,:] = scipy.log10(self.intersection[r,:])
