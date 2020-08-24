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

__doc__ = '''PySCeS plotting module'''

import numpy, scipy
import os
import itertools
from time import strftime
from getpass import getuser


class PyscesGPlot:
    '''Prototype plotting class which wraps gnuplot functions and adds some stuff'''

    __version__ = __version__
    save_html_header = 1
    save_html_footer = 1
    mode_gnuplot4 = 1

    def __init__(self):
        import scipy.sandbox.gplt as gnuplot

        self.gnuplot = gnuplot

    def plot2D(
        self, ginput, x, ylist, cmdout=0, name=None, fmt='w l', ykey=None, log=0
    ):
        """
        plot2D(ginput,x,ylist,cmdout=0,name=None,fmt='w l',ykey=None,log=0)

        Plot a 2D graph using gnuplot

        Arguments:
        =========
        ginput: an array containing data to be plotted
        x: the index of the x-axis
        ylist: a list of y-axis indices
        cmdout [default=0]: do not plot but instead output the equivalent gnuplot command
        name [default=None]: the array name to be used with cmdout=1 (otherwise ginput)
        fmt [default='w l']: gnuplot line format string 'w l'=lines, 'w p'=points
        ykey [default=None]: list of y-axis key names (otherwise line index is used)
        log [default=0]: set x and y axis to linear (0) or log (1) scales

        """
        assert ginput.shape[0] >= 1, 'ginputRow must be >= 1'
        assert ginput.shape[1] > 1, 'ginputCol must be > 1'
        assert x <= ginput.shape[1], 'X index should be < ginputCol index'
        assert len(ylist) >= 1, 'Plotlist must contain a value'

        if ykey != None and len(ylist) != len(ykey):
            ykey = None

        pltcmd = 'self.gnuplot.plot('
        for dep in range(len(ylist)):
            pltcmd += (
                'ginput[:,' + str(int(x)) + '],ginput[:,' + str(int(ylist[dep])) + '],'
            )
            if ykey != None:
                tstr = str(ykey[dep])
            else:
                tstr = str(ylist[dep])
            pltcmd += '\'t \"' + tstr + '\" ' + fmt + '\','
        pltcmd += ')'
        # print pltcmd
        if cmdout and name != None:
            return pltcmd.replace('ginput', str(name))
        elif cmdout:
            return pltcmd
        else:
            eval(pltcmd)
            if log:
                self.gnuplot.logx('on')
                self.gnuplot.logy('on')
            else:
                self.gnuplot.logx('off')
                self.gnuplot.logy('off')

    def plot3D(
        self,
        ginput,
        x,
        y,
        z,
        cmdout=0,
        name=None,
        fmt='w l',
        zkey=None,
        arrayout=0,
        log=0,
    ):
        """
        plot3D(ginput,x,y,z,cmdout=0,name=None,fmt='w l',zkey=None,arrayout=0,log=0)

        Plot a 3D surface with gnuplot

        Arguments:
        =========
        ginput: an array containing data to be plotted
        x: x-axis index
        y: y-axis index
        z: z-axis index
        cmdout [default=0]: do not plot but instead output the equivalent gnuplot command
        name [default=None]: the array name to be used with cmdout=1 (otherwise ginput)
        fmt [default='w l']: gnuplot line format string 'w l'=lines, 'w p'=points
        zkey [default=None]: list of z-axis key names (otherwise line index is used)
        arrayout [default=0]: output the 3 column array of data used by scipy.plot3d
        log [default=0]: set x, y, z axis to linear (0) or log (1) scales

        """
        assert ginput.shape[0] >= 1, 'ginputRow must be >= 1'
        assert ginput.shape[1] > 1, 'ginputCol must be > 1'
        assert x <= ginput.shape[1], 'X index should be < ginputCol index'
        assert y <= ginput.shape[1], 'Y index should be < ginputCol index'
        assert z <= ginput.shape[1], 'Z index should be < ginputCol index'

        if zkey != None:
            tstr = str(zkey)
        else:
            tstr = str(z)

        pltcmd = 'self.gnuplot.plot3d('
        pltcmd += 'ginput[:,:3],'
        pltcmd += '\'t \"' + tstr + '\" ' + fmt + '\','
        pltcmd += ')'

        plotfile = []
        plotfile.append(ginput[:, x])
        plotfile.append(ginput[:, y])
        plotfile.append(ginput[:, z])
        ginput = scipy.transpose(scipy.array(plotfile))

        if cmdout and name != None:
            if arrayout:
                return pltcmd.replace('ginput', str(name)), ginput
            else:
                return pltcmd.replace('ginput', str(name))
        elif cmdout:
            if arrayout:
                return pltcmd, ginput
            else:
                return pltcmd
        else:
            eval(pltcmd)
            if log:
                self.gnuplot.logx('on')
                self.gnuplot.logy('on')
            else:
                self.gnuplot.logx('off')
                self.gnuplot.logy('off')

    def plotX(self, ginput, cmdout=0, name=None, fmt='w l', log=0):
        """
        plotX(ginput,cmdout=0,name=None,fmt='w l',log=0)

        Quick and dirty 2D plot where all other columns are plotted against the first column

        Arguments:
        =========
        ginput: an array containing data to be plotted
        cmdout [default=0]: do not plot but instead output the equivalent gnuplot command
        name [default=None]: the array name to be used with cmdout=1 (otherwise ginput)
        fmt [default='w l']: gnuplot line format string 'w l'=lines, 'w p'=points
        log [default=0]: set x and y axis to linear (0) or log (1) scales

        """
        assert (
            ginput.shape[0] >= 1 and ginput.shape[1] > 1
        ), 'array must have at least dim(1,2)'

        row, col = ginput.shape
        pltcmd = 'self.gnuplot.plot('
        for x in range(1, col):
            pltcmd += 'ginput[:,0],ginput[:,' + str(x) + '],'
            pltcmd += '\'t \"' + str(x) + '\" ' + fmt + '\''
            if x < col - 1:
                pltcmd += ','
        pltcmd += ')'

        if cmdout and name != None:
            return pltcmd.replace('ginput', str(name))
        elif cmdout:
            return pltcmd
        else:
            eval(pltcmd)
            if log:
                self.gnuplot.logx('on')
                self.gnuplot.logy('on')
            else:
                self.gnuplot.logx('off')
                self.gnuplot.logy('off')

    def __save_command__(self, FileName):
        """
        __save_command__(FileName)

        Selects the correct gnuplot command depending on gnuplot version set with mode_gnuplot4

        Arguments:
        =========
        FileName: the name of the resulting output fil

        """
        if self.mode_gnuplot4:
            self.gnuplot.output(
                FileName,
                'png',
                options='medium size 640,480 xffffff x000000 x404040\
                              xff0000 xffa500 x66cdaa xcdb5cd\
                              xadd8e6 x0000ff xdda0dd x9500d3',
            )
        else:
            self.gnuplot.save(FileName)

    def save(self, filename, path=None):
        """
        save(filename,path=None)

        Save the current active plot to a file

        Arguments:
        =========
        filename: the name of the output file
        path [default=None]: the directory where the file should be saved

        """
        try:
            if filename[-4:] != '.png':
                filename += '.png'
                print('File saved as: ' + filename)
        except:
            pass

        if path != None:
            filename = os.path.join(path, filename)
        else:
            filename = os.path.join(os.getcwd(), filename)
        self.__save_command__(filename)

    def save_html(self, imagename, File=None, path=None, name=None, close_file=1):
        """
        save_html(imagename,File=None,path=None,name=None,close_file=1)

        Save the current plot to a file plus an HTML file that includes the resulting image

        Arguments:
        =========
        imagename: the output image name
        File [default=None]: the HTML filename defaults to imagename
        path [default=None]: the output directory
        name [default=None]: the HTML title of the image
        close_file [default=1]: close the HTML (1) or leave it open for further editing (0)

        """
        imagename = str(imagename)
        if name != None:
            name = str(name)
        else:
            name = ''
        try:
            if imagename[-4:] != '.png':
                imagename += '.png'
                print('File saved as: ' + imagename)
        except:
            pass

        if path != None:
            imagenameout = os.path.join(path, imagename)
            if File == None:
                File = open(os.path.join(path, imagename[:-3] + 'html'), 'w')
        else:
            imagenameout = os.path.join(os.getcwd(), imagename)
            if File == None:
                File = open(os.path.join(os.getcwd(), imagename[:-3] + 'html'), 'w')
        self.__save_command__(imagename)

        fname = (
            'PySCeS generated image - '
            + imagename
            + '" generated from model file: '
            + strftime("%H:%M:%S")
        )

        if File != None:
            assert (
                File != file
            ), 'WriteArray(input,File=None,Row=None,Col=None,close_file=0)'
            header = '\n'
            header += (
                '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">\n'
            )
            header += '<html>\n'
            header += '<head>\n'
            header += (
                '<title>PySCeS generated image - '
                + imagename
                + ' - '
                + strftime("%H:%M:%S (%Z)")
                + '</title>\n'
            )
            header += '<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1">\n'
            header += '</head>\n'
            header += '<body bgcolor="#FFFFCC">\n\n'
            header += (
                '<h4><a href="http://pysces.sourceforge.net">PySCeS</a> generated image - '
                + imagename
                + '</h4>\n\n'
            )

            if self.save_html_header:
                File.write(header)
            File.write('\n<!-- ' + imagename + '-->\n\n')

            File.write('<p align="center">\n')

            File.write(
                '<img src="' + imagename + '" width="640" height="480" border="2">\n'
            )

            if name != None:
                File.write('<br><font size="3">' + str(name) + '</font>\n')

            File.write('</p>\n\n')
            if self.save_html_footer:
                try:
                    File.write(
                        '<p><a href="http://pysces.sourceforge.net"><font size="3">PySCeS '
                        + self.__version__
                        + '</font></a><font size="2"> HTML output (image <i>'
                        + imagename
                        + '</i> produced at '
                        + strftime("%H:%M:%S")
                        + ' by <i>'
                        + getuser()
                        + '</i>)</font></p>\n'
                    )
                except:
                    File.write(
                        '<p><a href="http://pysces.sourceforge.net"><font size="3">PySCeS '
                        + __version__
                        + '</font></a><font size="2"> HTML output (model <i>'
                        + imagename
                        + '</i> produced at '
                        + strftime("%H:%M:%S - %Z")
                        + ')</font></p>\n'
                    )
                File.write('</body>\n')
                File.write('</html>\n')
            else:
                File.write('\n')

            if close_file:
                File.close()

    def gridon(self):
        """
        gridon()

        Enable grid for the current active plot

        Arguments:
        None

        """
        self.gnuplot.grid('on')

    def gridoff(self):
        """
        gridoff()

        Disable grid for the current active plot

        Arguments:
        None

        """
        self.gnuplot.grid('off')

    def logx(self):
        """
        logx()

        Set x logscale for active plot

        Arguments:
        None

        """
        self.gnuplot.logx('on')

    def logy(self):
        """
        logy()

        Set y logscale for active plot

        Arguments:
        None

        """
        self.gnuplot.logy('on')

    def logxy(self):
        """
        logxy()

        Set x and y logscale for active plot

        Arguments:
        None

        """
        self.gnuplot.logx('on')
        self.gnuplot.logy('on')

    def linx(self):
        """
        linx()

        Set x linear scale for active plot

        Arguments:
        None

        """
        self.gnuplot.logx('off')

    def liny(self):
        """
        liny()

        Set y linear scale for active plot

        Arguments:
        None

        """
        self.gnuplot.logy('off')

    def linxy(self):
        """
        linxy()

        Set x and y linear scale for active plot

        Arguments:
        None

        """
        self.gnuplot.logx('off')
        self.gnuplot.logy('off')

    def logxliny(self):
        """
        logxliny()

        Set logscale x and linear scale y for currently active plot

        Arguments:
        None

        """
        self.gnuplot.logx('on')
        self.gnuplot.logy('off')

    def linxlogy(self):
        """
        linxlogy()

        Set logscale y and linear scale x for currently active plot

        Arguments:
        None

        """
        self.gnuplot.logx('off')
        self.gnuplot.logy('on')

    def xrng(self, start, end):
        """
        xrng(start,end)

        Set the range of the x-axis for currently active plot

        Arguments:
        =========
        start: lower value
        end: upper value

        """
        self.gnuplot.xaxis((float(start), float(end)))

    def yrng(self, start, end):
        """
        yrng(start,end)

        Set the range of the y-axis for currently active plot

        Arguments:
        =========
        start: lower value
        end: upper value

        """
        self.gnuplot.yaxis((float(start), float(end)))

    def zrng(self, start, end):
        """
        zrng(start,end)

        Set the range of the z-axis for currently active plot

        Arguments:
        =========
        start: lower value
        end: upper value

        """
        self.gnuplot.zaxis((float(start), float(end)))

    def xlabel(self, l=''):
        """
        xlabel(l='')

        Set the x-axis label for the currently active plot

        Arguments:
        =========
        l [default=' ']: axis label string

        """
        self.gnuplot.xtitle(str(l))

    def ylabel(self, l=''):
        """
        ylabel(l='')

        Set the y-axis label for the currently active plot

        Arguments:
        =========
        l [default=' ']: axis label string

        """
        self.gnuplot.ytitle(str(l))

    def zlabel(self, l=''):
        """
        zlabel(l='')

        Set the z-axis label for the currently active plot

        Arguments:
        =========
        l [default=' ']: axis label string

        """
        self.gnuplot.ztitle(str(l))

    def title(self, l=''):
        """
        title(l='')

        Set the graph title for the currently active plot

        Arguments:
        =========
        l [default=' ']: graph label string

        """
        self.gnuplot.title(str(l))

    def ticslevel(self, x=0.0):
        """
        ticslevel(x=0.0)

        Set the gnuplot ticlevel for the currently active plot

        Arguments:
        =========
        x [default=0.0]: gnuplot ticlevel

        """
        self.gnuplot.ticlevel(float(x))


class LineObj(object):
    __id__ = None
    idx = None
    data = None
    prop = None

    def __init__(self, data, idx, **kwargs):
        self.idx = idx
        self.data = data
        self.prop = {
            'label': None,
            'linewidth': None,
            'colour': None,
            'style': None,
        }
        for k in list(kwargs.keys()):
            if k in self.prop:
                self.prop[k] = kwargs[k]


class AxisObj(object):
    __id__ = None
    idx = None
    data = None
    prop = None

    def __init__(self, data, idx, **kwargs):
        self.idx = idx
        self.data = data
        self.prop = {
            'label': None,
            'log': False,
        }
        for k in list(kwargs.keys()):
            if k in self.prop:
                self.prop[k] = kwargs[k]


class GraphicsObj(object):
    xaxis = None
    yaxis = None
    data = None
    data_shape = None
    data_type = None
    data_lines = None
    var_idx = None
    scalars = None

    def setData(self, data):
        assert type(data) == numpy.ndarray, '\nI need a numpy/scipy ndarray!'
        assert data.shape[1] > 1, '\nI need at least a 2D array!'
        self.data = data
        self.data_shape = data.shape
        self.data_type = data.dtype.char
        self.scalars = numpy.zeros((data.shape[0], 1), 'd')
        self.resetAxis()
        self.setDataLabels(['l' + str(l) for l in range(self.data_shape[1])])

    def setDataLabels(self, labels):
        assert len(labels) == self.data_shape[1], '\nList unacceptable length.'
        if self.data_lines != None:
            [self.__delattr__(n) for n in self.data_lines]
        self.data_lines = []
        for l in range(len(labels)):
            self.data_lines.append(labels[l])
            ##  lo = LineObj(self.data[:,l].reshape(self.data_shape[0], 1), l, label=labels[l])
            lo = LineObj(numpy.take(self.data, [l], axis=1), l, label=labels[l])
            setattr(self, labels[l], lo)

    def setAxis(self, **kwargs):
        ##  print kwargs
        for key in list(kwargs.keys()):
            if key in ['x', 'y'] and kwargs[key] in self.var_idx:
                if getattr(self, key + 'axis') != None:
                    self.var_idx.append(getattr(self, key + 'axis').idx)
                    self.var_idx.sort()
                    ##  print '\t', self.var_idx
                ##  ao = AxisObj(self.data[:,kwargs[key]].reshape(self.data_shape[0], 1), kwargs[key], label=key)
                ao = AxisObj(
                    numpy.take(self.data, [kwargs[key]], axis=1), kwargs[key], label=key
                )
                setattr(self, key + 'axis', ao)
                a = self.var_idx.pop(self.var_idx.index(kwargs[key]))
                ##  print '\tkey idx :', key, kwargs[key]
        ##  print 'var_idx', self.var_idx

    def resetAxis(self):
        self.xaxis = self.yaxis = None
        self.var_idx = list(range(self.data_shape[1]))

    def getIdx(self, label):
        assert (
            label in self.data_lines
        ), '\n%s is not a valid data line.\nSelect one of %s' % (label, self.data_lines)
        return self.data_lines.index(label)

    def setData2Scalar(self, idx):
        self.scalars = self.data.take([idx], axis=1)

    def dataMinMax(self, idx):
        return (numpy.min(self.data[:, idx]), numpy.max(self.data[:, idx]))

    def dataLog10(self, idx):
        self.data[:, idx] = numpy.log10(self.data[:, idx])

    def getVarData(self):
        return self.data.take(self.var_idx, axis=1)

    def getRawData(self):
        return {
            'data': self.data,
            'scalars': self.scalars,
            'x_idx': self.xaxis.idx,
            'y_idx': self.yaxis.idx,
            'var_idx': self.var_idx,
            'data_lines': self.data_lines,
        }


class Graphics2dObj(GraphicsObj):
    title = 'Graph'
    xtitle = 'x'
    ytitle = 'y'
    xlog = False
    ylog = False
    selected = None
    __dims__ = 2

    def __init__(self, data, xidx=None, labels=None, selected=None):
        self.setData(data)
        if xidx != None:
            self.setAxis(x=xidx)
        if labels != None:
            self.setDataLabels(labels)
        if selected != None:
            self.selected = selected
        else:
            selected = []

    # working on this ... think I'm almost finished now
    def getNewByIndex(self, x, lines):
        data = self.data.take([x] + lines, axis=1)
        G = Graphics2dObj(data, 0)
        G.setDataLabels([self.data_lines[x]] + [self.data_lines[i] for i in lines])
        return G

    # working on this ... think I'm almost finished now
    def getNewByName(self, x, lines):
        return self.getNewByIndex(self.getIdx(x), [self.getIdx(l) for l in lines])


class Graphics3dObj(Graphics2dObj):
    ztitle = 'z'
    zlog = False
    __dims__ = 3

    def __init__(self, data, xidx, yidx):
        self.setData(data)
        self.setAxis(x=xidx, y=yidx)


class PyscesGPlot2MPL:
    """
    Old 'gnuplot' plotting functions using a matplotlib backend
    """

    mode_interactive = True
    mode_create_new_figure = True
    new_figure_generator = None
    current_figure = None
    mode_hold = False
    # backwards compatibility
    save_html_header = 1
    save_html_footer = 1
    max_open_windows = 50

    def __init__(self):
        import matplotlib

        try:
            matplotlib.use('TKagg')
        except Exception as ex:
            print(ex)
            print(
                "\nPySCeS uses matplotlib's TKagg backend for interactive plotting please enable this backend when compiling matplotlib"
            )
        ## import matplotlib.axes3d
        import pylab

        ## self._matplotlib_axes3d_ = matplotlib.axes3d
        self.P = pylab
        if self.mode_interactive:
            self.P.ion()
        self.setNewFigureGenerator(self.max_open_windows)

    def setNewFigureGenerator(self, num):
        """
        Create a new_figure_generator with range [1,num+1]
        """
        self.new_figure_generator = itertools.cycle(list(range(1, num + 1)))
        self.max_open_windows = num

    def setActiveFigure(self, fignum):
        self.current_figure = self.P.figure(fignum)

    def setNewFigure(self):
        self.setActiveFigure(next(self.new_figure_generator))

    def closeAllPlotWindows(self):
        for x in range(self.max_open_windows):
            self.P.close()

    def plot2D(self, data, x, ylist=[], labels=None, style=None):
        assert type(ylist) == list, '\nylist must be a list'
        if self.mode_create_new_figure:
            self.setNewFigure()
        if not self.mode_hold:
            self.P.clf()
        D = Graphics2dObj(data, x)
        if labels != None:
            D.setDataLabels(labels)
        if len(ylist) == 0:
            print('ylist empty plotting all data vs (%s)' % x)
            ylist = list(range(data.shape[1]))
            ylist.pop(ylist.index(x))
        if len(ylist) > 0:
            ylt = numpy.array(ylist) < D.data_shape[1]
            assert (
                type(ylt) == numpy.ndarray
            ), '\nInvalid y list this can be caused by a forgotten keyword labels='
            assert ylt.all(), '\nInvalid y index'
            D = D.getNewByIndex(x, ylist)
        self.P.ioff()
        usedline = []
        for line in D.data_lines:
            L = getattr(D, line)
            ##  print line, D.xaxis.idx, L.idx
            if L.idx != D.xaxis.idx and line not in usedline:
                ##  print 'done'
                usedline.append(line)
                if style != None:
                    self.P.plot(D.xaxis.data, L.data, style, label=L.prop['label'])
                else:
                    self.P.plot(D.xaxis.data, L.data, label=L.prop['label'])
        if self.mode_interactive:
            self.P.ion()
        self.P.legend()

    def plot3D(self, *args, **kwargs):
        print("*****\nplot3D not implemented yet.\n*****")

    def plotX(self, data):
        if self.mode_create_new_figure:
            self.setNewFigure()
        if not self.mode_hold:
            self.P.clf()
        D = Graphics2dObj(data, None)
        ##  print type(D)
        ##  print D.var_idx
        ##  print D.data[1]
        self.P.plot(D.getVarData())

    def save(self, fname, path=None):
        fname += '.png'
        if path != None:
            fname = os.path.join(path, fname)
        self.P.savefig(fname)

    def legendOff(self):
        L = self.P.gca().get_legend()
        L._visible = False
        self.P.draw()

    def gridon(self):
        self.P.grid(True)

    def gridoff(self):
        self.P.grid(False)

    def logx(self):
        A = self.P.gca()
        A.set_xscale('log')
        self.P.draw()

    def logy(self):
        A = self.P.gca()
        A.set_yscale('log')
        if self.mode_interactive:
            self.P.draw()

    def linx(self):
        A = self.P.gca()
        A.set_xscale('linear')
        if self.mode_interactive:
            self.P.draw()

    def liny(self):
        A = self.P.gca()
        A.set_yscale('linear')
        if self.mode_interactive:
            self.P.draw()

    def linxy(self):
        A = self.P.gca()
        A.set_xscale('linear')
        A.set_yscale('linear')
        if self.mode_interactive:
            self.P.draw()

    def logxy(self):
        A = self.P.gca()
        A.set_xscale('log')
        A.set_yscale('log')
        if self.mode_interactive:
            self.P.draw()

    def logxliny(self):
        A = self.P.gca()
        A.set_xscale('log')
        A.set_yscale('linear')
        if self.mode_interactive:
            self.P.draw()

    def linxlogy(self):
        A = self.P.gca()
        A.set_xscale('linear')
        A.set_yscale('log')
        if self.mode_interactive:
            self.P.draw()

    def xrng(self, start, end):
        self.P.xlim(start, end)

    def yrng(self, start, end):
        self.P.ylim(start, end)

    def zrng(self, start, end):
        print("*****\n\tNot implemented yet.\n*****")

    def xlabel(self, l=''):
        self.P.xlabel(l)

    def ylabel(self, l=''):
        self.P.ylabel(l)

    def zlabel(self, l=''):
        print("*****\n\tNot implemented yet.\n*****")

    def title(self, l=''):
        self.P.title(l)

    def ticslevel(self, x=0.0):
        print("*****\n\tNot implemented yet.\n*****")

    def setLineWidth(self, w=1):
        'set line width for current axis'
        A = self.P.gca()
        [l.set_linewidth(w) for l in A.get_lines()]
        if self.mode_interactive:
            self.P.draw()

    def save_html(self, imagename, File=None, path=None, name=None, close_file=1):
        imagename = str(imagename)
        if name != None:
            name = str(name)
        else:
            name = ''
        try:
            if imagename[-4:] != '.png':
                imagename += '.png'
                print('File saved as: ' + imagename)
        except:
            pass

        if path != None:
            imagenameout = os.path.join(path, imagename)
            if File == None:
                File = open(os.path.join(path, imagename[:-3] + 'html'), 'w')
        else:
            imagenameout = os.path.join(os.getcwd(), imagename)
            if File == None:
                File = open(os.path.join(os.getcwd(), imagename[:-3] + 'html'), 'w')
        self.P.savefig(imagenameout, dpi=80)

        fname = (
            'PySCeS generated image - '
            + imagename
            + '" generated from model file: '
            + strftime("%H:%M:%S")
        )

        if File != None:
            assert (
                File != file
            ), 'WriteArray(input,File=None,Row=None,Col=None,close_file=0)'
            header = '\n'
            header += (
                '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">\n'
            )
            header += '<html>\n'
            header += '<head>\n'
            header += (
                '<title>PySCeS generated image - '
                + imagename
                + ' - '
                + strftime("%H:%M:%S (%Z)")
                + '</title>\n'
            )
            header += '<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1">\n'
            header += '</head>\n'
            header += '<body bgcolor="#FFFFCC">\n\n'
            header += (
                '<h4><a href="http://pysces.sourceforge.net">PySCeS</a> generated image - '
                + imagename
                + '</h4>\n\n'
            )

            if self.save_html_header:
                File.write(header)
            File.write('\n<!-- ' + imagename + '-->\n\n')

            File.write('<p align="center">\n')

            File.write(
                '<img src="' + imagename + '" width="640" height="480" border="2">\n'
            )

            if name != None:
                File.write('<br><font size="3">' + str(name) + '</font>\n')

            File.write('</p>\n\n')
            if self.save_html_footer:
                try:
                    File.write(
                        '<p><a href="http://pysces.sourceforge.net"><font size="3">PySCeS '
                        + self.__version__
                        + '</font></a><font size="2"> HTML output (image <i>'
                        + imagename
                        + '</i> produced at '
                        + strftime("%H:%M:%S")
                        + ' by <i>'
                        + getuser()
                        + '</i>)</font></p>\n'
                    )
                except:
                    File.write(
                        '<p><a href="http://pysces.sourceforge.net"><font size="3">PySCeS '
                        + __version__
                        + '</font></a><font size="2"> HTML output (model <i>'
                        + imagename
                        + '</i> produced at '
                        + strftime("%H:%M:%S - %Z")
                        + ')</font></p>\n'
                    )
                File.write('</body>\n')
                File.write('</html>\n')
            else:
                File.write('\n')

            if close_file:
                File.close()


class PyscesGPlot2MplExt(PyscesGPlot2MPL):
    def __init__(self):
        PyscesGPlot2MPL.__init__(self)

    def plotLines2DX(self, data, x, labels=None):
        if not self.mode_hold:
            self.P.clf()
        ##  F = self.P.figure(1)
        ##  S = self.P.subplot(1,1,1)
        D = Graphics2dObj(data, x)
        print(D.var_idx)
        print(D.data[1])
        self.P.plot(D.xaxis.data, D.getVarData())
