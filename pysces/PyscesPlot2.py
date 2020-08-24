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

from .version import __version__
from . import __SILENT_START__
import subprocess, os, shutil, time, math, itertools, copy
from getpass import getuser
import numpy, scipy

try:
    input = raw_input  # Py2 compatibility
except NameError:
    pass

__doc__ = '''
            PyscesPlot2
            -----------

            PyscesPlot2 is a new graphics susbsystem for PySCeS which will include a
            Unified Plotting Interface which can take advantage of different plotting
            backends via a common user interface.
            '''


class PlotBase(object):
    """
    Abstract class defining the Unified Plotting Interface methods. These methods should
    be overridden and the class extended by interface specific subclasses.
    """

    CommonStyleDefs = {'points': '', 'lines': ''}

    def plot(self, data, x, y, title='', format=''):
        """
        Plot a single line data[y] vs data[x] where:

        - *data* the data array
        - *x* x column index
        - *y* y column index
        - *title* is the line key
        - *format* is the XXX format string (default='')

        Format can also be the *CommonStyle* 'lines' or 'points'
        """
        raise NotImplementedError

    def plotLines(self, data, x, y=[], titles=[], formats=['']):
        """
        Plot a multiple lines data[y1, y2, ] vs data[x] where:

        - *data* the data array
        - *x* x column index
        - *y* is a list of line indexes, if empty all of y not including x is plotted
        - *titles* is a list of line keys, if empty Line1,Line2,Line3 is used
        - *formats* is a list (per line) of XXX format strings.

        If *formats* only contains a single item, this format is used for all lines
        and can also be the *CommonStyle* 'lines' or 'points'.
        """
        raise NotImplementedError

    def splot(self, data, x, y, z, titles='', format=''):
        """
        Plot a surface data[z] vs data[y] vs data[x] where:

        - *data* the data array
        - *x* x column index
        - *y* y column index
        - *z* z column index
        - *title* is the surface key
        - *format* is the XXX format string (default='')

        Format can also be the *CommonStyle* 'lines' or 'points'.
        """
        raise NotImplementedError

    def splotSurfaces(self, data, x, y, z=[], titles=[], formats=['']):
        """
        Plot data[z1, z2, ] vs data[y] vs data[x] where:

        - *data* the data array
        - *x* x column index
        - *y* y column index
        - *z* list of z column indexes, if empty all of z not including x, y are plotted
        - *titles* is a list of surface keys, if empty Surf1, Surf2, Surf3 is used
        - *formats* is a list (per line) of XXX format strings (default='').

        If *formats* only contains a single item, this format is used for all surfaces
        and can also be the *CommonStyle* 'lines' or 'points'.
        """
        raise NotImplementedError

    def setLogScale(self, axis):
        """
        Set axis to logscale where:

        - *axis* = x, y, z, xy, xz, yz, zyx
        """
        raise NotImplementedError

    def setNoLogScale(self, axis):
        """
        Set axis to a linear scale where:

        - *axis* = x, y, z, xy, xz, yz, zyx
        """
        raise NotImplementedError

    def setRange(self, axis, min=None, max=None):
        """
        Set axis range where

        - *axis* = x, y, z, xy, xz, yz, zyx
        - *min* = range(s) lower bound (default=None) autoscale
        - *max* = range(s) upper bound (default=None) autoscale
        """
        raise NotImplementedError

    def setGrid(self, value):
        """
        Display or remove graph grid.

        - *value* (boolean) True (on) or False (off)
        """
        raise NotImplementedError

    def setGraphTitle(self, title='PySCeS Plot'):
        """
        Set the graph title, unset if title=None

        - *title* (string, default='PySCeS Plot') the graph title
        """
        raise NotImplementedError

    def setKey(self, value=False):
        """
        Enable or disable the current plot key, no arguments removes key.

        - *value* boolean (default = False)
        """
        raise NotImplementedError

    def setAxisLabel(self, axis, label=''):
        """
        Set the axis label:

        - *axis* = x, y, z, xy, xz, yz, zyx
        - *label* = string (default='')

        Called with only the axis argument clears the axis label.
        """
        raise NotImplementedError

    def save(self, name, directory=None, dfmt='%.8e'):
        """
        Save the plot data and (optionally) XXX format file

        - *filename* the filename
        - *directory* optional (default = current working directory)
        - *dfmt* the data format string (default='%.8e')
        """
        raise NotImplementedError

    def export(self, name, directory=None, type='png'):
        """
        Export the current plot as a <format> file.

        - *filename* the filename
        - *directory* optional (default = current working directory)
        - *type* the file format (default='png').

        Currently only PNG is guaranteed to be available in all interfaces.
        """
        raise NotImplementedError

    def axisInputStringToList(self, input):
        """
        Extracts axis information from a string input, returns a boolean triple
        representing (x=True/False, y=True/False, z=True/False).

        - *input* the input string
        """
        out = [None, None, None]
        input = input.lower()
        if 'x' in input:
            out[0] = 'x'
        if 'y' in input:
            out[1] = 'y'
        if 'z' in input:
            out[2] = 'z'
        return out

    def wait(self, seconds=3):
        """
        Wait *seconds* (default = 3) or until enter is pressed (seconds = -1)
        """
        if seconds == -1:
            input("\nPress <enter> to continue ...\n")
        else:
            time.sleep(seconds)


class GnuPlotUPI(PlotBase):
    """
    PySCeS/GnuPlot is reborn, leaner and meaner than ever before. This class enables
    plotting with GnuPlot via a subprocess link:

    - *work_dir* optional argument setting directory for dat file(s)
    - *gnuplot_dir* optional argument specifying the location of pgnuplot.exe (win32) or gnuplot

    GnuPlot backend to the Unified Plotting Interface.
    """

    __G_SUBPROC__ = None
    __ECHO__ = False
    __DATA_FILE_NAME__ = '_gnuplot.dat'
    __WORK_DIR__ = None
    __DATA_FILE_PATH__ = None
    DATF_FORMAT = '%.8e'
    __LAST_COMMAND__ = ''
    __MULTIPLOT__ = False
    PAUSE_TIME = 0.1
    __GNUPLOT_EXE_PATH__ = None
    __DISP_TERMINAL__ = None
    ##  MULTIPLOT_CNTR = 0
    ##  MULTIPLOT_MAX = 3
    ##  datFbuffsize = 5
    ##  datFcache = None

    CommonStyleDefs = {'points': 'w p', 'lines': 'w l'}

    Terminals = {
        'x11': '',
        'windows': '',
        'png': 'medium size 800,600'
        # 'xffffff x000000 x404040 xff0000 xffa500 x66cdaa xcdb5cd xadd8e6 x0000ff xdda0dd x9500d3'
    }

    def __init__(self, work_dir=None, gnuplot_dir=None):
        if gnuplot_dir != None or not os.path.exists(gnuplot_dir):
            GNUPLOT_PATH = gnuplot_dir
        else:
            GNUPLOT_PATH = ''
        if os.sys.platform == 'win32':
            self.__GNUPLOT_EXE_PATH__ = os.path.join(GNUPLOT_PATH, 'pgnuplot.exe')
            self.__DISP_TERMINAL__ = 'windows'
        else:
            self.__GNUPLOT_EXE_PATH__ = os.path.join(GNUPLOT_PATH, 'gnuplot')
            self.__DISP_TERMINAL__ = 'x11'
        try:
            self.__G_SUBPROC__ = subprocess.Popen(
                self.__GNUPLOT_EXE_PATH__, stdin=subprocess.PIPE
            )
        except Exception as ex:
            print('\nGnuPlot load failure\n')
            print(ex)
        if work_dir == None or not os.path.exists(work_dir):
            work_dir = os.getcwd()
        self.__DATA_FILE_PATH__ = os.path.join(work_dir, self.__DATA_FILE_NAME__)
        self.__WORK_DIR__ = work_dir

    def set(self, key, value=''):
        """
        Send *set <key>* or optionally *set <key> <value>* to GnuPlot.
        """
        self.g_write("set %s %s" % (key, value))

    def unset(self, key, value=''):
        """
        Send *unset <key>* or optionally *unset <key> <value>* to GnuPlot.
        """
        self.g_write("unset %s %s" % (key, value))

    def g_write(self, cmd):
        """
        Write a command to the GnuPlot interpreter

        - *cmd* the GnuPlot command
        """
        cmd = cmd.replace('\\', '/')
        self.last_command = cmd
        if self.__ECHO__:
            print('SENT: %s' % cmd)

        self.__G_SUBPROC__.stdin.write("""%s\n""" % cmd)

    def g_file_write_array(self, arr, dfmt=None):
        """
        Write a normal (2D) dataset to temp file. Dumps the array to file
        using the format:

        - *arr* the array (r>0, c>1)
        - *fmt* default '%.8e'
        """
        if dfmt == None:
            dfmt = self.DATF_FORMAT
        numpy.savetxt(self.__DATA_FILE_PATH__, arr, fmt=dfmt, delimiter=' ')
        if self.__ECHO__:
            print('WROTE: %s' % self.__DATA_FILE_PATH__)
        self.g_pause()  # multiplot protection

    def g_file_write_array3D(self, arr, yaxis=1, dfmt=None):
        """
        Write a GnuPlot format 3D dataset. The *yaxis* argument
        specifies the column that should be used to split the dataset
        into GnuPlot slices.

        - *arr* the array (r>1, c>2)
        - *fmt* default '%.8e'
        - *yaxis* default 1
        """
        outlist = []
        set_idx = 0
        curr_y = arr[0, yaxis]
        for r in range(arr.shape[0]):
            if arr[r, yaxis] != curr_y:
                outlist.append(arr[set_idx:r, :])
                curr_y = arr[r, yaxis]
                set_idx = r
        if set_idx < arr.shape[0]:
            outlist.append(arr[set_idx:, :])

        if dfmt == None:
            dfmt = self.DATF_FORMAT
        if len(outlist) <= 1:
            self.g_file_write_array(arr, dfmt=dfmt)
        else:
            F = open(self.__DATA_FILE_PATH__, 'w+')
            for d in range(len(outlist)):
                numpy.savetxt(F, outlist.pop(0), fmt=dfmt)
                F.write(' \n')
            F.flush()
            F.close()
            if self.__ECHO__:
                print('WROTE: %s' % self.__DATA_FILE_PATH__)
        del outlist
        self.g_pause()

    def g_pause(self):
        """
        A small pause defined by *self.PAUSE_TIME* (multiplied by 2 when in multiplot).
        """
        if self.__MULTIPLOT__:
            time.sleep(self.PAUSE_TIME * 2.0)
        else:
            time.sleep(self.PAUSE_TIME)

    def plot(self, data, x, y, title='', format='w l'):
        """
        Plot a single line data[y] vs data[x] where:

        - *data* the data array
        - *x* x column index
        - *y* y column index
        - *title* is the line key
        - *format* is the GnuPlot format string (default='w l')

        Format can also be the *CommonStyle* 'lines' or 'points'.
        """
        if format in list(self.CommonStyleDefs.keys()):
            format = self.CommonStyleDefs[format]
        if title == '':
            title = []
        else:
            title = ['X', title]
        if format == '':
            format = 'w l'
        self.plotLines(data, x, [y], title, [format] * data.shape[1])

    def plotLines(self, data, x, y=[], titles=[], formats=['w l']):
        """
        Plot a multiple lines data[y1, y2, ] vs data[x] where:

        - *data* the data array
        - *x* x column index
        - *y* is a list of line indexes, if empty all of y not including x is plotted
        - *titles* is a list of line keys if empty Line1, Line2, Line3 is used
        - *formats* is a list (per line) of GnuPlot format strings (default='w l').

        If *formats* only contains a single item, this format is used for all lines
        and can also be the *CommonStyle* 'lines' or 'points'.
        """

        self.g_file_write_array(data, dfmt=self.DATF_FORMAT)
        if len(y) == 0:
            y = list(range(data.shape[1]))
            y.pop(x)
        if len(formats) != 1 and len(formats) != data.shape[1]:
            print(
                "len(titles) must be one or equal the number data columns (%s) using: %s"
                % (data.shape[1], formats[-1])
            )
            formats = [formats[-1]]
        elif formats[0] == '' or formats[0] == None:
            formats = ['w l']
        if len(formats) == 1 and data.shape[1] > 1:
            formats = data.shape[1] * formats
        cFormats = list(self.CommonStyleDefs.keys())
        for f in range(len(formats)):
            if formats[f] in cFormats:
                formats[f] = self.CommonStyleDefs[formats[f]]
        if len(titles) == 0:
            titles = ['col%s' % (l) for l in range(data.shape[1])]

        cmd = '''plot '''
        for yi in y:
            cmd += '''"%s" u %s:%s t "%s" %s, ''' % (
                self.__DATA_FILE_PATH__,
                x + 1,
                yi + 1,
                titles[yi],
                formats[yi],
            )
        cmd = cmd[:-2]
        self.g_write(cmd)
        self.setKey(True)

    def splot(self, data, x, y, z, titles='', format='w l'):
        """
        Plot a surface data[z] vs data[y] vs data[x] where:

        - *data* the data array
        - *x* x column index
        - *y* y column index
        - *z* z column index
        - *titles* is the surface key
        - *format* is the GnuPlot format string (default='w l')

        Format can also be the *CommonStyle* 'lines' or 'points'.
        """
        if format in list(self.CommonStyleDefs.keys()):
            format = self.CommonStyleDefs[format]
        if titles == '':
            titles = []
        else:
            assert (
                len(titles) == data.shape[1]
            ), '\nTitle list must match number of columns'
        if format == '':
            format = 'w l'
        self.splotSurfaces(data, x, y, [z], titles, [format] * data.shape[1])

    def splotSurfaces(self, data, x, y, z=[], titles=[], formats=['w l']):
        """
        Plot data[z1, z2, ] vs data[y] vs data[x] where:

        - *data* the data array
        - *x* x column index
        - *y* y column index
        - *z* list of z column indexes, if empty all of z not including x, y are plotted
        - *titles* is a list of surface keys, if empty Surf1, Surf2, Surf3 is used
        - *formats* is a list (per line) of GnuPlot format strings (default='w l').

        If *formats* only contains a single item, this format is used for all surface and
        can also be the *CommonStyle* 'lines' or 'points'.
        """
        self.g_file_write_array3D(data, yaxis=y, dfmt=self.DATF_FORMAT)

        if len(z) == 0:
            z = list(range(data.shape[1]))
        if len(formats) != 1 and len(formats) != data.shape[1]:
            print(
                "len(titles) must be one or equal the number data columns (%s) using: %s"
                % (data.shape[1], formats[-1])
            )
            formats = [formats[-1]]
        elif formats[0] == '' or formats[0] == None:
            formats = ['w l']
        if len(formats) == 1 and data.shape[1] > 1:
            formats = data.shape[1] * formats
        cFormats = list(self.CommonStyleDefs.keys())
        for f in range(len(formats)):
            if formats[f] in cFormats:
                formats[f] = self.CommonStyleDefs[formats[f]]
        if len(titles) == 0:
            titles = ['col%s' % (s + 1) for s in range(data.shape[1])]

        cmd = '''splot '''

        for zi in z:
            if zi not in (x, y):
                cmd += '''"%s" u %s:%s:%s t "%s" %s, ''' % (
                    self.__DATA_FILE_PATH__,
                    x + 1,
                    y + 1,
                    zi + 1,
                    titles[zi],
                    formats[zi],
                )
        cmd = cmd[:-2]
        print(cmd)
        self.g_write(cmd)
        self.setKey(True)

    def save(self, name, directory=None, dfmt=None):
        """
        Save the last plot as a GnuPlot file *name*.plt which references
        *name*.dat.

        - *name* the name of the GnuPlot plt and and datafile
        - *directory* (optional) the directory to use (defaults to working directory)
        - *dfmt* is ignored and uses the value of self.DATF_FORMAT
        """
        if directory != None:
            out_n = os.path.join(directory, name)
        else:
            out_n = os.path.join(self.__WORK_DIR__, name)
        self.g_write('save "%s.plt"' % out_n)
        shutil.copy(
            os.path.join(self.__WORK_DIR__, self.__DATA_FILE_NAME__), "%s.dat" % out_n
        )
        F = open('%s.plt' % out_n, 'r')
        fnew = F.read().replace('_gnuplot.dat', '%s.dat' % name)
        fnew = fnew.replace(
            'noequal_axes', ''
        )  # fixes the "noequal_axes" bug in gnuplot save
        fnew = fnew.replace(
            '#!/gnuplot',
            '#!/gnuplot\n#\n# Plot created using PySCeS %s (http://pysces.sourceforge.net)'
            % __version__,
        )  # fixes the "noequal_axes" bug in gnuplot save
        F.close()
        F = open('%s.plt' % out_n, 'w')
        F.write(fnew)
        F.flush()
        F.close()

    def replot(self):
        """
        Replot the current GnuPlot plot
        """
        self.g_write('replot')
        self.g_pause()

    def replotAndWait(self, seconds=0.5):
        """
        Replot the current GnuPlot plot and wait default (*seconds* = 0.5)
        or until enter is pressed (*seconds* = -1)
        """
        self.replot()
        self.wait(seconds)

    def setLogScale(self, axis):
        """
        Set axis to logscale where:

        - *axis* = x, y, z, xy, xz, yz, zyx
        """
        axout = self.axisInputStringToList(axis)
        axstr = ''
        for s in axout:
            if s != None:
                axstr += s
        self.set('logscale', axstr)

    def setNoLogScale(self, axis):
        """
        Set axis to a linear scale where:

        - *axis* = x, y, z, xy, xz, yz, zyx
        """
        axout = self.axisInputStringToList(axis)
        axstr = ''
        for s in axout:
            if s != None:
                axstr += s
        self.unset('logscale', axstr)

    def setRange(self, axis, min=None, max=None):
        """
        Set axis range where:

        - *axis* = x, y, z, xy, xz, yz, zyx
        - *min* = range(s) lower bound (default=None) autoscale
        - *max* = range(s) upper bound (default=None) autoscale

        If only the *axis* argument is provided, GnuPlot will autoscale the
        ranges to the data.
        """
        axout = self.axisInputStringToList(axis)
        if min == None:
            min = '*'
        else:
            min = '%s' % min
        if max == None:
            max = '*'
        else:
            max = '%s' % max
        for ax in axout:
            if ax != None:
                self.set('%srange' % ax, '[ %s : %s ]' % (min, max))

    def setGrid(self, value):
        """
        Display or remove graph grid.

        - *value* (boolean) True (on) or False (off)
        """
        if value:
            self.set('grid')
        else:
            self.unset('grid')

    def setGraphTitle(self, title='PySCeS Plot'):
        """
        Set the graph title, unset if title argument is None

        - *title* (string, default='PySCeS Plot') the graph title
        """
        if title == None:
            self.unset('title')
        else:
            self.set('title', '\"%s\"' % title)

    def setKey(self, value=False):
        """
        Enable or disable the current plot key, no arguments removes key.

        - *value* boolean (default = False)
        """
        if value:
            self.set('key')
        else:
            self.unset('key')

    def setAxisLabel(self, axis, label=''):
        """
        Set the axis label:

        - *axis* = x, y, z, xy, xz, yz, zyx
        - *label* = string (default='')

        Called with only the axis argument clears the axis label.
        """
        axout = self.axisInputStringToList(axis)
        for ax in axout:
            if ax != None:
                self.set('%slabel \"%s\"' % (ax, label))

    def export(self, name, directory=None, type='png'):
        """
        Export the current plot as a <format> file.

        - *filename* the filename
        - *directory* optional (default = current working directory)
        - *type* the file format (default='png').

        Currently only PNG is guaranteed to be available in all interfaces.
        """
        imgPath = None
        if type.lower() == 'png':
            if name[-4:] != '.%s' % type:
                name += '.%s' % type
            if directory != None and os.path.exists(directory):
                imgPath = os.path.join(directory, name)
            else:
                imgPath = os.path.join(os.getcwd(), name)

        if imgPath != None and type in list(self.Terminals.keys()):
            self.set('output', '\"%s\"' % imgPath)
            self.setTerminal(type, self.Terminals[type])
            self.replot()
            self.set('output')
            self.setTerminal(
                self.__DISP_TERMINAL__, self.Terminals[self.__DISP_TERMINAL__]
            )
            print('Image exported as \"%s\"' % imgPath)
        else:
            print('Image not exported as type \"%s\"' % type)

    def setDataFileNumberFormat(self, format='%.8e'):
        """
        Sets the format string for data written to file

        - *format* format string (default='%.8e')
        """
        print('GnuPlot array data will be written as %s' % format % 0.123456789)
        self.DATF_FORMAT = format

    def setTerminal(self, name, options=''):
        """
        Sets the terminal, gnuplot: set terminal *name* *options*
        """
        self.set('terminal', '%s %s' % (name, options))

    def setMultiplot(self):
        """
        Begin a multiplot session
        """
        self.__MULTIPLOT__ = True
        self.set("multiplot")

    def unsetMultiplot(self):
        """
        End a multiplot session.
        """
        self.__MULTIPLOT__ = False
        self.unset("multiplot")

    def setSize(self, width=1.0, height=1.0):
        """
        Set the size of the next plot relative to the GnuPlot canvas (e.g. screen)
        size which is defined to be 1. For example if ``width = height = 0.5`` the
        plot is 1/4 the size of the viewable canvas. If no arguments are supplied
        reset size to 1,1.

        - *width* of next plot (default = 1.0)
        - *height* of next plot (default = 1.0)
        """
        self.set("size", "%s,%s" % (width, height))

    def setOrigin(self, xpos=0, ypos=0):
        """
        Set the origin (lower left corner) of the next plot. Uses GnuPlot
        screen coordinates. If no arguments are supplied reset origin to 0,0.

        - *xpos* of next plot (default = 0)
        - *ypos* of next plot (default = 0)
        """
        self.set("origin", "%s,%s" % (xpos, ypos))

    def setSizeAndOrigin(self, width=1, height=1, xpos=0, ypos=0):
        """
        Set the size and origin of the next plot. If no arguments
        are supplied, reset the size to 1,1 and origin to 0.0

        - *width* of next plot (default = 1.0)
        - *height* of next plot (default = 1.0)
        - *xpos* of next plot (default = 0)
        - *ypos* of next plot (default = 0)
        """
        self.setSize(width, height)
        self.setOrigin(xpos, ypos)


class MatplotlibUPI(PlotBase):
    """
    Refactored Matplotlib backend to the Unified Plotting Interface

    - *work_dir* (optional) working directory
    """

    __MODE_INTERACTIVE__ = True
    __MODE_NEW_FIGURE__ = True
    __FIGURE_GENERATOR__ = None
    __CURRENT_FIGURE__ = None
    __MODE_HOLD__ = False
    MAX_OPEN_WINDOWS = 10
    __ARRAY_DATA__ = None
    __WORK_DIR__ = None
    pyplot = None
    # __ENABLE_HTML5__ = False

    CommonStyleDefs = {'points': 'o', 'lines': '-'}
    __BACKENDS__ = (
        'GTK',
        'GTKAgg',
        'GTKCairo',
        'FltkAgg',
        'MacOSX',
        'QtAgg',
        'Qt4Agg',
        'TkAgg',
        'WX',
        'WXAgg',
        'CocoaAgg',
        'agg',
        'nbAgg',
        'cairo',
        'emf',
        'gdk',
        'pdf',
        'ps',
        'svg',
        'template',
    )

    __INTERACTIVE_BACKENDS__ = [
        'GTK',
        'GTKAgg',
        'GTKCairo',
        'FltkAgg',
        'MacOSX',
        'QtAgg',
        'Qt4Agg',
        'TkAgg',
        'WX',
        'WXAgg',
        'CocoaAgg',
        'nbAgg',
    ]

    __BACKEND__ = None

    def __init__(self, work_dir=None, backend=None):
        if work_dir != None and os.path.exists(work_dir):
            self.__WORK_DIR__ = work_dir
        else:
            self.__WORK_DIR__ = os.getcwd()
        try:
            import matplotlib

            if self.isnotebook():
                backend = 'nbAgg'
                import pysces

                pysces.__MATPLOTLIB_BACKEND__ = 'nbAgg'
            if backend in self.__INTERACTIVE_BACKENDS__:
                matplotlib.use(backend)
                self.__BACKEND__ = backend
                if not __SILENT_START__:
                    print(('Matplotlib backend set to: \"{}\"'.format(backend)))
            else:
                if backend == 'native':
                    print(
                        (
                            'Using natively configured Matplotlib backend: \"{}\"'.format(
                                matplotlib.get_backend()
                            )
                        )
                    )
                    self.__BACKEND__ = matplotlib.get_backend()
                else:
                    matplotlib.use('TkAgg')
                    self.__BACKEND__ = 'TkAgg'
                    print(
                        (
                            'Matplotlib \"{}\" backend not set, defaulting to: \"{}\"'.format(
                                backend, 'TkAgg'
                            )
                        )
                    )

            # if self.__ENABLE_HTML5__:
            # matplotlib.use('module://mplh5canvas.backend_h5canvas') # HTML5
            # else:
            # matplotlib.use('TKagg', warn=False)
        except Exception as ex:
            print(ex)
            print(
                "\nPySCeS defaults to matplotlib's TKagg backend if not specified \
                     in the user configuration file, set \"matplotlib_backend = <backend>\" "
            )

        from matplotlib import pyplot
        from matplotlib import pylab

        ##  self.pyplot = pyplot
        self.pyplot = pylab
        if self.__MODE_INTERACTIVE__:
            self.pyplot.ion()
        self._setNewFigureGenerator(self.MAX_OPEN_WINDOWS)

    def isnotebook(self):
        try:
            shell = get_ipython().__class__.__name__
            if shell == 'ZMQInteractiveShell':
                return True  # Jupyter notebook or qtconsole
            elif shell == 'TerminalInteractiveShell':
                return False  # Terminal running IPython
            else:
                return False  # Other type (?)
        except NameError:
            return False  # Probably standard Python interpreter

    def _setNewFigureGenerator(self, num):
        """
        Create a figure_generator with range [1,num+1]
        """
        self.__FIGURE_GENERATOR__ = itertools.cycle(list(range(1, num + 1)))
        self.MAX_OPEN_WINDOWS = num

    def _setActiveFigure(self, fignum):
        self.__CURRENT_FIGURE__ = self.pyplot.figure(fignum)

    def _setNewFigure(self):
        self._setActiveFigure(next(self.__FIGURE_GENERATOR__))

    def closeAll(self):
        """
        Close all open matplotlib figures.
        """
        for x in range(self.MAX_OPEN_WINDOWS):
            self.pyplot.close()

    def hold(self, hold=False):
        """
        Enable plot holding where each new graph is plotted on top of
        the previous one.

        - *hold* boolean (default = False)
        """
        if hold:
            self.__MODE_NEW_FIGURE__ = False
            self.__MODE_HOLD__ = True
        else:
            self.__MODE_NEW_FIGURE__ = True
            self.__MODE_HOLD__ = False

    def plot(self, data, x, y, title='', format='-'):
        """
        Plot a single line data[y] vs data[x] where:

        - *data* the data array
        - *x* x column index
        - *y* y column index
        - *title* is the line key
        - *format* is the Matplotlib format string (default='-')

        Format can also be the *CommonStyle* 'lines' or 'points'.
        """
        if format in list(self.CommonStyleDefs.keys()):
            format = self.CommonStyleDefs[format]
        if title == '':
            title = []
        else:
            title = data.shape[1] * [title]
        # print(title)

        self.plotLines(data, x, [y], title, [format] * data.shape[1])

    def plotLines(self, data, x, y=[], titles=[], formats=['-']):
        """
        Plot a multiple lines data[y1, y2, ] vs data[x] where:

        - *data* the data array
        - *x* x column index
        - *y* is a list of line indexes
        - *titles* is a list of line keys
        - *formats* is a list (per line) of Matplotlib format strings.

        If *formats* only contains a single item, this format is used for all lines
        and can also be the *CommonStyle* 'lines' or 'points'.
        """
        self.__ARRAY_DATA__ = data
        if self.__MODE_NEW_FIGURE__:
            self._setNewFigure()
        if not self.__MODE_HOLD__:
            self.pyplot.clf()
        if len(y) == 0:
            y = list(range(self.__ARRAY_DATA__.shape[1]))
            y.pop(x)

        if len(formats) != 1 and len(formats) != data.shape[1]:
            print(
                "len(titles) must be one or equal the number data columns (%s) using: \"%s\""
                % (data.shape[1], formats[-1])
            )
            formats = [copy.copy(formats[-1])]
        elif formats[0] == '' or formats[0] == None:
            formats = ['-']
        if len(formats) == 1 and data.shape[1] > 1:
            formats = data.shape[1] * formats
        cFormats = list(self.CommonStyleDefs.keys())
        for f in range(len(formats)):
            if formats[f] in cFormats:
                formats[f] = self.CommonStyleDefs[formats[f]]
        if len(titles) == 0:
            titles = ['col%s' % (l) for l in range(data.shape[1])]

        self.pyplot.ioff()
        for yi in y:
            # print(self.__ARRAY_DATA__.take([x],axis=1))
            # print(self.__ARRAY_DATA__.take([yi],axis=1))
            self.pyplot.plot(
                self.__ARRAY_DATA__.take([x], axis=1),
                self.__ARRAY_DATA__.take([yi], axis=1),
                formats[yi],
                label=titles[yi],
            )
        if self.__MODE_INTERACTIVE__:
            self.pyplot.ion()
        self.pyplot.legend()
        # if self.__ENABLE_HTML5__:
        # self.pyplot.show()

    def setLogScale(self, axis):
        """
        Set axis to logscale where:

        - *axis* = x, y, z, xy, xz, yz, zyx
        """
        axout = self.axisInputStringToList(axis)
        A = self.pyplot.gca()
        for ax in axout:
            if ax != None and ax in ['x', 'y']:
                getattr(A, 'set_%sscale' % ax)('log')
        if self.__MODE_INTERACTIVE__:
            self.pyplot.draw()

    def setNoLogScale(self, axis):
        """
        Set axis to a linear scale where:

        - *axis* = x, y, z, xy, xz, yz, zyx
        """
        axout = self.axisInputStringToList(axis)
        A = self.pyplot.gca()
        for ax in axout:
            if ax != None and ax in ['x', 'y']:
                getattr(A, 'set_%sscale' % ax)('linear')
        if self.__MODE_INTERACTIVE__:
            self.pyplot.draw()

    def setRange(self, axis, min=None, max=None):
        """
        Set axis range where

        - *axis* = x, y, z, xy, xz, yz, zyx
        - *min* = range(s) lower bound (default=None) autoscale
        - *max* = range(s) upper bound (default=None) autoscale
        """
        axout = self.axisInputStringToList(axis)
        for ax in axout:
            if ax != None and ax in ['x', 'y']:
                getattr(self.pyplot, '%slim' % ax)(min, max)
        ##  if self.__MODE_INTERACTIVE__: self.pyplot.draw()

    def setGrid(self, value):
        """
        Display or remove graph grid.

        - *value* (boolean) True (on) or False (off)
        """
        self.pyplot.grid(value)

    def setGraphTitle(self, title='PySCeS Plot'):
        """
        Set the graph title, unset if title=None

        - *title* (string, default='PySCeS Plot') the graph title
        """
        self.pyplot.title(title)

    def setKey(self, value=False):
        """
        Enable or disable the current plot key, no arguments removes key.

        - *value* boolean (default = False)
        """
        L = self.pyplot.gca().get_legend()
        if value:
            L._visible = True
        else:
            L._visible = False
        if self.__MODE_INTERACTIVE__:
            self.pyplot.draw()

    def setAxisLabel(self, axis, label=''):
        """
        Set the axis label:

        - *axis* = x, y, z, xy, xz, yz, zyx
        - *label* = string (default='')

        Called with only the axis argument clears the axis label.
        """
        axout = self.axisInputStringToList(axis)
        for ax in axout:
            if ax in ['x', 'y']:
                getattr(self.pyplot, '%slabel' % ax)(label)

    def save(self, name, directory=None, dfmt='%.8e'):
        """
        Save the plot data to

        - *filename* the filename
        - *directory* optional (default = current working directory)
        - *dfmt* the data format string (default='%.8e')
        """
        if directory != None:
            out_n = os.path.join(directory, name)
        else:
            out_n = os.path.join(self.__WORK_DIR__, name)
        numpy.savetxt(out_n, self.__ARRAY_DATA__, fmt=dfmt, delimiter=' ')
        print('Data saved as \"%s\"' % out_n)

    def export(self, name, directory=None, type='png'):
        """
        Export the current plot as a <format> file.

        - *filename* the filename
        - *directory* optional (default = current working directory)
        - *type* the file format (default='png').

        Currently only PNG is guaranteed to be available in all interfaces.
        """
        imgPath = None
        terminals = ['png']
        if type.lower() in terminals:
            if name[-4:] != '.%s' % type:
                name += '.%s' % type
            if directory != None and os.path.exists(directory):
                imgPath = os.path.join(directory, name)
            else:
                imgPath = os.path.join(os.getcwd(), name)
            self.pyplot.savefig(imgPath)
            print('Image exported as \"%s\"' % imgPath)
        else:
            print('Image not exported as type \"%s\"' % type)

    def setLineWidth(self, width=1):
        """
        Sets the line width for current axis

        - *width* the line width
        """
        [l.set_linewidth(width) for l in self.pyplot.gca().get_lines()]
        if self.__MODE_INTERACTIVE__:
            self.pyplot.draw()


class PyscesUPI(PlotBase):
    """
    This is the frontend to the PySCeS Unified Plotting Interface (pysces.plt.*) that
    allows one to specify which backend should be used to plot when a
    UPI method is called. More than one interface can be active at the same time and
    so far the Matplotlib and GnuPlot backends are available for use.

    This is an experiment which must be refactored into a more general
    way of doing things. Basically, I want an instance of the abstract plotting
    class which will plot to one, any or all currently available backends.
    If anybody has an idea how I can generate this class automatically please let
    me know ;-)
    """

    g = None
    m = None
    __USE_GNUPLOT__ = False
    __USE_MATPLOTLIB__ = False
    __GINTERFACES__ = ('matplotlib', 'gnuplot')

    def p_setInterface(self, name, instance):
        """
        Add an interface to the backend selector

        - *name* the interface name currently one of ['matplotlib','gnuplot']
        - *instance* an instance of a PlotBase derived (UPI) interface
        """
        if name in self.__GINTERFACES__:
            self.__setattr__('__USE_%s__' % name.upper(), True)
            self.__setattr__('%s' % name[0], instance)

    def p_deactivateInterface(self, interface):
        """
        Deactivate the interface. This does not delete the interface and
        it is possible to reactivate the deactivated interface with **p_activateInterface**.

        - *interface* one of ['matplotlib','gnuplot']
        """
        if interface in self.__GINTERFACES__:
            self.__setattr__('__USE_%s__' % interface.upper(), False)

    def p_activateInterface(self, interface):
        """
        Activate an interface that has been set with **p_setInterface()** but deactivated with
        **p_deactivateInterface**

        - *interface* one of ['matplotlib','gnuplot']
        """
        if (
            interface in self.__GINTERFACES__
            and self.__getattribute__('%s' % interface[0]) != None
        ):
            self.__setattr__('__USE_%s__' % interface.upper(), True)
        else:
            print('Cannot activate \"%s\"interface' % interface)

    def __gselect__(self, func, *args, **kwargs):
        """
        This methods uses plots a UPI method to *any* activated plotting
        interface
        """
        # avoids fishy nested args bug (not really relevant anymore but here just in case)
        if self.__USE_GNUPLOT__ and self.__USE_MATPLOTLIB__:
            # args2 = tuple([copy.deepcopy(a) for a in args])
            ##  args2 = copy.deepcopy(args)
            ##  kwargs2 = copy.deepcopy(kwargs)
            args2 = args
            kwargs2 = kwargs
        else:
            args2 = args
            kwargs2 = kwargs
        if self.__USE_MATPLOTLIB__:
            try:
                self.m.__getattribute__(func)(*args, **kwargs)
            except NotImplementedError:
                print('\"%s\" method not available using Matplotlib' % func)

        if self.__USE_GNUPLOT__:
            try:
                self.g.__getattribute__(func)(*args2, **kwargs2)
            except NotImplementedError:
                print('\"%s\" method not available using GnuPlot ' % func)
        del args, kwargs, args2, kwargs2

    def plot(self, data, x, y, title='', format=''):
        """
        Plot a single line data[y] vs data[x] where:

        - *data* the data array
        - *x* x column index
        - *y* y column index
        - *title* is the line key
        - *format* is the backend format string (default='')
        """
        self.__gselect__('plot', data, x, y, title, format)

    def plotLines(self, data, x, y=[], titles=[], formats=['']):
        """
        Plot a multiple lines data[y1, y2, ] vs data[x] where:

        - *data* the data array
        - *x* x column index
        - *y* is a list of line indexes, if empty all of y not including x is plotted
        - *titles* is a list of line keys, if empty Line1,Line2,Line3 is used
        - *formats* is a list (per line) of XXX format strings.

        If *formats* only contains a single item, this format is used for all lines.
        """
        self.__gselect__('plotLines', data, x, y, titles, formats)

    def splot(self, data, x, y, z, titles='', format=''):
        """
        Plot a surface data[z] vs data[y] vs data[x] where:

        - *data* the data array
        - *x* x column index
        - *y* y column index
        - *z* z column index
        - *titles* is a list of surface keys whose len matches data columns
        - *format* is the XXX format string (default='')
        """
        self.__gselect__('splot', data, x, y, z, titles, format)

    def splotSurfaces(self, data, x, y, z=[], titles=[], formats=['']):
        """
        Plot data[z1, z2, ] vs data[y] vs data[x] where:

        - *data* the data array
        - *x* x column index
        - *y* y column index
        - *z* list of z column indexes, if empty all of z not including x, y are plotted
        - *titles* is a list of surface keys, if empty Surf1, Surf2, Surf3 is used
        - *formats* is a list (per line) of XXX format strings (default='').

        If *formats* only contains a single item, this format is used for all surfaces.
        """
        self.__gselect__('splotSurfaces', data, x, y, z, titles, formats)

    def setLogScale(self, axis):
        """
        Set axis to logscale where:

        - *axis* = x, y, z, xy, xz, yz, zyx
        """
        self.__gselect__('setLogScale', axis)

    def setNoLogScale(self, axis):
        """
        Set axis to a linear scale where:

        - *axis* = x, y, z, xy, xz, yz, zyx
        """
        self.__gselect__('setNoLogScale', axis)

    def setRange(self, axis, min=None, max=None):
        """
        Set axis range where

        - *axis* = x, y, z, xy, xz, yz, zyx
        - *min* = range(s) lower bound (default=None) autoscale
        - *max* = range(s) upper bound (default=None) autoscale
        """
        self.__gselect__('setRange', axis, min, max)

    def setGrid(self, value):
        """
        Display or remove graph grid.

        - *value* (boolean) True (on) or False (off)
        """
        self.__gselect__('setGrid', value)

    def setGraphTitle(self, title='PySCeS Plot'):
        """
        Set the graph title, unset if title=None

        - *title* (string, default='PySCeS Plot') the graph title
        """
        self.__gselect__('setGraphTitle', title)

    def setKey(self, value=False):
        """
        Enable or disable the current plot key, no arguments removes key.

        - *value* boolean (default = False)
        """
        self.__gselect__('setKey', value)

    def setAxisLabel(self, axis, label=''):
        """
        Set the axis label:

        - *axis* = x, y, z, xy, xz, yz, zyx
        - *label* = string (default=None)

        Called with only the axis argument clears the axis label.
        """
        self.__gselect__('setAxisLabel', axis, label)

    def save(self, name, directory=None, dfmt='%.8e'):
        """
        Save the plot data and (optionally) XXX format file

        - *filename* the filename
        - *directory* optional (default = current working directory)
        - *dfmt* the data format string (default='%.8e')
        """
        self.__gselect__('save', name, directory, dfmt)

    def export(self, name, directory=None, type='png'):
        """
        Export the current plot as a <format> file.

        - *filename* the filename
        - *directory* optional (default = current working directory)
        - *type* the file format (default='png').

        Currently only PNG is guaranteed to be available in all interfaces.
        """
        self.__gselect__('export', name, directory, type)

    def replot(self):
        """
        Replot the current figure for all active interfaces
        """
        if self.__USE_GNUPLOT__:
            self.g.replot()
        if self.__USE_MATPLOTLIB__:
            self.m.pyplot.draw()

    def closeAll(self):
        """Close all active Matplolib figures"""
        if self.__USE_MATPLOTLIB__:
            self.m.closeAll()


class FIFOBuffer:
    """
    Simple fixed size FIFO buffer.
    """

    def __init__(self, size):
        self.data = [None for i in range(size)]

    def add(self, x):
        self.data.pop(0)
        self.data.append(x)

    def get(self):
        return self.data


# The following code will be adapted used in future for higher level
# graphics operations.
'''
class LineObj2(object):
    """
    Class describing a "graph line" or series of data.

    - **data** a ndarray column of data shape = (r, 1)
    - **label** a string label describing the column
    - **\*\*kwargs** (optional) key=value pairs of property descriptors

    Properties are stored as LineObj2 attributes
    """
    data = None
    label = None
    type = None
    axis = None
    IS_SELECTED = False

    __properties__ = ['index', 'linewidth', 'pointsize', 'colour', 'style',\
                    'logscale', 'visible']

    def __init__(self, data, label, **kwargs):
        self.axis = AxisObj2()
        self.setData(data)
        self.label = label
        self.type = data.dtype.char
        self.setProperties(**kwargs)

    def setProperties(self, **kwargs):
        for k in kwargs.keys():
            if k in self.__properties__:
                setattr(self, k, kwargs[k])

    def getData(self):
        return self.data

    def setData(self, data):
        self.axis.setRange(min(data), max(data))
        self.data = data

    def getLabel(self):
        return self.label

    def getType(self):
        return self.type

    def setSelected(self):
        self.IS_SELECTED = True

    def unsetSelected(self):
        self.IS_SELECTED = False

    def getProperties(self):
        out = {}
        for p in self.__properties__:
            if hasattr(self, p):
                out.update({p : getattr(self, p)})
        return out

    def unsetProperty(self, prop):
        if hasattr(self, prop):
            self.__delattr__(prop)

    def unsetAllProperties(self):
        for p in self.__properties__:
             self.unsetProperty(p)


class AxisObj2(object):
    """
    Class holding axis information (usually found as the *axis* attribute of
    a LineObj2):

    - **IS_AXIS** boolean
    - **IS_LOG** boolean
    - **axis** axis name 'x' or 'y'
    - **label** axis label such as 'Vdem'
    - **range_min** axis viewport minimum
    - **range_max** axis viewport maximum
    """
    IS_AXIS = False
    IS_LOG = False
    axis = None
    label = None
    units = None
    range_min = None
    range_max = None


    def setLog(self):
        """
        Set *IS_LOG = True* and convert range_(min/max) to log10
        """
        if not self.IS_LOG and self.range_min != None and self.range_max != None:
            self.setRange(numpy.log10(self.range_min), numpy.log10(self.range_max))
        self.IS_LOG = True

    def unsetLog(self):
        """
        Set *IS_LOG = False* and unlog range_(min/max)
        """
        if self.IS_LOG and self.range_min != None and self.range_max != None:
            self.setRange(10.0**self.range_min, 10.0**self.range_max)
        self.IS_LOG = False

    def setRange(self, min, max):
        """
        Set the axis range:

        - **min** lower bound
        - **max** upper bound
        """
        self.range_min = min
        self.range_max = max

    def setAxis(self, axis):
        """
        Enable this axis so that *IS_AXIS = True* where:

        - **axis** is either 'x' or 'y'

        If no *label* has been defined it is set to the axis name.
        """
        if axis in ['x','X']:
            self.IS_AXIS = True
            self.axis = 'x'
            if self.label == None or self.label in ['y','Y']:
                self.label = 'x'
        elif axis in ['y','Y']:
            self.IS_AXIS = True
            self.axis = 'y'
            if self.label == None or self.label in ['x','X']:
                self.label = 'y'
        else:
            print 'setAxis: %s is not a valid axis name' % axis

    def unsetAxis(self):
        """
        Disable this axis so that *IS_AXIS = False* and *axis = None*.
        """
        self.IS_AXIS = False
        self.axis = None


class GraphicsObj2(object):
    """
    Class describing a "graph" or collection of series.

    - **arr** an ndarray of data shape = (r, c>=1)
    - **labels** (optional) a list of column labels

    If labels are not supplied they are created as line1, line2 etc.
    Data series (LineObj2) objects are stored as GraphicsObj2
    instance attributes.
    """

    series = None
    scalar = None
    name = None

    def __init__(self, arr, labels=None):
        self.addSeriesFromArray(arr, labels)

    def setName(self, name):
        self.name = name

    def getName(self):
        return self.name

    def setAxis(self, axis, label):
        if hasattr(self, label):
            self.__getattribute__(label).axis.setAxis(axis)
        else:
            print 'setAxis: %s is not a valid label' % label

    def getAxis(self, axis, label=False):
        out = None
        for s in self.series:
            if s.axis.IS_AXIS and s.axis.axis == axis:
                if label:
                    out =  (s.getData(), s.getLabel())
                else:
                    out = s.getData()
                break
        return out

    def addOneSeries(self, series, label):
        L = LineObj2(series, label)
        self.series.append(L)
        self.__setattr__(label, L)

    def addSeriesFromArray(self, arr, labels=None):
        assert type(arr) == numpy.ndarray, '\nI need an ndarray'
        assert arr.shape[1] >= 1, '\nndarray must be 2D'
        self.series = []

        if labels == None or len(labels) != arr.shape[1]:
            labels = ['line%s' % c for c in range(arr.shape[1])]

        for c in range(arr.shape[1]):
            self.addOneSeries(arr.take([c], axis=1), labels[c])

    def hasSeries(self):
        return [l.getLabel() for l in self.series]

    def updateOneSeries(self, series, label):
        if hasattr(self, label):
            self.__getattribute__(label).setData(series)
        else:
            print '%s is an invalid series name.' % label

    def updateSeriesFromArray(self, arr):
        labels = self.hasSeries()
        if arr.shape[1] != len(labels):
            print 'updateSeriesFromArray: incorrect array shape for series (%s)' % self.hasSeries()
        else:
            for c in range(arr.shape[1]):
                self.updateOneSeries(arr.take([c], axis=1), labels[c])

    def delOneSeries(self, label):
        currS = self.hasSeries()
        if label in currS:
            self.series.pop(currS.index(label))
            self.__delattr__(label)

    def getAllSeriesAsArray(self):
        return numpy.hstack([s.getData() for s in self.series])

    def getAllSeriesAsDict(self):
        out = {}
        for s in self.series:
            out.update({s.getLabel() : s.getData()})
        return out

    def getAllSeriesAndProperties(self):
        out = {}
        for s in self.series:
            out.update({s.getLabel() : {'data' : s.getData(),
                                         'properties' : s.getProperties()}
                        })
        return out

    def setSeriesProperties(self, label, **kwargs):
        if hasattr(self, label):
            self.__getattribute__(label).setProperties(**kwargs)

    def setAllSeriesProperties(self, **kwargs):
        for s in self.hasSeries():
            self.setSeriesProperties(s, **kwargs)

    def unsetAllSeriesProperties(self):
        for s in self.hasSeries():
            self.__getattribute__(s).unsetAllProperties()
'''
