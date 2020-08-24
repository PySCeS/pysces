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

__doc__ = '''network and internet oriented utilities'''

from time import strftime
from getpass import getuser


class PyscesHTML:
    """PySCeS HTML formatting class: contains some basic html elements that can be used in generated reports."""

    __version__ = __version__

    def HTML_header(self, File):
        """
        HTML_header(File)

        Write an HTML page header to file (use with HTML_footer)

        Arguments:
        =========
        File: an open, writable Python file object

        """
        header = '\n'
        header += '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">\n'
        header += '<html>\n'
        header += '<head>\n'
        header += (
            '<title>PySCeS data generated at '
            + strftime("%H:%M:%S (%Z)")
            + '</title>\n'
        )
        header += (
            '<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1">\n'
        )
        header += '</head>\n'
        header += '<body>\n\n'
        header += '<h4><a href="http://pysces.sourceforge.net">PySCeS</a></h4>\n\n'
        File.write(header)
        File.write('<!-- PySCeS data generated at ' + strftime("%H:%M:%S") + '-->\n\n')
        return File

    def HTML_footer(self, File):
        """
        HTML_footer(File)

        Write an HTML page footer to file (use with HTML_header)

        Arguments:
        =========
        File: an open, writable Python file object

        """
        File.write(
            '\n<p><a href="http://pysces.sourceforge.net"><font size="3">PySCeS '
            + __version__
            + '</font></a><font size="2"> output\n generated at '
            + strftime("%H:%M:%S")
            + ' by <i>'
        )
        try:
            File.write(getuser())
        except:
            File.write('PySCeS')
        File.write('</i>)</font></p>\n')
        File.write('</body>\n')
        File.write('</html>\n\n')
        return File

    def par(self, str, File=None, align='l', txtout=0):
        """
        par(str,File=None,align='l',txtout=0)

        Format <par> text and write it to a file (or string)

        Arguments:
        =========
        str: the string of text to be written
        File [default=None]: an open, writable, Python file object
        align [default='l']: HTML alignment attribute ('l','c','r')
        txtout [default=0]: do not write to file (1) return formatted HTML string

        """
        if not txtout:
            assert type(File) == file, 'The 2nd argument needs to be an open file'
        if align == 'l':
            align = 'left'
        elif align == 'r':
            align == 'right'
        elif align == 'c':
            align = 'center'
        else:
            align = ''

        strout = '\n<p align="' + align + '">'
        cntr = 0
        max_str_len = 75
        seeker_active = 0
        for x in range(len(str)):
            cntr += 1
            strout += str[x]
            if seeker_active:
                if str[x] == ' ' or str[x] == '.' or str[x] == ',':
                    cntr = max_str_len
                    seeker_active = 0
            if cntr >= max_str_len:
                if str[x] == ' ' or str[x] == '.' or str[x] == ',':
                    strout += '\n '
                else:
                    seeker_active = 1
                cntr = 0

        strout += '\n</p>\n'
        if txtout:
            return strout
        else:
            File.write(strout)
        del str
        del strout

    def h1(self, str, File=None, align='l', txtout=0):
        """
        h1(str,File=None,align='l',txtout=0)

        Format <h1> text and write it to a file (or string)

        Arguments:
        =========
        str: the string of text to be written
        File [default=None]: an open, writable, Python file object
        align [default='l']: HTML alignment attribute ('l','c','r')
        txtout [default=0]: do not write to file (1) return formatted HTML string

        """
        if not txtout:
            assert type(File) == file, 'The 2nd argument needs to be an open file'
        if align == 'l':
            align = 'left'
        elif align == 'r':
            align == 'right'
        elif align == 'c':
            align = 'center'
        else:
            align = ''

        strout = '\n<h1 align="' + align + '">'

        cntr = 0
        max_str_len = 75
        seeker_active = 0
        for x in range(len(str)):
            cntr += 1
            strout += str[x]
            if seeker_active:
                if str[x] == ' ' or str[x] == '.' or str[x] == ',':
                    cntr = max_str_len
                    seeker_active = 0
            if cntr >= max_str_len:
                print(str[x])
                if str[x] == ' ' or str[x] == '.' or str[x] == ',':
                    strout += '\n '
                else:
                    seeker_active = 1
                cntr = 0

        strout += '\n</h1>\n'
        if txtout:
            return strout
        else:
            File.write(strout)
        del str
        del strout

    def h2(self, str, File=None, align='l', txtout=0):
        """
        h2(str,File=None,align='l',txtout=0)

        Format <h2> text and write it to a file (or string)

        Arguments:
        =========
        str: the string of text to be written
        File [default=None]: an open, writable, Python file object
        align [default='l']: HTML alignment attribute ('l','c','r')
        txtout [default=0]: do not write to file (1) return formatted HTML string

        """
        if not txtout:
            assert type(File) == file, 'The 2nd argument needs to be an open file'
        if align == 'l':
            align = 'left'
        elif align == 'r':
            align == 'right'
        elif align == 'c':
            align = 'center'
        else:
            align = ''

        strout = '\n<h2 align="' + align + '">'
        cntr = 0
        max_str_len = 75
        seeker_active = 0
        for x in range(len(str)):
            cntr += 1
            strout += str[x]
            if seeker_active:
                if str[x] == ' ' or str[x] == '.' or str[x] == ',':
                    cntr = max_str_len
                    seeker_active = 0
            if cntr >= max_str_len:
                print(str[x])
                if str[x] == ' ' or str[x] == '.' or str[x] == ',':
                    strout += '\n '
                else:
                    seeker_active = 1
                cntr = 0

        strout += '\n</h2>\n'
        if txtout:
            return strout
        else:
            File.write(strout)
        del str
        del strout

    def h3(self, str, File=None, align='l', txtout=0):
        """
        h3(str,File=None,align='l',txtout=0)

        Format <h3> text and write it to a file (or string)

        Arguments:
        =========
        str: the string of text to be written
        File [default=None]: an open, writable, Python file object
        align [default='l']: HTML alignment attribute ('l','c','r')
        txtout [default=0]: do not write to file (1) return formatted HTML string

        """
        if not txtout:
            assert type(File) == file, 'The 2nd argument needs to be an open file'
        if align == 'l':
            align = 'left'
        elif align == 'r':
            align == 'right'
        elif align == 'c':
            align = 'center'
        else:
            align = ''

        strout = '\n<h3 align="' + align + '">'
        cntr = 0
        max_str_len = 75
        seeker_active = 0
        for x in range(len(str)):
            cntr += 1
            strout += str[x]
            if seeker_active:
                if str[x] == ' ' or str[x] == '.' or str[x] == ',':
                    cntr = max_str_len
                    seeker_active = 0
            if cntr >= max_str_len:
                print(str[x])
                if str[x] == ' ' or str[x] == '.' or str[x] == ',':
                    strout += '\n'
                else:
                    seeker_active = 1
                cntr = 0

        strout += '\n</h3>\n'
        if txtout:
            return strout
        else:
            File.write(strout)
        del str
        del strout


import email
import email.utils
import mimetypes

import smtplib
from email.mime.text import MIMEText

from time import sleep, strftime
from getpass import getuser

import os


class PyscesSMTP:
    """A purely experimental class that extends PySCeS with SMTP mailer capabilities. Initialise with
    sender address and local mail server name."""

    __smtp_active = 0

    def __init__(self, fromadd, server):
        self.server = server
        try:
            self.userstr = getuser()
        except:
            self.userstr = 'PySCeS '
        self.msgintro = ''
        self.fromhead = self.userstr + ' <' + fromadd + '>'
        self.signature = (
            3 * '\n'
            + '---\nSent using PySCeS 0.2.2 (http://pysces.sourceforge.net/)\n '
        )

        # auto-open connection now closed
        # self.SMTPOpen()

    def GenericMail(self, toadd, msgtxt, subj='PySCeS generated email'):
        """
        GenericMail( toadd, msgtxt, subj='PySCeS generated email')

        Generate and send a text (non-mime) email message

        Arguments:
        =========
        toadd: recipient address
        msgtxt: the message body as a string
        subj [default='PySCeS generated email']: message subject line

        """
        assert type(msgtxt) == str, '\nMessage text must be a string'
        assert self.__smtp_active, 'SMTP Server not active\n'

        msgtxt = self.msgintro + msgtxt

        msgtxt += self.signature
        outer = MIMEText(msgtxt)

        outer['Subject'] = subj
        outer['To'] = toadd
        outer['From'] = self.fromhead
        outer['Date'] = email.Utils.formatdate(localtime='true')
        outer.epilogue = ' '

        if self.CheckGo():
            try:
                self.__SMTPserver.sendmail(self.fromhead, toadd, outer.as_string())
            except SMTPServerDisconnected as e:
                print(e)
                self.SMTPOpen()
                self.__SMTPserver.sendmail(self.fromhead, toadd, outer.as_string())
            sleep(0.2)
        else:
            print('\nEmail send aborted')

    def CheckGo(self):
        """
        CheckGo()

        Do you want to continue yes or no?
        Returns 1 or 0

        Arguments:
        None

        """
        GO = 1
        while GO:
            resp = input('\nDo you want to continue (yes/no): ')
            if resp.lower() == 'yes':
                print('OK.')
                GO = 0
                return 1
            elif resp.lower() == 'no':
                print('Skipped.')
                GO = 0
                return 0
            else:
                print('\nyes to continue, no to exit')

    ##	def GenericMailHTML(self, toadd, msgtxt, htmltxt, subj='PySCeS generated email'):
    ##        """
    ##        GenericMailHTML( toadd, msgtxt, htmltxt, subj='PySCeS generated email')
    ##
    ##        Generate a mime-compliant HTML email message
    ##
    ##        Arguments:
    ##        =========
    ##        toadd: recipient address
    ##        msgtxt: text only message string
    ##        htmltxt: html formatted message string
    ##        subj [default='PySCeS generated email']: the subject line
    ##
    ##        """
    ##		assert type(msgtxt) == str, '\nMessage text must be a string'
    ##		assert self.__smtp_active, 'SMTP Server not active\n'
    ##		# Create the enclosing (outer) message
    ##		outer = email.MIMEMultipart.MIMEMultipart()
    ##		outer['Subject'] = subj
    ##		outer['To'] = toadd
    ##		outer['From'] = self.fromhead
    ##		outer['Date'] = email.Utils.formatdate(localtime='true')
    ##		outer.preamble = ' \n'
    ##		outer.epilogue = '\n---\nGenerated by PySCeS 0.2.2\n '
    ##
    ##		msgtxt += self.signature
    ##		msg = email.MIMEText.MIMEText(msgtxt)
    ##		msg.add_header('Content-Disposition', 'inline')
    ##		outer.attach(msg)
    ##
    ##		self.__SMTPserver.sendmail(self.fromhead,toadd,outer.as_string())
    ##
    ##		ctype='text/plain'
    ##		maintype, subtype = ctype.split('/', 1)
    ##		fp = open(infile, 'r')
    ##		att = email.MIMEBase.MIMEBase(maintype, subtype)
    ##		att.set_payload(fp.read())
    ##		fp.close()
    ##		# Encode the payload using Base64
    ##		#email.Encoders.encode_base64(att)
    ##		# Set the filename parameter
    ##		att.add_header('Content-Disposition', 'attachment', filename=infile)
    ##		outer.attach(att)
    ##
    ##		SMTPserver.sendmail(fromhead,toadd,outer.as_string())
    ##
    ##		sleep(0.2)      #seconds

    def SMTPOpen(self):
        """
        SMTPOpen()

        Start client and connect to an SMTP server

        Arguments:
        None

        """
        self.__SMTPserver = smtplib.SMTP(self.server)
        self.__smtp_active = 1
        print('\nSMTP server connection opened\n')

    def SMTPClose(self):
        """
        SMTPClose()

        Close connection to SMTP server

        Arguments:
        None

        """
        self.__SMTPserver.close()
        self.__smtp_active = 0
        print('\nSMTP server connection closed\n')


if __name__ == '__main__':
    replyTo = 'bgoli@sun.ac.za'
    server = 'mail.sun.ac.za'
    print('Reply to:', replyTo)
    print('SMTP server:', server)
    smtp = PyscesSMTP(replyTo, server)
    smtp.GenericMail(
        'bgoli@sun.ac.za',
        'This test message created: ' + strftime("%a, %d %b %Y %H:%M:%S"),
    )
    # smtp.GenericMail('jr@sun.ac.za','This test message created: '+ strftime("%a, %d %b %Y %H:%M:%S"))
    smtp.SMTPClose()
