# Copyright 2019 Pascal Audet & Andrew Schaeffer
#
# This file is part of SplitPy.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""

Module containing the main GUI widget classes based on :mod:`~PyQt5`:

- :class:`~splitpy.gui.Pick`
- :class:`~splitpy.gui.Keep`
- :class:`~splitpy.gui.Save`
- :class:`~splitpy.gui.Repeat`

These classes simply open a Qt5 message box and contain a boolean attribute
corresponding the anwer to simple yes/no buttons. 

"""

# -*- coding: utf-8 -*-
from PyQt5.QtWidgets import QWidget, QPushButton, QMessageBox
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import pyqtSlot


class Pick(QWidget):
    """
    GUI message box for picking new SKS analysis window

    """

    def __init__(self):
        super().__init__()
        self.title = 'PyQt5 messagebox - pythonspot.com'
        self.left = 10
        self.top = 10
        self.width = 520
        self.height = 200
        self.initUI()
        
    def initUI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
        buttonReply = QMessageBox.question(self, 'Pick window', "Re/Pick window? (Click no to continue)", QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        self.reply = buttonReply == QMessageBox.Yes
        self.show()
        self.close()

class Keep(QWidget):
    """
    GUI message box for keeping current split estimates

    """

    def __init__(self):
        super().__init__()
        self.title = 'PyQt5 messagebox - pythonspot.com'
        self.left = 10
        self.top = 10
        self.width = 520
        self.height = 200
        self.initUI()
        
    def initUI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
        buttonReply = QMessageBox.question(self, 'Keep estimates', "Keep estimates?", QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        self.reply = buttonReply == QMessageBox.Yes
        self.show()
        self.close()

class Save(QWidget):
    """
    GUI message box for overwriting saved split

    """

    def __init__(self):
        super().__init__()
        self.title = 'PyQt5 messagebox - pythonspot.com'
        self.left = 10
        self.top = 10
        self.width = 520
        self.height = 200
        self.initUI()
        
    def initUI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
        buttonReply = QMessageBox.question(self, 'Save Split Result', "Split Result Exists. Overwrite?", QMessageBox.Yes | QMessageBox.No, QMessageBox.Yes)
        self.reply = buttonReply == QMessageBox.Yes
        self.show()
        self.close()

class Repeat(QWidget):
    """
    GUI message box for repeating a saved result

    """

    def __init__(self):
        super().__init__()
        self.title = 'PyQt5 messagebox - pythonspot.com'
        self.left = 10
        self.top = 10
        self.width = 520
        self.height = 200
        self.initUI()
        
    def initUI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
        buttonReply = QMessageBox.question(self, 'Repeat Saved Split', "Split Results Exist. Repeat?", QMessageBox.Yes | QMessageBox.No, QMessageBox.Yes)
        self.reply = buttonReply == QMessageBox.Yes
        self.show()
        self.close()
