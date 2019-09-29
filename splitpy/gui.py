'''
SUBMODULE gui.py

Module containing the main GUI widget classes

'''

from PyQt5.QtWidgets import QWidget, QPushButton, QMessageBox
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import pyqtSlot


class Pick(QWidget):
    """pick: GUI for picking"""

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
    """keep: GUI for saving estimates"""

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
    """save: GUI for overwriting saved split"""

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
    """save: GUI for repeating a saved result"""

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
