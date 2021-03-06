# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'c:\Users\manuel.collongues\AppData\Roaming\QGIS\QGIS3\profiles\default\python\plugins\interpolation_shpm\interpolation_shpm_dialog_base.ui'
#
# Created by: PyQt5 UI code generator 5.15.4
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_SHPMDialogBase(object):
    def setupUi(self, SHPMDialogBase):
        SHPMDialogBase.setObjectName("SHPMDialogBase")
        SHPMDialogBase.resize(696, 958)
        self.button_box = QtWidgets.QDialogButtonBox(SHPMDialogBase)
        self.button_box.setGeometry(QtCore.QRect(600, 920, 81, 32))
        self.button_box.setOrientation(QtCore.Qt.Horizontal)
        self.button_box.setStandardButtons(QtWidgets.QDialogButtonBox.Close)
        self.button_box.setObjectName("button_box")
        self.pushshp_rivg = QtWidgets.QPushButton(SHPMDialogBase)
        self.pushshp_rivg.setGeometry(QtCore.QRect(570, 20, 41, 23))
        self.pushshp_rivg.setObjectName("pushshp_rivg")
        self.lineshp_rg = QtWidgets.QLineEdit(SHPMDialogBase)
        self.lineshp_rg.setGeometry(QtCore.QRect(230, 20, 321, 20))
        self.lineshp_rg.setObjectName("lineshp_rg")
        self.label_select_rg = QtWidgets.QLabel(SHPMDialogBase)
        self.label_select_rg.setGeometry(QtCore.QRect(60, 20, 161, 16))
        self.label_select_rg.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label_select_rg.setObjectName("label_select_rg")
        self.label_select_rd = QtWidgets.QLabel(SHPMDialogBase)
        self.label_select_rd.setGeometry(QtCore.QRect(70, 60, 151, 16))
        self.label_select_rd.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label_select_rd.setObjectName("label_select_rd")
        self.pushshp_rivd = QtWidgets.QPushButton(SHPMDialogBase)
        self.pushshp_rivd.setGeometry(QtCore.QRect(570, 60, 41, 23))
        self.pushshp_rivd.setObjectName("pushshp_rivd")
        self.lineshp_rd = QtWidgets.QLineEdit(SHPMDialogBase)
        self.lineshp_rd.setGeometry(QtCore.QRect(230, 60, 321, 20))
        self.lineshp_rd.setObjectName("lineshp_rd")
        self.label_select_lidar = QtWidgets.QLabel(SHPMDialogBase)
        self.label_select_lidar.setGeometry(QtCore.QRect(70, 100, 151, 16))
        self.label_select_lidar.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label_select_lidar.setObjectName("label_select_lidar")
        self.line_raster_lidar = QtWidgets.QLineEdit(SHPMDialogBase)
        self.line_raster_lidar.setGeometry(QtCore.QRect(230, 100, 321, 20))
        self.line_raster_lidar.setObjectName("line_raster_lidar")
        self.push_raster_lidar = QtWidgets.QPushButton(SHPMDialogBase)
        self.push_raster_lidar.setGeometry(QtCore.QRect(570, 100, 41, 23))
        self.push_raster_lidar.setObjectName("push_raster_lidar")
        self.line_maillage = QtWidgets.QLineEdit(SHPMDialogBase)
        self.line_maillage.setGeometry(QtCore.QRect(230, 370, 321, 20))
        self.line_maillage.setObjectName("line_maillage")
        self.label_valeur_maillage = QtWidgets.QLabel(SHPMDialogBase)
        self.label_valeur_maillage.setGeometry(QtCore.QRect(10, 370, 211, 20))
        self.label_valeur_maillage.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label_valeur_maillage.setObjectName("label_valeur_maillage")
        self.line_largeur_vallee = QtWidgets.QLineEdit(SHPMDialogBase)
        self.line_largeur_vallee.setGeometry(QtCore.QRect(230, 410, 321, 20))
        self.line_largeur_vallee.setObjectName("line_largeur_vallee")
        self.label_valeur_largeur_vallee = QtWidgets.QLabel(SHPMDialogBase)
        self.label_valeur_largeur_vallee.setGeometry(QtCore.QRect(10, 410, 211, 20))
        self.label_valeur_largeur_vallee.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label_valeur_largeur_vallee.setObjectName("label_valeur_largeur_vallee")
        self.label_nom_couche_profils = QtWidgets.QLabel(SHPMDialogBase)
        self.label_nom_couche_profils.setGeometry(QtCore.QRect(10, 450, 211, 20))
        self.label_nom_couche_profils.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label_nom_couche_profils.setObjectName("label_nom_couche_profils")
        self.line_couche_profils = QtWidgets.QLineEdit(SHPMDialogBase)
        self.line_couche_profils.setGeometry(QtCore.QRect(230, 450, 321, 20))
        self.line_couche_profils.setObjectName("line_couche_profils")
        self.label_nom_couche_branche = QtWidgets.QLabel(SHPMDialogBase)
        self.label_nom_couche_branche.setGeometry(QtCore.QRect(10, 490, 211, 20))
        self.label_nom_couche_branche.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label_nom_couche_branche.setObjectName("label_nom_couche_branche")
        self.line_couche_branche = QtWidgets.QLineEdit(SHPMDialogBase)
        self.line_couche_branche.setGeometry(QtCore.QRect(230, 490, 321, 20))
        self.line_couche_branche.setObjectName("line_couche_branche")
        self.push_calcul_interpolation = QtWidgets.QPushButton(SHPMDialogBase)
        self.push_calcul_interpolation.setGeometry(QtCore.QRect(300, 530, 181, 23))
        self.push_calcul_interpolation.setObjectName("push_calcul_interpolation")
        self.text_log = QtWidgets.QTextBrowser(SHPMDialogBase)
        self.text_log.setGeometry(QtCore.QRect(20, 590, 651, 321))
        self.text_log.setLocale(QtCore.QLocale(QtCore.QLocale.French, QtCore.QLocale.France))
        self.text_log.setObjectName("text_log")
        self.label_nom_couche_lidar = QtWidgets.QLabel(SHPMDialogBase)
        self.label_nom_couche_lidar.setGeometry(QtCore.QRect(10, 330, 211, 20))
        self.label_nom_couche_lidar.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label_nom_couche_lidar.setObjectName("label_nom_couche_lidar")
        self.line_couche_lidar = QtWidgets.QLineEdit(SHPMDialogBase)
        self.line_couche_lidar.setGeometry(QtCore.QRect(230, 330, 321, 20))
        self.line_couche_lidar.setObjectName("line_couche_lidar")
        self.line_couche_rg = QtWidgets.QLineEdit(SHPMDialogBase)
        self.line_couche_rg.setGeometry(QtCore.QRect(230, 250, 321, 20))
        self.line_couche_rg.setObjectName("line_couche_rg")
        self.label_nom_couche_rd = QtWidgets.QLabel(SHPMDialogBase)
        self.label_nom_couche_rd.setGeometry(QtCore.QRect(10, 290, 211, 20))
        self.label_nom_couche_rd.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label_nom_couche_rd.setObjectName("label_nom_couche_rd")
        self.line_couche_rd = QtWidgets.QLineEdit(SHPMDialogBase)
        self.line_couche_rd.setGeometry(QtCore.QRect(230, 290, 321, 20))
        self.line_couche_rd.setObjectName("line_couche_rd")
        self.label_nom_couche_rg = QtWidgets.QLabel(SHPMDialogBase)
        self.label_nom_couche_rg.setGeometry(QtCore.QRect(10, 250, 211, 20))
        self.label_nom_couche_rg.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label_nom_couche_rg.setObjectName("label_nom_couche_rg")
        self.label = QtWidgets.QLabel(SHPMDialogBase)
        self.label.setGeometry(QtCore.QRect(310, 200, 111, 16))
        self.label.setObjectName("label")
        self.label_select_semis_interp = QtWidgets.QLabel(SHPMDialogBase)
        self.label_select_semis_interp.setGeometry(QtCore.QRect(60, 140, 161, 20))
        self.label_select_semis_interp.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label_select_semis_interp.setObjectName("label_select_semis_interp")
        self.lineshp_semi_interp = QtWidgets.QLineEdit(SHPMDialogBase)
        self.lineshp_semi_interp.setGeometry(QtCore.QRect(230, 140, 321, 20))
        self.lineshp_semi_interp.setObjectName("lineshp_semi_interp")
        self.pushshp_semis_interp = QtWidgets.QPushButton(SHPMDialogBase)
        self.pushshp_semis_interp.setGeometry(QtCore.QRect(570, 140, 41, 23))
        self.pushshp_semis_interp.setObjectName("pushshp_semis_interp")

        self.retranslateUi(SHPMDialogBase)
        self.button_box.accepted.connect(SHPMDialogBase.accept)
        self.button_box.rejected.connect(SHPMDialogBase.reject)
        QtCore.QMetaObject.connectSlotsByName(SHPMDialogBase)

    def retranslateUi(self, SHPMDialogBase):
        _translate = QtCore.QCoreApplication.translate
        SHPMDialogBase.setWindowTitle(_translate("SHPMDialogBase", "Interpolation SHPM"))
        self.pushshp_rivg.setText(_translate("SHPMDialogBase", "..."))
        self.label_select_rg.setText(_translate("SHPMDialogBase", "Charger shape rive gauche"))
        self.label_select_rd.setText(_translate("SHPMDialogBase", "Charger shape rive droite"))
        self.pushshp_rivd.setText(_translate("SHPMDialogBase", "..."))
        self.label_select_lidar.setText(_translate("SHPMDialogBase", "Charger LIDAR"))
        self.push_raster_lidar.setText(_translate("SHPMDialogBase", "..."))
        self.line_maillage.setText(_translate("SHPMDialogBase", "10"))
        self.label_valeur_maillage.setText(_translate("SHPMDialogBase", "Maillage (-1 r??cup??re valeur de la branche)"))
        self.line_largeur_vallee.setText(_translate("SHPMDialogBase", "500"))
        self.label_valeur_largeur_vallee.setText(_translate("SHPMDialogBase", "Largeur de la vall??e"))
        self.label_nom_couche_profils.setText(_translate("SHPMDialogBase", "Nom de la couche des profils"))
        self.line_couche_profils.setText(_translate("SHPMDialogBase", "profiles"))
        self.label_nom_couche_branche.setText(_translate("SHPMDialogBase", "Nom de la couche des branches"))
        self.line_couche_branche.setText(_translate("SHPMDialogBase", "branchs"))
        self.push_calcul_interpolation.setText(_translate("SHPMDialogBase", "Calcul interpolation"))
        self.text_log.setPlaceholderText(_translate("SHPMDialogBase", "hjkljkj"))
        self.label_nom_couche_lidar.setText(_translate("SHPMDialogBase", "Couche LIDAR"))
        self.line_couche_lidar.setText(_translate("SHPMDialogBase", "interpolation_lidar"))
        self.line_couche_rg.setText(_translate("SHPMDialogBase", "interpolation_rive_gauche"))
        self.label_nom_couche_rd.setText(_translate("SHPMDialogBase", "Couche rive droite"))
        self.line_couche_rd.setText(_translate("SHPMDialogBase", "interpolation_rive_droite"))
        self.label_nom_couche_rg.setText(_translate("SHPMDialogBase", "Couche rive gauche"))
        self.label.setText(_translate("SHPMDialogBase", "Param??tres"))
        self.label_select_semis_interp.setText(_translate("SHPMDialogBase", "Semis de point interpol?? [sortie]"))
        self.pushshp_semis_interp.setText(_translate("SHPMDialogBase", "..."))
