#-*-coding: UTF-8-*-
#########################################################################
#  Author: Junjun Guo
#  E-mail: jjguo2@bjtu.edu.cn/guojj_ce@163.com
#  Environemet: Successfully excucted in python 3.11
#########################################################################
#import necessary modules
import ctypes
ll=ctypes.cdll.LoadLibrary
try:
	import os
	currentPath=os.getcwd()
	os.add_dll_directory(currentPath+"\\auxiliaryModules\\CalculateIMs")
	lib_responseSpect=ll("responseSpect.dll")
except:
	lib_responseSpect=ll("responseSpect.dll")
lib_responseSpect.responseSpect.restype=None
lib_responseSpect.responseSpect.argtypes=(ctypes.c_double,ctypes.c_double,ctypes.POINTER(ctypes.c_double),
	ctypes.POINTER(ctypes.c_double),ctypes.c_int,ctypes.c_int,ctypes.c_double,ctypes.POINTER(ctypes.c_double),
	ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double))
def returnSaSvSd(beta,dt,Tarray,accArray,accLength,Tsizes,m,SaArray,SvArray,SdArray):
	"""
	Warper function for calculate acceleration,velocity and displacement response spectra based on C++ dll
	:param beta: damping ratio
	:param dt:time interval(s)
	:param Tarray: periods list
	:param accArray: acceleration time history (g)
	:param sizes: number of response spectra points
	:param m: mass,default=1
	:param SaArray: store acc spectra and initial with zero list
	:param SvArray:store vel spectra and initial with zero list
	:param SdArray:store disp spectra and initial with zero list
	:return:Sa-acceleration spectra list(g)
			Sv-velocity spectra list(cm/s)
			Sd-displacement spectra list(cm)
	"""
	Tarray_size=len(Tarray)
	Tarray=(ctypes.c_double*Tarray_size)(*Tarray)
	accArray_size=len(accArray)
	accArray=(ctypes.c_double*accArray_size)(*accArray)
	SaArray_size=len(SaArray)
	SaArray=(ctypes.c_double*SaArray_size)(*SaArray)
	SvArray_size=len(SvArray)
	SvArray=(ctypes.c_double*SvArray_size)(*SvArray)
	SdArray_size=len(SdArray)
	SdArray=(ctypes.c_double*SdArray_size)(*SdArray)
	lib_responseSpect.responseSpect(beta,dt,Tarray,accArray,accLength,Tsizes,m,SaArray,SvArray,SdArray)
	return SaArray,SvArray,SdArray
