######################################################################################
#  Author: Junjun Guo
#  E-mail: jjguo2@bjtu.edu.cn/guojj_ce@163.com
#    Date: 1/16/2021
#  Environemet: Successfully excucted in python 3.11
######################################################################################
from sectionFiberMain import circleSection
name="circleHole" #section name
outD = 2  # the diameter of the outside circle
coverThick = 0.06  # the thinckness of the cover concrete
outbarD = 0.03  # outside bar diameter
outbarDist = 0.15  # outside bar space
coreSize = 0.1  # the size of core concrete fiber
coverSize = 0.1  # the size of cover concrete fiber
plotState = True  # plot the fiber or not plot=True or False
inD =1 # the diameter of the inside circle
inBarD=0.03 # inside bar diameter
inBarDist=0.15 # inside bar space
corFiber, coverFiber, barFiber = circleSection(name,outD, coverThick, outbarD, outbarDist, coreSize, coverSize,
                                                   plotState,inD,inBarD,inBarDist)