# -*- coding: utf-8 -*-
from __future__ import division
#Declarations
#The dictionary of parameters
#name,bname,type,family,unit,value,mode,description,group,min,max,list,enable,iscombocheckbox,isused
parameterDict = {}
try:
	if Parameter:
		pass
except NameError:
	class Parameter:
		def __init__(self, **d):
			pass
#Type:String
#Mode:In
#Description:ARC/Venus tool  type
#List:ARC475/ARC675/Venus475/Venus675
tool_type = u"ARC475"
parameterDict.update({'tool_type' : Parameter(name='tool_type',bname='',type='String',family='',unit='',value='ARC475',mode='In',description='ARC/Venus tool  type',group='',min='',max='',list='ARC475/ARC675/Venus475/Venus675',enable='True',iscombocheckbox='False',isused='True')})
#Type:String
#Mode:In
#Description:Mud Type
#List:WBM/OBM
mud_type = u"OBM"
parameterDict.update({'mud_type' : Parameter(name='mud_type',bname='',type='String',family='',unit='',value='OBM',mode='In',description='Mud Type',group='',min='',max='',list='WBM/OBM',enable='True',iscombocheckbox='False',isused='True')})
#Type:Number
#Unit:in
bit_size_unit = u"in"
#Mode:In
#Description:Bit Size
#Minimum:
#Maximum:
#List:
bit_size = 8.5
parameterDict.update({'bit_size' : Parameter(name='bit_size',bname='',type='Number',family='',unit='in',value='8.5',mode='In',description='Bit Size',group='',min='',max='',list='',enable='True',iscombocheckbox='False',isused='True')})
#Type:Variable
#BName:dept
#Family:Measured Depth
#Unit:FT
#Mode:In
#Description:Measured depth
DEPT = Variable("Zabazaba - 1X", "ARC5_MWD_010PDP_60B", "TDEP", u"Measured Depth", u"FT")
parameterDict.update({'DEPT' : Parameter(name='DEPT',bname='dept',type='Variable',family='Measured Depth',unit='FT',value='Zabazaba - 1X.ARC5_MWD_010PDP_60B.TDEP',mode='In',description='Measured depth',group='',min='',max='',list='',enable='True',iscombocheckbox='False',isused='True')})
#Type:Variable
#BName:a10h_unc
#Family:Resistivity - Attenuation
#Unit:ohm.m
#Mode:In
#Description:Required ARC channel
A10H_UNC = Variable("Zabazaba - 1X", "ARC5_MWD_010PDP_60B", "A10H_UNC", u"Resistivity - Attenuation", u"ohm.m")
parameterDict.update({'A10H_UNC' : Parameter(name='A10H_UNC',bname='a10h_unc',type='Variable',family='Resistivity - Attenuation',unit='ohm.m',value='Zabazaba - 1X.ARC5_MWD_010PDP_60B.A10H_UNC',mode='In',description='Required ARC channel',group='',min='',max='',list='',enable='True',iscombocheckbox='False',isused='True')})
#Type:Variable
#BName:a16h_unc
#Family:Resistivity - Attenuation
#Unit:ohm.m
#Mode:In
#Description:Required ARC channel
A16H_UNC = Variable("Zabazaba - 1X", "ARC5_MWD_010PDP_60B", "A16H_UNC", u"Resistivity - Attenuation", u"ohm.m")
parameterDict.update({'A16H_UNC' : Parameter(name='A16H_UNC',bname='a16h_unc',type='Variable',family='Resistivity - Attenuation',unit='ohm.m',value='Zabazaba - 1X.ARC5_MWD_010PDP_60B.A16H_UNC',mode='In',description='Required ARC channel',group='',min='',max='',list='',enable='True',iscombocheckbox='False',isused='True')})
#Type:Variable
#BName:a22h_unc
#Family:Resistivity - Attenuation
#Unit:ohm.m
#Mode:In
#Description:Required ARC channel
A22H_UNC = Variable("Zabazaba - 1X", "ARC5_MWD_010PDP_60B", "A22H_UNC", u"Resistivity - Attenuation", u"ohm.m")
parameterDict.update({'A22H_UNC' : Parameter(name='A22H_UNC',bname='a22h_unc',type='Variable',family='Resistivity - Attenuation',unit='ohm.m',value='Zabazaba - 1X.ARC5_MWD_010PDP_60B.A22H_UNC',mode='In',description='Required ARC channel',group='',min='',max='',list='',enable='True',iscombocheckbox='False',isused='True')})
#Type:Variable
#BName:a28h_unc
#Family:Resistivity - Attenuation
#Unit:ohm.m
#Mode:In
#Description:Required ARC channel
A28H_UNC = Variable("Zabazaba - 1X", "ARC5_MWD_010PDP_60B", "A28H_UNC", u"Resistivity - Attenuation", u"ohm.m")
parameterDict.update({'A28H_UNC' : Parameter(name='A28H_UNC',bname='a28h_unc',type='Variable',family='Resistivity - Attenuation',unit='ohm.m',value='Zabazaba - 1X.ARC5_MWD_010PDP_60B.A28H_UNC',mode='In',description='Required ARC channel',group='',min='',max='',list='',enable='True',iscombocheckbox='False',isused='True')})
#Type:Variable
#BName:a34h_unc
#Family:Resistivity - Attenuation
#Unit:ohm.m
#Mode:In
#Description:Required ARC channel
A34H_UNC = Variable("Zabazaba - 1X", "ARC5_MWD_010PDP_60B", "A34H_UNC", u"Resistivity - Attenuation", u"ohm.m")
parameterDict.update({'A34H_UNC' : Parameter(name='A34H_UNC',bname='a34h_unc',type='Variable',family='Resistivity - Attenuation',unit='ohm.m',value='Zabazaba - 1X.ARC5_MWD_010PDP_60B.A34H_UNC',mode='In',description='Required ARC channel',group='',min='',max='',list='',enable='True',iscombocheckbox='False',isused='True')})
#Type:Variable
#BName:p10h_unc
#Family:Resistivity - Phase Shift
#Unit:ohm.m
#Mode:In
#Description:Required ARC channel
P10H_UNC = Variable("Zabazaba - 1X", "ARC5_MWD_010PDP_60B", "P10H_UNC", u"Resistivity - Phase Shift", u"ohm.m")
parameterDict.update({'P10H_UNC' : Parameter(name='P10H_UNC',bname='p10h_unc',type='Variable',family='Resistivity - Phase Shift',unit='ohm.m',value='Zabazaba - 1X.ARC5_MWD_010PDP_60B.P10H_UNC',mode='In',description='Required ARC channel',group='',min='',max='',list='',enable='True',iscombocheckbox='False',isused='True')})
#Type:Variable
#BName:p16h_unc
#Family:Resistivity - Phase Shift
#Unit:ohm.m
#Mode:In
#Description:Required ARC channel
P16H_UNC = Variable("Zabazaba - 1X", "ARC5_MWD_010PDP_60B", "P16H_UNC", u"Resistivity - Phase Shift", u"ohm.m")
parameterDict.update({'P16H_UNC' : Parameter(name='P16H_UNC',bname='p16h_unc',type='Variable',family='Resistivity - Phase Shift',unit='ohm.m',value='Zabazaba - 1X.ARC5_MWD_010PDP_60B.P16H_UNC',mode='In',description='Required ARC channel',group='',min='',max='',list='',enable='True',iscombocheckbox='False',isused='True')})
#Type:Variable
#BName:p22h_unc
#Family:Resistivity - Phase Shift
#Unit:ohm.m
#Mode:In
#Description:Required ARC channel
P22H_UNC = Variable("Zabazaba - 1X", "ARC5_MWD_010PDP_60B", "P22H_UNC", u"Resistivity - Phase Shift", u"ohm.m")
parameterDict.update({'P22H_UNC' : Parameter(name='P22H_UNC',bname='p22h_unc',type='Variable',family='Resistivity - Phase Shift',unit='ohm.m',value='Zabazaba - 1X.ARC5_MWD_010PDP_60B.P22H_UNC',mode='In',description='Required ARC channel',group='',min='',max='',list='',enable='True',iscombocheckbox='False',isused='True')})
#Type:Variable
#BName:p28h_unc
#Family:Resistivity - Phase Shift
#Unit:ohm.m
#Mode:In
#Description:Required ARC channel
P28H_UNC = Variable("Zabazaba - 1X", "ARC5_MWD_010PDP_60B", "P28H_UNC", u"Resistivity - Phase Shift", u"ohm.m")
parameterDict.update({'P28H_UNC' : Parameter(name='P28H_UNC',bname='p28h_unc',type='Variable',family='Resistivity - Phase Shift',unit='ohm.m',value='Zabazaba - 1X.ARC5_MWD_010PDP_60B.P28H_UNC',mode='In',description='Required ARC channel',group='',min='',max='',list='',enable='True',iscombocheckbox='False',isused='True')})
#Type:Variable
#BName:p34h_unc
#Family:Resistivity - Phase Shift
#Unit:ohm.m
#Mode:In
#Description:Required ARC channel
P34H_UNC = Variable("Zabazaba - 1X", "ARC5_MWD_010PDP_60B", "P34H_UNC", u"Resistivity - Phase Shift", u"ohm.m")
parameterDict.update({'P34H_UNC' : Parameter(name='P34H_UNC',bname='p34h_unc',type='Variable',family='Resistivity - Phase Shift',unit='ohm.m',value='Zabazaba - 1X.ARC5_MWD_010PDP_60B.P34H_UNC',mode='In',description='Required ARC channel',group='',min='',max='',list='',enable='True',iscombocheckbox='False',isused='True')})
#Type:Variable
#BName:a10l_unc
#Family:Resistivity - Attenuation
#Unit:ohm.m
#Mode:In
#Description:Required ARC channel
A10L_UNC = Variable("Zabazaba - 1X", "ARC5_MWD_010PDP_60B", "A10L_UNC", u"Resistivity - Attenuation", u"ohm.m")
parameterDict.update({'A10L_UNC' : Parameter(name='A10L_UNC',bname='a10l_unc',type='Variable',family='Resistivity - Attenuation',unit='ohm.m',value='Zabazaba - 1X.ARC5_MWD_010PDP_60B.A10L_UNC',mode='In',description='Required ARC channel',group='',min='',max='',list='',enable='True',iscombocheckbox='False',isused='True')})
#Type:Variable
#BName:a16l_unc
#Family:Resistivity - Attenuation
#Unit:ohm.m
#Mode:In
#Description:Required ARC channel
A16L_UNC = Variable("Zabazaba - 1X", "ARC5_MWD_010PDP_60B", "A16L_UNC", u"Resistivity - Attenuation", u"ohm.m")
parameterDict.update({'A16L_UNC' : Parameter(name='A16L_UNC',bname='a16l_unc',type='Variable',family='Resistivity - Attenuation',unit='ohm.m',value='Zabazaba - 1X.ARC5_MWD_010PDP_60B.A16L_UNC',mode='In',description='Required ARC channel',group='',min='',max='',list='',enable='True',iscombocheckbox='False',isused='True')})
#Type:Variable
#BName:a22l_unc
#Family:Resistivity - Attenuation
#Unit:ohm.m
#Mode:In
#Description:Required ARC channel
A22L_UNC = Variable("Zabazaba - 1X", "ARC5_MWD_010PDP_60B", "A22L_UNC", u"Resistivity - Attenuation", u"ohm.m")
parameterDict.update({'A22L_UNC' : Parameter(name='A22L_UNC',bname='a22l_unc',type='Variable',family='Resistivity - Attenuation',unit='ohm.m',value='Zabazaba - 1X.ARC5_MWD_010PDP_60B.A22L_UNC',mode='In',description='Required ARC channel',group='',min='',max='',list='',enable='True',iscombocheckbox='False',isused='True')})
#Type:Variable
#BName:a28l_unc
#Family:Resistivity - Attenuation
#Unit:ohm.m
#Mode:In
#Description:Required ARC channel
A28L_UNC = Variable("Zabazaba - 1X", "ARC5_MWD_010PDP_60B", "A28L_UNC", u"Resistivity - Attenuation", u"ohm.m")
parameterDict.update({'A28L_UNC' : Parameter(name='A28L_UNC',bname='a28l_unc',type='Variable',family='Resistivity - Attenuation',unit='ohm.m',value='Zabazaba - 1X.ARC5_MWD_010PDP_60B.A28L_UNC',mode='In',description='Required ARC channel',group='',min='',max='',list='',enable='True',iscombocheckbox='False',isused='True')})
#Type:Variable
#BName:a34l_unc
#Family:Resistivity - Attenuation
#Unit:ohm.m
#Mode:In
#Description:Required ARC channel
A34L_UNC = Variable("Zabazaba - 1X", "ARC5_MWD_010PDP_60B", "A34L_UNC", u"Resistivity - Attenuation", u"ohm.m")
parameterDict.update({'A34L_UNC' : Parameter(name='A34L_UNC',bname='a34l_unc',type='Variable',family='Resistivity - Attenuation',unit='ohm.m',value='Zabazaba - 1X.ARC5_MWD_010PDP_60B.A34L_UNC',mode='In',description='Required ARC channel',group='',min='',max='',list='',enable='True',iscombocheckbox='False',isused='True')})
#Type:Variable
#BName:p10l_unc
#Family:Resistivity - Phase Shift
#Unit:ohm.m
#Mode:In
#Description:Required ARC channel
P10L_UNC = Variable("Zabazaba - 1X", "ARC5_MWD_010PDP_60B", "P10L_UNC", u"Resistivity - Phase Shift", u"ohm.m")
parameterDict.update({'P10L_UNC' : Parameter(name='P10L_UNC',bname='p10l_unc',type='Variable',family='Resistivity - Phase Shift',unit='ohm.m',value='Zabazaba - 1X.ARC5_MWD_010PDP_60B.P10L_UNC',mode='In',description='Required ARC channel',group='',min='',max='',list='',enable='True',iscombocheckbox='False',isused='True')})
#Type:Variable
#BName:p16l_unc
#Family:Resistivity - Phase Shift
#Unit:ohm.m
#Mode:In
#Description:Required ARC channel
P16L_UNC = Variable("Zabazaba - 1X", "ARC5_MWD_010PDP_60B", "P16L_UNC", u"Resistivity - Phase Shift", u"ohm.m")
parameterDict.update({'P16L_UNC' : Parameter(name='P16L_UNC',bname='p16l_unc',type='Variable',family='Resistivity - Phase Shift',unit='ohm.m',value='Zabazaba - 1X.ARC5_MWD_010PDP_60B.P16L_UNC',mode='In',description='Required ARC channel',group='',min='',max='',list='',enable='True',iscombocheckbox='False',isused='True')})
#Type:Variable
#BName:p22l_unc
#Family:Resistivity - Phase Shift
#Unit:ohm.m
#Mode:In
#Description:Required ARC channel
P22L_UNC = Variable("Zabazaba - 1X", "ARC5_MWD_010PDP_60B", "P22L_UNC", u"Resistivity - Phase Shift", u"ohm.m")
parameterDict.update({'P22L_UNC' : Parameter(name='P22L_UNC',bname='p22l_unc',type='Variable',family='Resistivity - Phase Shift',unit='ohm.m',value='Zabazaba - 1X.ARC5_MWD_010PDP_60B.P22L_UNC',mode='In',description='Required ARC channel',group='',min='',max='',list='',enable='True',iscombocheckbox='False',isused='True')})
#Type:Variable
#BName:p28l_unc
#Family:Resistivity - Phase Shift
#Unit:ohm.m
#Mode:In
#Description:Required ARC channel
P28L_UNC = Variable("Zabazaba - 1X", "ARC5_MWD_010PDP_60B", "P28L_UNC", u"Resistivity - Phase Shift", u"ohm.m")
parameterDict.update({'P28L_UNC' : Parameter(name='P28L_UNC',bname='p28l_unc',type='Variable',family='Resistivity - Phase Shift',unit='ohm.m',value='Zabazaba - 1X.ARC5_MWD_010PDP_60B.P28L_UNC',mode='In',description='Required ARC channel',group='',min='',max='',list='',enable='True',iscombocheckbox='False',isused='True')})
#Type:Variable
#BName:p34l_unc
#Family:Resistivity - Phase Shift
#Unit:ohm.m
#Mode:In
#Description:Required ARC channel
P34L_UNC = Variable("Zabazaba - 1X", "ARC5_MWD_010PDP_60B", "P34L_UNC", u"Resistivity - Phase Shift", u"ohm.m")
parameterDict.update({'P34L_UNC' : Parameter(name='P34L_UNC',bname='p34l_unc',type='Variable',family='Resistivity - Phase Shift',unit='ohm.m',value='Zabazaba - 1X.ARC5_MWD_010PDP_60B.P34L_UNC',mode='In',description='Required ARC channel',group='',min='',max='',list='',enable='True',iscombocheckbox='False',isused='True')})
#Type:Variable
#BName:rt_in
#Family:Resistivity - Attenuation
#Unit:ohm.m
#Mode:In
#Description:Known Rt for mode 1&2
RT_IN = Variable("Zabazaba - 1X", "ARC5_MWD_010PDP_60B", "A34H_UNC", u"Resistivity - Attenuation", u"ohm.m")
parameterDict.update({'RT_IN' : Parameter(name='RT_IN',bname='rt_in',type='Variable',family='Resistivity - Attenuation',unit='ohm.m',value='Zabazaba - 1X.ARC5_MWD_010PDP_60B.A34H_UNC',mode='In',description='Known Rt for mode 1&2',group='',min='',max='',list='',enable='True',iscombocheckbox='False',isused='True')})
#Type:Variable
#BName:rm_in
#Family:Mud Resistivity
#Unit:ohm.m
#Mode:In
#Description:Known Rm for mode 2,3,5
RM_IN = Variable("Zabazaba - 1X", "ARC5_MWD_010PDP_60B", "BHRM", u"Mud Resistivity", u"ohm.m")
parameterDict.update({'RM_IN' : Parameter(name='RM_IN',bname='rm_in',type='Variable',family='Mud Resistivity',unit='ohm.m',value='Zabazaba - 1X.ARC5_MWD_010PDP_60B.BHRM',mode='In',description='Known Rm for mode 2\,3\,5',group='',min='',max='',list='',enable='True',iscombocheckbox='False',isused='True')})
#Type:Variable
#BName:rt_invert2
#Family:Array Resistivity
#Unit:OHMM
#Mode:Out
#Description:inverted formation resistivity
#Format:auto
RT_INVERT2 = Variable("Zabazaba - 1X", "ARC5_MWD_010PDP_60B", "RT_INVERT2", u"Array Resistivity", u"OHMM")
RT_INVERT2.setGroupName("ARCeCaliper")
parameterDict.update({'RT_INVERT2' : Parameter(name='RT_INVERT2',bname='rt_invert2',type='Variable',family='Array Resistivity',unit='OHMM',value='Zabazaba - 1X.ARC5_MWD_010PDP_60B.RT_INVERT2',mode='Out',description='inverted formation resistivity',group='',min='',max='',list='',enable='True',iscombocheckbox='False',isused='True')})
#Type:Variable
#BName:dh_invert2
#Family:Borehole Radius
#Unit:IN
#Mode:Out
#Description:inverted hole diameter
#Format:auto
DH_INVERT2 = Variable("Zabazaba - 1X", "ARC5_MWD_010PDP_60B", "DH_INVERT2", u"Borehole Radius", u"IN")
DH_INVERT2.setGroupName("ARCeCaliper")
parameterDict.update({'DH_INVERT2' : Parameter(name='DH_INVERT2',bname='dh_invert2',type='Variable',family='Borehole Radius',unit='IN',value='Zabazaba - 1X.ARC5_MWD_010PDP_60B.DH_INVERT2',mode='Out',description='inverted hole diameter',group='',min='',max='',list='',enable='True',iscombocheckbox='False',isused='True')})
#Type:Variable
#BName:rm_invert2
#Family:Array Resistivity
#Unit:OHMM
#Mode:Out
#Description:inverted mud resistivity
#Format:auto
RM_INVERT2 = Variable("Zabazaba - 1X", "ARC5_MWD_010PDP_60B", "RM_INVERT2", u"Array Resistivity", u"OHMM")
RM_INVERT2.setGroupName("ARCeCaliper")
parameterDict.update({'RM_INVERT2' : Parameter(name='RM_INVERT2',bname='rm_invert2',type='Variable',family='Array Resistivity',unit='OHMM',value='Zabazaba - 1X.ARC5_MWD_010PDP_60B.RM_INVERT2',mode='Out',description='inverted mud resistivity',group='',min='',max='',list='',enable='True',iscombocheckbox='False',isused='True')})
#Type:Variable
#BName:epsr_invert2
#Family:Permittivity
#Unit:unitless
#Mode:Out
#Description:inverted tool eccentricity
#Format:auto
EPSR_INVERT2 = Variable("Zabazaba - 1X", "ARC5_MWD_010PDP_60B", "EPSR_INVERT2", u"Permittivity", u"unitless")
EPSR_INVERT2.setGroupName("ARCeCaliper")
parameterDict.update({'EPSR_INVERT2' : Parameter(name='EPSR_INVERT2',bname='epsr_invert2',type='Variable',family='Permittivity',unit='unitless',value='Zabazaba - 1X.ARC5_MWD_010PDP_60B.EPSR_INVERT2',mode='Out',description='inverted tool eccentricity',group='',min='',max='',list='',enable='True',iscombocheckbox='False',isused='True')})
#Type:Variable
#BName:err_invert2
#Family:
#Unit:%
#Mode:Out
#Description:inverted error in %
#Format:auto
ERR_INVERT2 = Variable("Zabazaba - 1X", "ARC5_MWD_010PDP_60B", "ERR_INVERT2", u"", u"%")
ERR_INVERT2.setGroupName("ARCeCaliper")
parameterDict.update({'ERR_INVERT2' : Parameter(name='ERR_INVERT2',bname='err_invert2',type='Variable',family='',unit='%',value='Zabazaba - 1X.ARC5_MWD_010PDP_60B.ERR_INVERT2',mode='Out',description='inverted error in %',group='',min='',max='',list='',enable='True',iscombocheckbox='False',isused='True')})
#Type:Variable
#BName:qc_invert2
#Family:
#Unit:
#Mode:Out
#Description:inversion QC
#Format:auto
QC_INVERT2 = Variable("Zabazaba - 1X", "ARC5_MWD_010PDP_60B", "QC_INVERT2", u"", u"")
QC_INVERT2.setGroupName("ARCeCaliper")
parameterDict.update({'QC_INVERT2' : Parameter(name='QC_INVERT2',bname='qc_invert2',type='Variable',family='',unit='',value='Zabazaba - 1X.ARC5_MWD_010PDP_60B.QC_INVERT2',mode='Out',description='inversion QC',group='',min='',max='',list='',enable='True',iscombocheckbox='False',isused='True')})
#Type:Variable
#BName:rt_invert4
#Family:Array Resistivity
#Unit:OHMM
#Mode:Out
#Description:inverted formation resistivity
#Format:auto
RT_INVERT4 = Variable("Zabazaba - 1X", "ARC5_MWD_010PDP_60B", "RT_INVERT4", u"Array Resistivity", u"OHMM")
RT_INVERT4.setGroupName("ARCeCaliper")
parameterDict.update({'RT_INVERT4' : Parameter(name='RT_INVERT4',bname='rt_invert4',type='Variable',family='Array Resistivity',unit='OHMM',value='Zabazaba - 1X.ARC5_MWD_010PDP_60B.RT_INVERT4',mode='Out',description='inverted formation resistivity',group='',min='',max='',list='',enable='True',iscombocheckbox='False',isused='True')})
#Type:Variable
#BName:dh_invert4
#Family:Borehole Radius
#Unit:IN
#Mode:Out
#Description:inverted hole diameter
#Format:auto
DH_INVERT4 = Variable("Zabazaba - 1X", "ARC5_MWD_010PDP_60B", "DH_INVERT4", u"Borehole Radius", u"IN")
DH_INVERT4.setGroupName("ARCeCaliper")
parameterDict.update({'DH_INVERT4' : Parameter(name='DH_INVERT4',bname='dh_invert4',type='Variable',family='Borehole Radius',unit='IN',value='Zabazaba - 1X.ARC5_MWD_010PDP_60B.DH_INVERT4',mode='Out',description='inverted hole diameter',group='',min='',max='',list='',enable='True',iscombocheckbox='False',isused='True')})
#Type:Variable
#BName:rm_invert4
#Family:Array Resistivity
#Unit:OHMM
#Mode:Out
#Description:inverted mud resistivity
#Format:auto
RM_INVERT4 = Variable("Zabazaba - 1X", "ARC5_MWD_010PDP_60B", "RM_INVERT4", u"Array Resistivity", u"OHMM")
RM_INVERT4.setGroupName("ARCeCaliper")
parameterDict.update({'RM_INVERT4' : Parameter(name='RM_INVERT4',bname='rm_invert4',type='Variable',family='Array Resistivity',unit='OHMM',value='Zabazaba - 1X.ARC5_MWD_010PDP_60B.RM_INVERT4',mode='Out',description='inverted mud resistivity',group='',min='',max='',list='',enable='True',iscombocheckbox='False',isused='True')})
#Type:Variable
#BName:epsr_invert4
#Family:Permittivity
#Unit:unitless
#Mode:Out
#Description:inverted tool eccentricity
#Format:auto
EPSR_INVERT4 = Variable("Zabazaba - 1X", "ARC5_MWD_010PDP_60B", "EPSR_INVERT4", u"Permittivity", u"unitless")
EPSR_INVERT4.setGroupName("ARCeCaliper")
parameterDict.update({'EPSR_INVERT4' : Parameter(name='EPSR_INVERT4',bname='epsr_invert4',type='Variable',family='Permittivity',unit='unitless',value='Zabazaba - 1X.ARC5_MWD_010PDP_60B.EPSR_INVERT4',mode='Out',description='inverted tool eccentricity',group='',min='',max='',list='',enable='True',iscombocheckbox='False',isused='True')})
#Type:Variable
#BName:err_invert4
#Family:
#Unit:%
#Mode:Out
#Description:inverted error in %
#Format:auto
ERR_INVERT4 = Variable("Zabazaba - 1X", "ARC5_MWD_010PDP_60B", "ERR_INVERT4", u"", u"%")
ERR_INVERT4.setGroupName("ARCeCaliper")
parameterDict.update({'ERR_INVERT4' : Parameter(name='ERR_INVERT4',bname='err_invert4',type='Variable',family='',unit='%',value='Zabazaba - 1X.ARC5_MWD_010PDP_60B.ERR_INVERT4',mode='Out',description='inverted error in %',group='',min='',max='',list='',enable='True',iscombocheckbox='False',isused='True')})
#Type:Variable
#BName:qc_invert4
#Family:
#Unit:
#Mode:Out
#Description:inversion QC
#Format:auto
QC_INVERT4 = Variable("Zabazaba - 1X", "ARC5_MWD_010PDP_60B", "QC_INVERT4", u"", u"")
QC_INVERT4.setGroupName("ARCeCaliper")
parameterDict.update({'QC_INVERT4' : Parameter(name='QC_INVERT4',bname='qc_invert4',type='Variable',family='',unit='',value='Zabazaba - 1X.ARC5_MWD_010PDP_60B.QC_INVERT4',mode='Out',description='inversion QC',group='',min='',max='',list='',enable='True',iscombocheckbox='False',isused='True')})
#DeclarationsEnd
import sys
import os
import ctypes
from time import *

start = time()
from ctypes import *
from _ctypes import FreeLibrary

# DLL directory, name and path
#dllName = "ARCEcaliper_2014.dll"
#dllName = "vdiinv_0515_sf.DLL.dll" 
dllName = "vdiinv_0827_sfdll.dll" 
if "64 bit" in sys.version:
	print "64 bit system"
	Loaded_DLL = windll.LoadLibrary(db.dirUser()+"\\External_DLLs\\x64\\" + dllName)
#else:
	#Loaded_DLL = windll.LoadLibrary(db.dirUser()+"\\External_DLLs\\Win32\\" + dllName)
#print "(db.dirUser()+"\\External_DLLs\\x64\\" + dllName)"
# It's 64 bit system now, so no choice needed, this is where to find DLL
#if "64 bit" in sys.version:   
#	Loaded_DLL = windll.LoadLibrary(db.dirUser()+"\\External_DLLs\\x64\\" + dllName)
#else:
#Loaded_DLL = windll.LoadLibrary(db.dirUser()+"\\External_DLLs\\Win32\\" + dllName)
print Loaded_DLL
print Loaded_DLL._handle

# this is where to find the lookup table
#Data_FilePath = db.dirUser()+"\\Dll_ARCEcaliper_new\\ARCeCaliper_data\\"
Data_FilePath = db.dirUser()+"\\External_DLLs\\Venus_diinv_2015\\Vdiinv_table\\"
pathSize = len(Data_FilePath)
print "Data file path: %s" %(Data_FilePath)

# read tool/ mud/ channel from input
toolList = {'ARC475':4, 'Venus475':5, 'Venus675':6, 'ARC675':7}
type = toolList[tool_type]
print tool_type, type
	
if mud_type == "WBM":
	mudtype = 0
	mode = 2
else:
	mudtype = 1
	mode = 1
print "Mud is: ", mud_type
print "mode is: ", mode

channels = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
#if type in (4,5): 
ChannelList = ['P10H','P16H', 'P22H','P28H', 'P34H', 'A10H', 'A16H', 'A22H', 'A28H', 'A34H', \
               'P10L','P16L', 'P22L','P28L', 'P34L', 'A10L', 'A16L', 'A22L', 'A28L', 'A34L']
#elif type in (6,7):
#ChannelList = ['P16H', 'P22H','P28H', 'P34H', 'P40H', 'A16H', 'A22H', 'A28H', 'A34H', 'A40H', \
               #'P16L', 'P22L','P28L', 'P34L', 'P40L', 'A16L', 'A22L', 'A28L', 'A34L', 'A40L']
#else:
#print "Tool not available"	              	
nch=20     # number of channels

nch_desire = 20
i_auto=0
misfit_epsr = 0.001  # just kept it, no real purpose

# variables
idtool      = c_int(type)
isOBM       = c_int(mudtype)
nx          = c_int(mode)
bs          = c_float(bit_size)
mf_epsr      = c_float(misfit_epsr)
path_c = c_char_p(Data_FilePath)
path_size = c_int(pathSize)
nchannel = c_int(nch)
nch_desire_c = c_int(nch_desire)             
i_auto_c = c_int(i_auto)


FChannel = c_float * nch
IChannel = c_int * nch
# define c_type index_channel and corresponding weight
weight_channel = FChannel()
index_channel = IChannel()
for c_index in range(nch):
	index_channel[c_index] = c_int(channels[c_index])
	weight_channel[c_index] = c_float(1)
	print "index of selected channel is {0}, weight is {1}".format(index_channel[c_index],weight_channel[c_index])
print index_channel[:], weight_channel[:]

# find total number of measurements (from depth)
numberOfmeas = DEPT.referenceSize()   # total number is already pointer
# print DEPT.referenceValue()
print "total number is {}".format(numberOfmeas) 
totalnumber = c_int(numberOfmeas)

# create types and pointers
F20 = c_float*20 
f100k   = c_float * numberOfmeas
f100k_c = f100k * nch

## create all the 1D, 2D arrays
rad2_c = F20()
depth_c = f100k()
Measure  = f100k_c()
rt_input_c = f100k()
rm_input_c = f100k()
dh_input_c = f100k()
epsr_input_c = f100k()
inv_rt_c2 = f100k()
inv_dh_c2 = f100k()
inv_epsr_c2 = f100k()
inv_rm_c2 = f100k()
inv_err_c2 = f100k()
inv_qc_c2 = f100k()
inv_rt_c4 = f100k()
inv_dh_c4 = f100k()
inv_epsr_c4 = f100k()
inv_rm_c4 = f100k()
inv_err_c4 = f100k()
inv_qc_c4 = f100k()

# assign memory space for each row of output 2d matrix Measure
indexOfrows = 0
for indexOfrows in range(nch):
	Measure[indexOfrows] = f100k()

#Loop to read in the data	
tn=0
### Begin Automatic Generation Loop ###
loopSize = DEPT.referenceSize()
for loopIterator in xrange(loopSize):
	datasetIterator = loopIterator
	dept = DEPT.value(loopIterator)
	a10h_unc = A10H_UNC.value(loopIterator)
	a16h_unc = A16H_UNC.value(loopIterator)
	a22h_unc = A22H_UNC.value(loopIterator)
	a28h_unc = A28H_UNC.value(loopIterator)
	a34h_unc = A34H_UNC.value(loopIterator)
	p10h_unc = P10H_UNC.value(loopIterator)
	p16h_unc = P16H_UNC.value(loopIterator)
	p22h_unc = P22H_UNC.value(loopIterator)
	p28h_unc = P28H_UNC.value(loopIterator)
	p34h_unc = P34H_UNC.value(loopIterator)
	a10l_unc = A10L_UNC.value(loopIterator)
	a16l_unc = A16L_UNC.value(loopIterator)
	a22l_unc = A22L_UNC.value(loopIterator)
	a28l_unc = A28L_UNC.value(loopIterator)
	a34l_unc = A34L_UNC.value(loopIterator)
	p10l_unc = P10L_UNC.value(loopIterator)
	p16l_unc = P16L_UNC.value(loopIterator)
	p22l_unc = P22L_UNC.value(loopIterator)
	p28l_unc = P28L_UNC.value(loopIterator)
	p34l_unc = P34L_UNC.value(loopIterator)
	rt_in = RT_IN.value(loopIterator)
	rm_in = RM_IN.value(loopIterator)
	rt_invert2 = MissingValue
	dh_invert2 = MissingValue
	rm_invert2 = MissingValue
	epsr_invert2 = MissingValue
	err_invert2 = MissingValue
	qc_invert2 = MissingValue
	rt_invert4 = MissingValue
	dh_invert4 = MissingValue
	rm_invert4 = MissingValue
	epsr_invert4 = MissingValue
	err_invert4 = MissingValue
	qc_invert4 = MissingValue
	### Automatic Generation Loop End ###
	#rad2 and rps2 have a list of 5 values per depth
	#print "error"
	#if type in (4,5): 
	rad2_c = [p10h_unc,p16h_unc, p22h_unc, p28h_unc, p34h_unc, a10h_unc, a16h_unc, a22h_unc, a28h_unc, a34h_unc, \
	        p10l_unc,p16l_unc, p22l_unc, p28l_unc, p34l_unc, a10l_unc, a16l_unc, a22l_unc, a28l_unc, a34l_unc]
	#elif type in (6,7):
	#rad2_c = [p16h_unc, p22h_unc, p28h_unc, p34h_unc, p40h_unc, a16h_unc, a22h_unc, a28h_unc, a34h_unc, a40h_unc, \
	        #p16l_unc, p22l_unc, p28l_unc, p34l_unc, p40l_unc, a16l_unc, a22l_unc, a28l_unc, a34l_unc, a40l_unc]
	depth_c[tn] = c_float(dept)
	rt_input_c[tn]= c_float(rt_in)
	rm_input_c[tn]= c_float(rm_in)      #   IF rm_in is an input channel
	if isOBM.value == 1:
		rm_input_c[tn]= 1000.0
	dh_input_c[tn]= c_float(bit_size)
	epsr_input_c[tn]= 10.0
	for k in range(nch):
		IndexofPick = index_channel[k]
		Measure[k][tn] = c_float(rad2_c[IndexofPick-1])
	#print Measure[0][tn],Measure[16][tn],Measure[17][tn],Measure[18][tn],Measure[19][tn],Measure[10][tn]
	tn = tn+1
	#if tn>5:
		#break

#Check the mud type
	### Begin Automatic Generation EndLoop ###
	RT_INVERT2.setValue(loopIterator, rt_invert2)
	DH_INVERT2.setValue(loopIterator, dh_invert2)
	RM_INVERT2.setValue(loopIterator, rm_invert2)
	EPSR_INVERT2.setValue(loopIterator, epsr_invert2)
	ERR_INVERT2.setValue(loopIterator, err_invert2)
	QC_INVERT2.setValue(loopIterator, qc_invert2)
	RT_INVERT4.setValue(loopIterator, rt_invert4)
	DH_INVERT4.setValue(loopIterator, dh_invert4)
	RM_INVERT4.setValue(loopIterator, rm_invert4)
	EPSR_INVERT4.setValue(loopIterator, epsr_invert4)
	ERR_INVERT4.setValue(loopIterator, err_invert4)
	QC_INVERT4.setValue(loopIterator, qc_invert4)
RT_INVERT2.save(True)
DH_INVERT2.save(True)
RM_INVERT2.save(True)
EPSR_INVERT2.save(True)
ERR_INVERT2.save(True)
QC_INVERT2.save(True)
RT_INVERT4.save(True)
DH_INVERT4.save(True)
RM_INVERT4.save(True)
EPSR_INVERT4.save(True)
ERR_INVERT4.save(True)
QC_INVERT4.save(True)
### Automatic Generation EndLoop End ###
if isOBM.value == 0:
	print "Mud Type: WBM"
else:
	print "Mud Type: OBM"

# check if the inputs are loading correctly
print "tool_type: %s, bit_size = %s in." % (tool_type, bs.value)
if nx.value == 1: 
	print "Only invert for Rt, Epsr. Please select input channel for RT_IN and value for Hd, Rm"
elif nx.value == 2:
	print "Invert for Dh, Rm, Rt and Epsr. Please select input channel for RT_IN"
#elif nx.value == 3:
	#print "Invert for Dh, Rt, Rm"
#elif nx.value == 4:
	#print "Invert for Dh, Rt. Please select input channel for Rm"
#elif nx.value == 5:
	#print "Invert for Dh, Rm, Rt, Ecc."
#else:
	#print "Invert for Dh, Rt, Rt, EPsr, Please select input channel for Rm"
print totalnumber,nch_desire_c,i_auto_c
print "Program is running ..."
#print " rt_input_c",rt_input_c[:]
#print " dh_input_c, rm_input_c, rt_input_c, epsr_input_c",dh_input_c[:], rm_input_c[:], rt_input_c[:], epsr_input_c[:]
#print nchannel
#print totalnumber
##sys.exit()
#print rad2_c[1],rad2_c[19]
#print rad2_c[1][1],rad2_c[1][2],rad2_c[1][3],rad2_c[1][4],rad2_c[1][5],rad2_c[1][6],rad2_c[1][7],rad2_c[1][8],rad2_c[1][9]
#print rad2_c[2][1],rad2_c[2][2],rad2_c[2][3],rad2_c[2][4],rad2_c[2][5],rad2_c[2][6],rad2_c[2][7],rad2_c[2][8],rad2_c[2][9]
## call the function from dll
#print idtool,bs,nx,nchannel
#print index_channel[:],weight_channel[:]
##print bs,nx,nchannel,index_channel,weight_channel

##print depth_c[:]
#print path_c,path_size,isOBM   
#,mf_epsr
#print Measure[19][:]
h = Loaded_DLL._handle
Loaded_DLL.vdie_inv_v003(byref(idtool),byref(bs), byref(nx), byref(nchannel), index_channel, weight_channel, \
                       byref(totalnumber), depth_c, Measure, byref(nch_desire_c), byref(i_auto_c), dh_input_c, rm_input_c, rt_input_c, epsr_input_c,\
                       path_c, byref(path_size), byref(isOBM),byref(mf_epsr), inv_dh_c2, inv_rm_c2, inv_rt_c2,\
                       inv_epsr_c2, inv_err_c2,inv_qc_c2, inv_dh_c4, inv_rm_c4, inv_rt_c4,inv_epsr_c4,\
                       inv_err_c4,inv_qc_c4,Data_FilePath)


#Read in results
for k in range(tn):
	RT_INVERT2.setValue(k, inv_rt_c2[k])
	DH_INVERT2.setValue(k, inv_dh_c2[k])
	EPSR_INVERT2.setValue(k, inv_epsr_c2[k])
	RM_INVERT2.setValue(k, inv_rm_c2[k])
	ERR_INVERT2.setValue(k, inv_err_c2[k])
	QC_INVERT2.setValue(k, inv_qc_c2[k])
	RT_INVERT4.setValue(k, inv_rt_c4[k])
	DH_INVERT4.setValue(k, inv_dh_c4[k])
	EPSR_INVERT4.setValue(k, inv_epsr_c4[k])
	RM_INVERT4.setValue(k, inv_rm_c4[k])
	ERR_INVERT4.setValue(k, inv_err_c4[k])
	QC_INVERT4.setValue(k, inv_qc_c4[k])
	#if (inv_dh_c[k] > bit_size + 5) | (inv_dh_c[k] < bit_size - 10):
		#DH_INVERT.setValue(k, bit_size)
	#if (k < 4):
		#print 'k  inv_rt  inv_dh', k, inv_rt_c[k], inv_dh_c[k]
#print ('Data has been filtered')
RT_INVERT2.save()
DH_INVERT2.save()
EPSR_INVERT2.save()
RM_INVERT2.save()
ERR_INVERT2.save()
QC_INVERT2.save()
RT_INVERT4.save()
DH_INVERT4.save()
EPSR_INVERT4.save()
RM_INVERT4.save()
ERR_INVERT4.save()
QC_INVERT4.save()

FreeLibrary(h)
end = time()

if end - start >= 1:
	print "Duration (hh:mm:ss):", strftime("%H:%M:%S", gmtime(end - start))
else:
	print "Duration (s):", end - start

__doc__ = """Full Table Dielectric Inversion

History
Version 1 diinv_2015_v1    use   fortran compiled  vdiinv_0515_sf.dll
mode 1: invert for Rt and Epsr only, for OBM, or WBM
mode 2: invert for Rt, Epsr, Hd and RM, for WBM only, no OBM
two frequency sets 2MHz and 400kHz
Inversion QC
 

Inversion processing to obtain Rt and Epsr in WBM or OBM using uncorrected ARC (2 MHz and 400kHz) PS and AT measurements. 

Assumptions:
Circular borehole, No invasion, vertical well, tool in homogeneous formation, mud dielectric constant=80,  no eccentricity

Output channels are written to variable group ARCeCaliper (can be modified via Output variable options tab).
#####################################################
INPUTS:

Tool_type:  ARC 5, 6  or Venus 5 or 6  (Four different tables can be looked up from)

mud_type: to define WBM or OBM. 

Bit_size:  Hole size needs to be input as reference, not necessarily to be precise. Use largest element in BHA (BS or Reamer size), this is needed especially for mode 1, because it will be regarded as accurate. The inaccuracy will cause errors for other variables' inversion

Channels to double check:

DEPT: depth channel
ARC uncorrected attenuation and phase shift channels (10 ATT and 10 PS). Channel names listed under "Value" must be exactly the same as in the data set. 
A16H_UNC, A22H_UNC, A28H_UNC, A34H_UNC, A40H_UNC
P16H_UNC, P22H_UNC, P28H_UNC, P34H_UNC, P40H_UNC
A16L_UNC, A22L_UNC, A28L_UNC, A34L_UNC, A40L_UNC
P16L_UNC, P22L_UNC, P28L_UNC, P34L_UNC, P40L_UNC

If for Arc 5 or Venus 5, A10H_UNC,A10L_UNC,P10H_UNC,P10L_UNC

Rm_IN: Mud resistivity data channel.   BHRM or user defined. Needed for mode 1
Hd_IN:  Bit size input. Needed for mode 1

#####################################################
OUTPUTS

RT_INVERT: inverted formation resistivity (ohm.m)
HD_INVERT: inverted hole diameter (in)
RM_INVERT: inverted mud resistivity (ohm.m)
ERR_INVERT: inversion error (%)
EPSR_INVERT: inversion dielectric constant
QC_INVERT: inversion QC (unitless, 0-1)
"""
__author__ = """Tina Zhao(tzhao)"""
__date__ = """2015-06-25"""
__version__ = """1.0"""
__group__ = """ARCeCaliper"""
__suffix__ = """"""
__prefix__ = """"""