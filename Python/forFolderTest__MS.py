# for MS windows
romanNumeralMap = (('M',1000),
	('CM',900),
	('D',500),
	('CD',400),
	('C',100),
	('XC',90),
	('L',50),
	('XL',40),
	('X',10),
	('IX',9),
	('V',5),
	('IV',4),
	('I',1))
def toRoman(n):
	result = ""
	for numeral, integer in romanNumeralMap:
		while n >= integer:
		 	result += numeral
		 	n -= integer
	return result
def fromRoman(s):
	result = 0
	index = 0
	for numeral, integer in romanNumeralMap:
		while s[index:index+len(numeral)] == numeral:
		 	result += integer
		 	index += len(numeral)
	return result
# print toRoman(1356)
# print fromRoman('MCMLXXII')

import os

os.system('E:\\TMP\\Vid\\zipFile.bat') # zip all wanted folders
flog = open('C:\Users\SWUNQ\Desktop\pylog2.txt','a')
flog.write('===============================================\n')
flog.write('Original file name\t\t->\tNew file name\n')
flog.write('-----------------------------------------------\n')
finalPath = 'E:\\TMP\\Vid\\Test\\'
os.chdir(finalPath)
fileLs=os.listdir('.')
for oldName in fileLs:
	if not(oldName == 'renameFile.py' or 
		oldName == 'ziplog.txt' or oldName == 'zipFile.bat'):
		fName = oldName[:-7]    # file name without type name
		typeName = oldName[-7:] # file type, in this case, it's like .7z.00*
		flog.write(oldName)
		if len(fName)<8:
			flog.write('\t\t')
		flog.write('\t\t->\t')
		pre_f = fName[:-3]
		pre_f = pre_f.replace('-','_')
		pre_f = pre_f.replace(' ','_')
		pre_f = pre_f.replace('(720P)','HD_')		
		if not (pre_f.find('FHD') == -1):
			pre_f = pre_f[:5] + 'xx' + pre_f[5:]
		elif not (pre_f.find('HD') == -1):
			pre_f = pre_f[:4] + 'xx' + pre_f[4:]
		else:
			pre_f = pre_f[:1] + 'xx' + pre_f[1:]
		sNum = fName[-3:]
		s = sNum[-1]
		if s.isdigit():
			suf_f = toRoman(int(sNum))
			newName = pre_f + suf_f + typeName
			os.rename(os.path.join(finalPath,oldName), os.path.join(finalPath,newName))
			flog.write(newName)
		else:
			print fName, "No number here!\n"
		flog.write('\n')
			
flog.close()