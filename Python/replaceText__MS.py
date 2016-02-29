# for MS windows
import os
import re
from glob import glob
srcNm=glob("savedrecs*.txt")
for i in range(len(srcNm)):
	print srcNm[i]
	f1=file(srcNm[i],"r")
	gNm="goal_try.txt"
	f2=file(gNm,"w")
	for s in f1.readlines():
		f2.write(s.replace("Web of Science", "Web of Knowledge"))

	f1.close()
	os.remove(srcNm[i])
	f2.close()
	os.rename(gNm, srcNm[i])
print "Text replaced!"