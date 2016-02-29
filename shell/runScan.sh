#!/bin/sh
echo scan prepare log > scan_prepare.log
pwd >> scan_prepare.log
echo copy source files \(no log if success\) >> scan_prepare.log
scanDir="05 10 15"
echo without diamag term >> scan_prepare.log
for n in $scanDir; do
	cp -r source noDia_n$n >> scan_prepare.log
	cd noDia_n$n
	pwd >> scan_prepare.log
	sed -i "s/ZPERIOD\ \=\ 15/ZPERIOD\ \=\ $n/" ./data/BOUT.inp \
		>> scan_prepare.log
	sed -i 's/diamag\ \=\ true/diamag\ \=\ false/' ./data/BOUT.inp \
		>> scan_prepare.log
	bsub < bsubscript
	cd ..
done
echo -------------------- >> scan_prepare.log
echo with diamag term >> scan_prepare.log
for n in $scanDir; do
	cp -r source wDia_n$n >> scan_prepare.log
	cd wDia_n$n
	pwd >> scan_prepare.log
	sed -i "s/ZPERIOD\ \=\ 15/ZPERIOD\ \=\ $n/" ./data/BOUT.inp \
		>> scan_prepare.log
	bsub < bsubscript
	cd ..
done
echo finished! check scan_prepare.log for more infomations.
