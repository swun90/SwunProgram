@echo off
:chooseFunc
echo =============== MENU ===============
echo   1. connect free IP
echo   2. connect global/charge IP
echo   3. disconnect
echo   4. disconnect all
echo   q. quit
echo ===================================
set input=
set /p input=Please select a number: 
if /i '%input%'=='1' goto freeIP
if /i '%input%'=='2' goto globalIP
if /i '%input%'=='3' goto disconnect
if /i '%input%'=='4' goto disconnectAll
if /i '%input%'=='q' goto endD
echo wrong number, try again!
goto chooseFunc

:freeIP
echo linking free IP...
python pkuipgw.py -c --account=pkuipgwrc
goto end
:globalIP
echo set time (eg: 10m, 2h, 1:30(=1:30:00)):
set /p t=
if /i '%t:~-1,1%'=='m' (
	if %t:~,-1% LSS 10 (
		set time=00:0%t:~,-1%:00
	) else (
	set time=00:%t:~,-1%:00
	)
)
if /i '%t:~-1,1%'=='h' set time=%t:~,-1%:00:00
if '%t:~-3,1%'==':' (
	if '%t:~-6,1%'==':' (
		set time=%t%
	) else (
	set time=%t%:00
	)
)

echo linking global/charged IP...
echo And, do NOT close this window before switched!
python pkuipgw.py -ac --account=pkuipgwrc --time %time%
goto close
:disconnect
echo disconnecting...
python pkuipgw.py -d --account=pkuipgwrc
goto chooseFunc
:disconnectAll
echo disconnecting...
python pkuipgw.py -ad --account=pkuipgwrc
goto chooseFunc

rem python pkuipgw.py -ac --account=pkuipgwrc --time 0:10:00

goto end
:no
:end
echo Goodbye, work hard and have a nice day!
timeout 5

:endD
rem echo Goodbye, work hard and have a nice day!
:close