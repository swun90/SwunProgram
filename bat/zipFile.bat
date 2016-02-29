@echo off
set /a StartS=%time:~6,2%
set /a StartM=%time:~3,2%
echo start_time: %time% >> ziplog.txt
for /d %%i in (*) do 7z a E:\TMP\Vid\Shion\%%i %%i -v1000m -mx0 >> ziplog.txt
REM 7z a G:\Testing\Tes1 AB2 -v1000m -mx0 >> log.txt
set /a EndS=%time:~6,2%
set /a EndM=%time:~3,2%
echo ending_time: %time% >> ziplog.txt
set /a diffS_=%EndS%-%StartS%
set /a diffM_=%EndM%-%StartM%
echo whole_time: %diffM_%min-%diffS_%sec >> ziplog.txt
echo -------------------------------------------- >> ziplog.txt
echo -------------------------------------------- >> ziplog.txt
echo= >> ziplog.txt
