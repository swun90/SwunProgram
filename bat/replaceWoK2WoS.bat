@echo off
chcp 65001
setlocal enabledelayedexpansion
for /f "tokens=* delims=" %%i in (%1) do (
set var=%%i
echo !var:Web of Scienc=Web of Knowledg!>>new.txt
)
::move new.txt %1
echo successfully replaced words
pause
