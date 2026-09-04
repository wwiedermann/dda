@echo off
cd /d "%~dp0"
git push --force-with-lease origin main
echo.
echo Done! Press any key to close.
pause >nul
