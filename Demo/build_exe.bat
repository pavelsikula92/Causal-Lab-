@echo off
rem Builds the standalone demo Demo\dist\CausalPhysics.exe (Windows; one file, no Python needed to run it).
rem Run from anywhere: the script moves to the repository root. Needs Python 3.9+ with tkinter (the python.org installer has it).
cd /d "%~dp0\.."
python Demo\causal_physics.py --selftest || (echo selftest FAILED - nothing built & exit /b 1)
python -m pip install --upgrade pyinstaller numpy
python -m PyInstaller --onefile --windowed --noconfirm --name CausalPhysics --distpath Demo\dist --workpath build\pyinstaller --specpath build\pyinstaller Demo\causal_physics.py
if exist Demo\dist\CausalPhysics.exe (echo Built: Demo\dist\CausalPhysics.exe) else (echo Build FAILED & exit /b 1)
