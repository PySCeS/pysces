"%PYTHON%" setup.py build --compiler=mingw32 install --record=record.txt
if errorlevel 1 exit 1
