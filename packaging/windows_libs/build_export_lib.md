# Create export library

If F2Py is failing to find an msvcr90.dll

 cd <workdir>
 copy C:\WINDOWS\winsxs\amd64_microsoft.vc90.crt_1fc8b3b9a1e18e3b_9.0.30729.1_none_99b61f5e8371c1d4\msvcr90.dll .
 
 
 gendef.exe msvcr90.dll

dlltool.exe -D msvcr90.dll -l libmsvcr90.a -d msvcr90.def

copy libmsvcr90.a MINGW64_PATH\lib