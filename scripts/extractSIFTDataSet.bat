for %%a IN (%1\*.jpg) do convert %%a %1\%%~na.pgm
for %%b IN (%1\*.pgm) do "..\..\siftDemoV4\siftWin32.exe" <%%b >%1\%%~nbsift.txt