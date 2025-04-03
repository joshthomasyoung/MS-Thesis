@echo off
setlocal enabledelayedexpansion

REM Define paths
set "input_folder=%cd%\Fairbanks INP"
set "output_folder=%cd%\Fairbanks OUT"
set "smarts_exe=%cd%\smarts295bat.exe"

REM Loop through files 2880 to 5784
for /L %%i in (2880,1,5784) do (
    REM Format the file number to match the naming convention (e.g., 2880.txt)
    set "filename=%%i"
    set "input_file=%input_folder%\!filename!.txt"
    set "output_file_base=%output_folder%\!filename!"

    REM Check if input file exists
    if exist "!input_file!" (
        echo Processing file: !input_file!

        REM Copy the input file as smarts295.inp.txt in the working directory
        copy /Y "!input_file!" smarts295.inp.txt

        REM Run SMARTS in batch mode
        "%smarts_exe%" > nul

        REM Move the output files to the output folder with the correct names
        move /Y smarts295.out.txt "!output_file_base!out.txt"
        move /Y smarts295.ext.txt "!output_file_base!ext.txt"

        REM Clean up
        del /Q smarts295.inp.txt
    ) else (
        echo Skipping file !input_file! (not found)
    )
)

echo All specified iterations completed successfully.
endlocal