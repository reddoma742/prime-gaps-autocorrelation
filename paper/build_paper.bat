@echo off
REM ===================================================================
REM  build_paper.bat  -  compile manuscript.tex -> manuscript.pdf
REM
REM  REQUIRES a LaTeX distribution (MiKTeX or TeX Live).
REM  It is NOT installed on this machine.  Install MiKTeX first:
REM      https://www.miktex.org/download/
REM  then run this file from inside the  paper\  folder.
REM ===================================================================
setlocal
cd /d "%~dp0"

where pdflatex >nul 2>&1
if errorlevel 1 (
  echo.
  echo   pdflatex introuvable.  Installez MiKTeX ou TeX Live d'abord.
  echo.
  pause
  exit /b 1
)

echo === Passe 1 : manuscript ===
pdflatex -interaction=nonstopmode -halt-on-error manuscript.tex

echo === BibTeX ===
bibtex manuscript

echo === Passe 2 : manuscript ===
pdflatex -interaction=nonstopmode -halt-on-error manuscript.tex

echo === Passe 3 : references ===
pdflatex -interaction=nonstopmode -halt-on-error manuscript.tex

echo.
if exist manuscript.pdf (
  echo   OK -> manuscript.pdf
) else (
  echo   ECHEC - lisez manuscript.log
)
echo.
pause
