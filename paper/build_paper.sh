#!/bin/bash
# ===================================================================
#  build_paper.sh  -  compile manuscript.tex -> manuscript.pdf
#
#  Pour Overleaf :  ce fichier n'est pas necessaire, Overleaf compile
#  automatiquement.  Utilisez-le en local (Linux / macOS / WSL).
# ===================================================================
set -e
cd "$(dirname "$0")"

command -v pdflatex >/dev/null 2>&1 || {
    echo "pdflatex introuvable : installez TeX Live ou MiKTeX." >&2
    exit 1
}

echo "=== Passe 1 : manuscript ==="
pdflatex -interaction=nonstopmode -halt-on-error manuscript.tex

echo "=== BibTeX ==="
bibtex manuscript

echo "=== Passe 2 : manuscript ==="
pdflatex -interaction=nonstopmode -halt-on-error manuscript.tex

echo "=== Passe 3 : references croisees ==="
pdflatex -interaction=nonstopmode -halt-on-error manuscript.tex

if [ -f manuscript.pdf ]; then
    echo "OK -> manuscript.pdf"
else
    echo "ECHEC - lisez manuscript.log" >&2
    exit 1
fi
