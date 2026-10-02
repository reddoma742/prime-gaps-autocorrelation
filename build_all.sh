#!/usr/bin/env bash
# build_all.sh — reproduce everything from scratch.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

echo "=== 1/5  Compiling C++ pipelines ==="
cd src/cpp
g++ -O3 -std=c++17 -march=native -pthread \
    wheel_analysis_v9_1_complete.cpp -o v9_1
g++ -O3 -std=c++17 -march=native -pthread \
    wheel_analysis_v10_supplement_v2.cpp -o v10_v2
g++ -O3 -std=c++17 theta_powers.cpp -o theta_powers

echo "=== 2/5  Running C++ pipelines (~22 h) ==="
# Uncomment only if you need to regenerate data
# ./v9_1
# ./v10_v2
# ./theta_powers

echo "=== 3/5  Running Python analysis ==="
cd "$ROOT/src/python"
python fit_c.py
python delta_mean_verify.py
python compute_psi.py
python make_figures.py
python make_appendix_b.py

echo "=== 4/5  Assembling LaTeX package ==="
cd "$ROOT"
mkdir -p build
cp paper/manuscript.tex build/
cp paper/refs.bib       build/
cp paper/appendix_b_table.tex build/ 2>/dev/null || \
    cp src/python/appendix_b_table.tex build/
cp paper/cover_letter.tex build/
cp paper/figures/*.pdf  build/

echo "=== 5/5  Compiling PDF ==="
cd build
pdflatex -interaction=nonstopmode manuscript.tex
bibtex   manuscript
pdflatex -interaction=nonstopmode manuscript.tex
pdflatex -interaction=nonstopmode manuscript.tex

echo "Done. PDF at build/manuscript.pdf"