#!/bin/bash

echo "Running the script for question 1"
python3 1_poisson_distribution.py > poisson_distribution.txt

echo "Running the script for question 2"
python3 2_vandermonde_matrix.py > vandermonde_matrix.txt

echo "Generating the pdf"

pdflatex pranger.tex
bibtex pranger.aux
pdflatex pranger.tex
pdflatex pranger.tex