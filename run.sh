#!/bin/bash

echo "Run handin template -- this is an example, make sure to edit this (e.g. rename template.tex to [YourLastName].tex)!"

# echo "Clearing/creating the plotting directory"
# if [ ! -d "plots" ]; then
#   mkdir plots
# fi
# rm -rf plots/*

# echo "Check if the sine movie exist"
# if [ -e sinemovie.mp4 ]; then
#   echo "Remove mp4 file"
#   rm sinemovie.mp4
# fi

# echo "Download image for in report..."
# if [ ! -e sine.png ]; then
#   wget home.strw.leidenuniv.nl/~daalen/Handin_files/sine.png 
# fi

# # Script that returns a plot
# echo "Run the first script ..."
# python3 sine.py

# # Script that pipes output to a file
# echo "Run the second script ..."
# python3 helloworld.py > helloworld.txt

# # Script that saves data to a file
# echo "Run the third script ..."
# python3 cos.py

# # Script that generates movie frames
# echo "Run the fourth script ..."
# python3 sinemovie.py

# # code that makes a movie of the movie frames
# ffmpeg -framerate 25 -pattern_type glob -i "plots/snap*.png" -s:v 640x480 -c:v libx264 -profile:v high -level 4.0 -crf 10 -tune animation -preset slow -pix_fmt yuv420p -r 25 -threads 0 -f mp4 sinemovie.mp4


echo "Running the script for question 1"
python3 1_poisson_distribution.py > poisson_distribution.txt

echo "Running the script for question 2"
python3 2_vandermonde_matrix.py > vandermonde_matrix.txt

echo "Generating the pdf"

pdflatex pranger.tex
bibtex pranger.aux
pdflatex pranger.tex
pdflatex pranger.tex


