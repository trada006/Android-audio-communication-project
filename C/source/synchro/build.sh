cd code/data/
mpost -interaction=nonstopmode plot.mp
mpost -interaction=nonstopmode benchplot.mp 
cd ../..

pdflatex async.tex
bibtex async
pdflatex async.tex
pdflatex async.tex
