COMPILER=latexmk

all: $(NAME).tex
	$(COMPILER) -pdf -pdflatex="pdflatex -shell-escape -enable-write18" \
		-use-make -auxdir=file -outdir=file $(NAME) ;
	shopt -s extglob ;
	mv file/*.pdf file/pdf

init: 
	mkdir file ;
	mkdir file/pdf 
