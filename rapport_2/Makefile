MAIN_TEX = rapport

TEX = pdflatex
TEX_ARG = -interaction=nonstopmode -file-line-error

all: $(MAIN_TEX).pdf

clean:
	rm -f *.dvi *.eps *.log *.aux *.toc *.bbl *.blg
	rm -f $(MAIN_TEX).pdf

.SUFFIXES: .pdf .tex

$(MAIN_TEX).pdf:
	$(TEX) $(TEX_ARG) $(MAIN_TEX).tex
