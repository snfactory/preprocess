PDF = UserManual.pdf

all : $(PDF)

.PHONY : clean

%.pdf : %.dvi
	dvipdf $<

%.dvi : %.tex
	latex $<
	latex $<

clear: 
	rm -f *.dvi *.log *.toc *.aux *.pdf

UserManual.dvi : TitlePage.tex Executables.tex
