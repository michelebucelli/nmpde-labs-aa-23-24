.PHONY: all clean distclean

all: lab-09.pdf lab-09-solution.pdf

lab-09.pdf: lab-09.tex lab-09-common.tex ../../common/preamble.tex
	latexmk -pdf $<

lab-09-solution.pdf: lab-09-solution.tex lab-09-common.tex ../../common/preamble.tex ../../common/preamble-solution.tex
	latexmk -pdf $<

clean:
	rm -f *.aux *.fls *.tdo *.blg *.bbl *.out *.log *.fdb_latexmk

distclean: clean
	rm *.pdf
