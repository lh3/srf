.SUFFIXES: .gp .tex .eps .pdf .eps.gz

.eps.pdf:
		epstopdf --outfile $@ $<

.eps.gz.pdf:
		gzip -dc $< | epstopdf --filter > $@

.pdf.eps:
		pdftops -eps $< $@

all:srf.pdf

srf.pdf:srf.tex srf.bib
		pdflatex srf; bibtex srf; pdflatex srf; pdflatex srf;

clean:
		rm -fr *.toc *.aux *.bbl *.blg *.idx *.log *.out *~ srf.pdf
