all: linearfit.pdf linearfit_notes_andrew.pdf

%.pdf: %.ps
	ps2pdf $^ $@

%.gz: %
	gzip $^

%.ps: %.dvi
	dvips -Ppdf -o $@ $^

%.dvi: %.tex
	latex $^
	latex $^

%.tar: %.tex $(wildcard %.eps)
	tar -zcf $@ $^

