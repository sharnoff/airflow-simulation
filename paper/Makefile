SPARSEFIGS := $(shell find figs/sparse -name '*.build.json' | sed -e 's/\.build\.json$$/.png/g')
AUTOPLOTS := $(shell find figs -name '*.autoplot.tex' | sed -e 's/\.autoplot\.tex$$/.generated.tex/g')

LATEXMK_ARGS := -f -pdf -interaction=nonstopmode

main.pdf: $(SPARSEFIGS) $(AUTOPLOTS) figs/*.tex *.tex ref.bib
	latexmk $(LATEXMK_ARGS) main.tex

figs/sparse/%.png: figs/sparse/%.build.json figs/sparse/generate.py
	figs/sparse/generate.py $<

figs/%.generated.tex: figs/%.autoplot.tex figs/mkplot.py data
	figs/mkplot.py $<

data: ../data.tar.zst
	tar -C .. -xf ../data.tar.zst
	ln -s ../data data

clean:
	latexmk -c
	rm -f main.pdf
	rm -f $(SPARSEFIGS)
	rm -f figs/sparse/*.png
	rm -f figs/*.generated.tex
