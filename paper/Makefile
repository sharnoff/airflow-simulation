SPARSEFIGS := $(shell find figs/sparse -name '*.build.json' | sed -e 's/\.build\.json$$/.png/g')

LATEXMK_ARGS := -f -pdf -interaction=nonstopmode

main.pdf: $(SPARSEFIGS) $(TIKZPDFS) *.tex
	latexmk $(LATEXMK_ARGS) main.tex

figs/sparse/%.png: figs/sparse/%.build.json figs/sparse/generate.py
	figs/sparse/generate.py $<

clean:
	latexmk -c
	rm -f main.pdf
	rm -f $(TIKZFIGS) $(SPARSEFIGS)
	rm -f figs/sparse/*.png
