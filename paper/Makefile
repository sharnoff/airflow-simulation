main: figs/sparse/*.tikz.tex
	latexmk -pdf -interaction=nonstopmode main.tex

figs/sparse/%.tikz.tex: figs/sparse/%.build.json
	figs/sparse/generate.py $<

clean:
	latexmk -c
	rm figs/sparse/*.tikz.tex
