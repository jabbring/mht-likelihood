rfile = replication

texopt = -interaction=batchmode

figobjects = fig1.csv fig2.csv fig3hist.csv fig3invlap.csv 
    # mhtellherr.csv mhteproberr.csv mchist.csv mcinvlap.csv

texobjects = $(rfile).aux $(rfile).brf $(rfile).log $(rfile).out

all: figures $(rfile).pdf

figures: $(figobjects)

fig1.csv: figure1.m strkdur.asc
	matlab -batch figure1

fig2.csv: figure2.m
	matlab -batch figure2

fig3hist.csv fig3invlap.csv: figure3.m 
	matlab -batch figure3 > figure3.out

#tables: tab1.tex
#	$(here come the tables)

replication.pdf: $(figobjects) replication.tex
	pdflatex $(texopt) $(rfile).tex
	pdflatex $(texopt) $(rfile).tex
	open $(rfile).pdf
	rm -f $(texobjects)

clean:
	rm -f $(figobjects) figure3.out $(texobjects) $(rfile).pdf

git:
	git add .
	git commit -m"Lazy add, commit, and push (via make)"
	git push
	
.PHONY: figures clean git
