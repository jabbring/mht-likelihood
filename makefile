# Makefile for The Likelihood of Mixed Hitting Times: Replication Package
# Type 'make help' for instructions

.PHONY: all help figures tables clean git

# variables
rfile = replication
texopt = -interaction=batchmode
figobjects = fig1.csv fig2.csv fig3hist.csv fig3invlap.csv fig4.csv # mhtellherr.csv mhteproberr.csv mchist.csv mcinvlap.csv mhthazard.txt
tabobjects = tab1.tex
texobjects = $(rfile).aux $(rfile).brf $(rfile).log $(rfile).out

# main rules
all: figures tables $(rfile).pdf

help:
	@echo "----------------------------------------------------------------"
	@echo "Replication of 'The Likelihood of Mixed Hitting Times'          "
	@echo "  Jaap H. Abbring and Tim Salimans                              "
	@echo "----------------------------------------------------------------"
	@echo ""
	@echo " make all   --  replicate all tables & figures and display them "
	@echo "                  in replication.pdf (default)                  "
	@echo " make clean --  remove all output files                         "
	@echo ""

# replication figures
figures: $(figobjects)

fig1.csv: figure1.m strkdur.asc
	matlab -batch figure1

fig2.csv: figure2.m
	matlab -batch figure2

fig3hist.csv fig3invlap.csv: figure3.m 
	matlab -batch figure3

fig4.csv: figure4.m fig4.mat

# replication tables
tables: $(tabobjects)

tab1.tex fig4.mat: table1.m strkdur.asc
	matlab -batch table1

# run pdfLaTeX to create output file
replication.pdf: $(figobjects) $(tabobjects) replication.tex
	pdflatex $(texopt) $(rfile).tex
	pdflatex $(texopt) $(rfile).tex
	open $(rfile).pdf
	rm -f $(texobjects)

# remove output
clean:
	rm -f $(figobjects) fig4.mat $(tabobjects) $(texobjects) $(rfile).pdf

# Lazy one-command add, commit, and push to Github
git:
	git add .
	git commit -m"Lazy add, commit, and push (via make)"
	git push
