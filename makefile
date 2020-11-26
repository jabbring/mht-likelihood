# Makefile for The Likelihood of Mixed Hitting Times: Replication Package
# Type 'make help' for instructions

.PHONY: all help figures clean git

# variables
rfile = replication
texopt = -interaction=batchmode
figobjects = fig1.csv fig2.csv fig3hist.csv fig3invlap.csv # mhtellherr.csv mhteproberr.csv mchist.csv mcinvlap.csv
texobjects = $(rfile).aux $(rfile).brf $(rfile).log $(rfile).out

# main rules
all: figures $(rfile).pdf

help:
	@echo "----------------------------------------------------------------"
	@echo "Replicate 'The Likelihood of Mixed Hitting Times'               "
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
	matlab -batch figure3 > figure3.out

# replication tables
#tables: tab1.tex
#	$(here come the tables)

# run pdfLaTeX to create output file
replication.pdf: $(figobjects) replication.tex
	pdflatex $(texopt) $(rfile).tex
	pdflatex $(texopt) $(rfile).tex
	open $(rfile).pdf
	rm -f $(texobjects)

# remove output
clean:
	rm -f $(figobjects) figure3.out $(texobjects) $(rfile).pdf

# Lazy one-command add, commit, and push to Github
git:
	git add .
	git commit -m"Lazy add, commit, and push (via make)"
	git push