# Makefile for The Likelihood of Mixed Hitting Times: Replication Package
# Type 'make help' for instructions

.PHONY: all help figures tables other clean git

# variables
rfile = replication
texopt = -interaction=batchmode
figobjects = fig1.csv fig1times.tex fig2.csv fig3hist.csv fig3invlap.csv fig4.csv fig4.tex weibullmph.mat 
# mhtellherr.csv mhteproberr.csv mchist.csv mcinvlap.csv mhthazard.txt
tabobjects = tab1.tex tab1times.tex tab1.mat
othobjects = chckgrad.tex
texobjects = $(rfile).aux $(rfile).brf $(rfile).log $(rfile).out $(rfile).synctex.gz

# main rules
all: paper other $(rfile).pdf
paper: figures tables

help:
	@echo "----------------------------------------------------------------"
	@echo "Replication of 'The Likelihood of Mixed Hitting Times'          "
	@echo "  Jaap H. Abbring and Tim Salimans                              "
	@echo "----------------------------------------------------------------"
	@echo ""
	@echo " make all   --  replicate all tables & figures and display them "
	@echo "                  in replication.pdf (default)                  "
	@echo " make paper --  only generate results for paper
	@echo " make clean --  remove all output files                         "
	@echo ""

# replication figures
figures: $(figobjects)

fig1.csv fig1times.tex: figure1.m strkdur.asc migaussmle.m lhmigauss.m numinvlap.m
	matlab -batch figure1

fig2.csv: figure2.m numinvlap2.m
	matlab -batch figure2

fig3hist.csv fig3invlap.csv: figure3.m simmht.m numinvlap.m 
	matlab -batch figure3

fig4.csv fig4.tex weibullmph.mat: figure4.m tab1.mat igausscdf.m igausspdf.m weibullcdf.m weibullpdf.m
	matlab -batch figure4

# replication tables
tables: $(tabobjects)

tab1.tex tab1times tab1.mat: table1.m strkdur.asc mhtmle.m
	matlab -batch table1

# run pdfLaTeX to create output file
replication.pdf: $(figobjects) $(tabobjects) $(othobjects) replication.tex
	pdflatex $(texopt) $(rfile).tex
	pdflatex $(texopt) $(rfile).tex
	open $(rfile).pdf
	rm -f $(texobjects)

# other checks 
other: $(othobjects)

chckgrad.tex: checkgradient.m strkdur.asc weibullmph.mat numgrad.m mhtobj.m nllhmph.m 
	matlab -batch checkgradient

# remove output
clean:
	rm -f $(figobjects) $(tabobjects) $(othobjects) $(texobjects)

# Lazy one-command add, commit, and push to Github
git:
	git add .
	git commit -m"Lazy add, commit, and push (via make)"
	git push
