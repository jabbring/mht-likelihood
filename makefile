# Makefile for The Likelihood of Mixed Hitting Times: Replication Package
# Type 'make help' for instructions

.PHONY: all help figures tables other clean git

# variables
rfile = replication
texopt = -interaction=batchmode
figobjects = fig1.csv fig1times.tex fig2.csv fig3hist.csv fig3invlap.csv fig4.csv fig4.tex weibullmph.mat 
tabobjects = tab1.tex tab1times.tex tab1.mat tab1mig.tex tab1migtimes.tex tab1lowM.tex tab1lowMtimes.tex 
othobjects = chckgrad.tex
texobjects = $(rfile).aux $(rfile).brf $(rfile).log $(rfile).out $(rfile).synctex.gz

specs = gammagamma.m gammapoint.m pointgamma.m pointpoint.m
numinvlap = numinvlap.m $(specs)
mhtobj = mhtobj.m $(numinvlap)
mhtgrad = mhtgrad.m $(numinvlap)
mhtmle = $(mhtobj) $(mhtgrad) $(numinvlap) numjac.m
numinvlap2 = numinvlap2.m $(specs)
mhtobj2 = mhtobj2.m $(numinvlap2)
mhtgrad2 = mhtgrad2.m $(numinvlap2)
mhtmle = $(mhtobj2) $(mhtgrad2) $(numinvlap2) numjac.m

igspecs = igausscdf.m igausspdf.m
nllhmigauss = nllhmigauss.m $(igspecs)
lhmigauss = lhmigauss.m $(igspecs)
migaussmle = migaussmle.m $(nllhmigauss)

mphspecs = weibullcdf.m weibullpdf.m
nllhmph = nllhmph.m $(mphspecs)

simmht = simmht.m randraw.m

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

fig1.csv fig1times.tex: figure1.m strkdur.asc $(migaussmle) $(lhmigauss) $(numinvlap) $(numinvlap2)
	matlab -batch figure1

fig2.csv: figure2.m $(numinvlap2)
	matlab -batch figure2

fig3hist.csv fig3invlap.csv: figure3.m $(simmht) $(numinvlap) 
	matlab -batch figure3

fig4.csv fig4.tex weibullmph.mat: figure4.m tab1.mat $(igspecs) $(mphspecs) $(nllhmph) 
	matlab -batch figure4

# replication tables
tables: $(tabobjects)

tab1.tex tab1times tab1.mat: table1.m strkdur.asc $(mhtmle)
	matlab -batch table1

tab1mig.tex tab1migtimes: table1mig.m strkdur.asc $(migaussmle)
	matlab -batch table1mig

tab1lowM.tex tab1lowMtimes: table1lowM.m strkdur.asc $(mhtmle2)
	matlab -batch table1lowM

# run pdfLaTeX to create output file
replication.pdf: $(figobjects) $(tabobjects) $(othobjects) replication.tex
	pdflatex $(texopt) $(rfile).tex
	pdflatex $(texopt) $(rfile).tex
	open $(rfile).pdf
	rm -f $(texobjects)

# other checks 
other: $(othobjects)

chckgrad.tex: checkgradient.m strkdur.asc weibullmph.mat numgrad.m $(mhtobj) $(nllhmph) 
	matlab -batch checkgradient

# remove output
clean:
	rm -f $(figobjects) $(tabobjects) $(othobjects) $(texobjects)

# Lazy one-command add, commit, and push to Github
git:
	git add .
	git commit -m"Lazy add, commit, and push (via make)"
	git push
