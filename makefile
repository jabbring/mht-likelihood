# Makefile for The Likelihood of Mixed Hitting Times: Replication Package
# Type 'make help' for instructions

.PHONY: all help figures tables other clean git

# variables
rfile = replication
texopt = -interaction=batchmode
datafile = strkdur.asc
figobjects = fig1.csv fig1times.tex fig2.csv fig3hist.csv fig3invlap.csv fig4.csv fig4.tex weibullmph.mat 
tabobjects = tab1.tex tab1times.tex tab1.mat tab1lowM.tex tab1lowMtimes.tex tab1mig.tex tab1migtimes.tex
othobjects = chckgrad.tex
texobjects = $(rfile).aux $(rfile).brf $(rfile).log $(rfile).out $(rfile).synctex.gz

specs = gammagamma.m gammapoint.m pointgamma.m pointpoint.m
mhtobj = mhtobj.m numinvlap.m
mhtgrad = mhtgrad.m numinvlap.m
mhtmle = mhtmle.m $(specs) $(mhtobj) $(mhtgrad) numinvlap.m numjac.m
mhtobj2 = mhtobj2.m numinvlap2.m
mhtgrad2 = mhtgrad2.m numinvlap2.m
mhtmle2 = mhtmle2.m $(specs) $(mhtobj2) $(mhtgrad2) numinvlap2.m numjac.m

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

fig1.csv fig1times.tex: figure1.m $(datafile) $(migaussmle) $(lhmigauss) numinvlap.m numinvlap2.m pointpoint.m
	matlab -batch figure1

fig2.csv: figure2.m numinvlap2.m pointpoint.m
	matlab -batch figure2

fig3hist.csv fig3invlap.csv: figure3.m $(simmht) numinvlap.m $(specs)
	matlab -batch figure3

fig4.csv fig4.tex weibullmph.mat: figure4.m $(datafile) tab1.mat $(igspecs) $(mphspecs) $(nllhmph) 
	matlab -batch figure4

# replication tables
tables: $(tabobjects)

tab1.tex tab1times tab1.mat: table1.m $(datafile) $(mhtmle)
	matlab -batch table1

tab1lowM.tex tab1lowMtimes: table1lowM.m $(datafile) $(mhtmle2)
	matlab -batch table1lowM

tab1mig.tex tab1migtimes: table1mig.m $(datafile) $(migaussmle)
	matlab -batch table1mig

# run pdfLaTeX to create output file
replication.pdf: $(figobjects) $(tabobjects) $(othobjects) replication.tex
	pdflatex $(texopt) $(rfile).tex
	pdflatex $(texopt) $(rfile).tex
	open $(rfile).pdf
	rm -f $(texobjects)

# other checks 
other: $(othobjects)

chckgrad.tex: checkgradient.m $(datafile) numgrad.m $(mhtobj) $(specs) weibullmph.mat $(nllhmph)
	matlab -batch checkgradient

# remove output
clean:
	rm -f $(figobjects) $(tabobjects) $(othobjects) $(texobjects)

# Lazy one-command add, commit, and push to Github
git:
	git add .
	git commit -m"Lazy add, commit, and push (via make)"
	git push
