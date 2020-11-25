all: figures tables

figures: fig1.csv fig2.csv fig3hist.csv fig3invlap.csv 
    # mhtellherr.csv mhteproberr.csv mchist.csv mcinvlap.csv 

fig1.csv: figure1.m strkdur.asc
	matlab -batch figure1

fig2.csv: figure2.m
	matlab -batch figure2

fig3hist.csv fig3invlap.csv: figure3.m 
	matlab -batch figure3 > figure3.out

tables: tab1.tex
	$(here come the tables)

clean:
	rm -f *.csv *.out

git:
	git add .
	git commit -m"Lazy add, commit, and push (via make)"
	git push

.PHONY: clean git
