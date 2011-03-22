PDFFILES = 
TEXFILES = 

SUBDIRS =

all: all-recursive dune-fem-functionals.bbl dune-fem-functionals.pdf

clean: clean-recursive
	@rm -f *.pdf
	@rm -f *.aux *.bbl *.blg *.log

PRECIOUS: dune-fem-functionals.aux

%.aux: %.tex $(TEXFILES) $(PDFFILES) Makefile
	@echo "Running latex...";
	@pdflatex --interaction=nonstopmode $< > /dev/null

%.bbl: %.aux %.bib Makefile
	@echo "Running bibtex...";
	@bibtex $< > /dev/null 2>/dev/null

%.pdf: %.tex %.bbl $(TEXFILES) $(PDFFILES) Makefile
	@echo "Rerun" > $*.log;
	@while grep "Rerun" $*.log > /dev/null ; do \
	  echo "Running pdflatex..." ; \
	  pdflatex --interaction=nonstopmode $< > /dev/null 2> /dev/null; \
	done

all-recursive:
	@for DIR in $(SUBDIRS) ; do (cd $$DIR ; $(MAKE) -s all) ; done

clean-recursive:
	@for DIR in $(SUBDIRS) ; do (cd $$DIR ; $(MAKE) -s clean) ; done
