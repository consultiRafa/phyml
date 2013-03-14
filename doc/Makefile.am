docdir = $(datadir)/doc

docfiles = phyml-manual.tex

if HAVE_PDFLATEX
docfiles += phyml-manual.pdf
endif

MANNAME = phyml-manual
MANTEXSRC = $(MANNAME).tex

CLEANFILES = $(MANNAME).pdf $(MANNAME).log $(MANNAME).idx $(MANNAME).out \
 $(MANNAME).toc $(MANNAME).aux

phyml-manual.pdf: $(MANTEXSRC)
	pdflatex $< > /dev/null


all-am: intro phyml-manual.pdf
	@echo ""
	@echo "Done."

intro: 
	@echo ""
	@echo ""
	@echo ".: Building documentation :."
	@echo ""
	@echo ""

dist_doc_DATA = $(docfiles)
