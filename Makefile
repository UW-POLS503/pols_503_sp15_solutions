RSCRIPT = Rscript
RMD_FILES = $(wildcard *.Rmd)
HTML_FILES = $(RMD_FILES:%.Rmd=%.html)
MD_FILES = $(RMD_FILES:%.Rmd=%.md)
PDF_FILES = $(RMD_FILES:%.Rmd=%.pdf)

all: build

build: html pdf

html: $(HTML_FILES)

pdf: $(PDF_FILES)

%.html: %.Rmd
	$(RSCRIPT) -e 'rmarkdown::render("$^", output_file="$@", runtime="static", quiet=TRUE, output_format="html_document")'

%.pdf: %.Rmd
	$(RSCRIPT) -e 'rmarkdown::render("$^", output_file="$@", runtime="static", quiet=TRUE, output_format="pdf_document")'

clean:
	-rm $(HTML_FILES) $(PDF_FILES) $(MD_FILES) *_cache *_files

.PHONY: all build html pdf
