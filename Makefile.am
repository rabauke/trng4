AUTOMAKE_OPTIONS = gnu
SUBDIRS = trng src examples
EXTRA_DIST = doc/trng.pdf

doc/trng.pdf:
	cd doc; pdflatex trng; bibtex trng; pdflatex trng; pdflatex trng
