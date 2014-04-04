AUTOMAKE_OPTIONS = gnu
ACLOCAL_AMFLAGS = -I m4
SUBDIRS = trng src examples
DIST_SUBDIRS = trng src examples
EXTRA_DIST = doc/trng.pdf

doc/trng.pdf:
	cd doc; make trng.pdf

examples: all
	cd examples && $(MAKE) examples

.PHONY: examples
