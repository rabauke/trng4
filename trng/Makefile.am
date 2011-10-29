AM_CXXFLAGS = @BOOST_CXXFLAGS@ @CXXFLAGS@

header-links: remove-links
	HEADERLIST="$(top_srcdir)/src/*.hpp"; \
	for h in $$HEADERLIST; do \
	  BASENAME=`basename $$h`; \
	  test -r $$BASENAME || $(LN_S) $$h $$BASENAME; \
	done

remove-links: 
	rm -f *.hpp

all: all-am header-links

clean-local: clean-am remove-links

distclean-local: distclean-am remove-links
	-rm -f Makefile
