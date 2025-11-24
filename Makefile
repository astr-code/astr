# ====================================
# Top-level Makefile
# ====================================

.PHONY: all astr pastr clean

all: astr pastr

astr:
	$(MAKE) -f Makefile.astr

pastr:
	$(MAKE) -f Makefile.pastr

astr_clean:
	$(MAKE) -f Makefile.astr clean

pastr_clean:
	$(MAKE) -f Makefile.pastr clean

clean:
	$(MAKE) -f Makefile.astr clean
	$(MAKE) -f Makefile.pastr clean