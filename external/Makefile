.PHONY: all nulibclean NuLib EOSsuper EOSsuperclean

all:

#########
# NuLib #
#########
NuLib:
	touch NuLib/make.inc # trick nulib into thinking make.inc exists there
	$(MAKE) -j1 -C NuLib

nulibclean:
	$(MAKE) -C NuLib clean
	rm -f NuLib/make.inc

############
# EOSsuper #
############
EOSsuper:
	$(MAKE) -C EOSsuper

EOSsuperclean:
	$(MAKE) -C EOSsuper clean
