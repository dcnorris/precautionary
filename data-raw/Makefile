.PHONY: sysdata tabfiles

sysdata : tabfiles make_sysdata_AUb.R
	Rscript make_sysdata_AUb.R

tabfiles : ../exec/prolog/esc.pl
	swipl ../exec/prolog/esc.pl
#	scryer-prolog ../exec/prolog/esc.pl
