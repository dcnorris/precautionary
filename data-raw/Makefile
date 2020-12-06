.PHONY: swi scryer swi_tabs scryer_tabs sysdata

swi: swi_tabs sysdata

scryer: scryer_tabs sysdata

sysdata : T2.tab T3.tab T4.tab T5.tab T6.tab T7.tab T8.tab make_sysdata_AUb.R
	Rscript make_sysdata_AUb.R

swi_tabs : ../exec/prolog/swi.pl ../exec/prolog/esc.pl
	swipl ../exec/prolog/swi.pl ../exec/prolog/esc.pl

scryer_tabs : ../exec/prolog/scryer.pl ../exec/prolog/esc.pl
	scryer-prolog ../exec/prolog/scryer.pl ../exec/prolog/esc.pl
