.PHONY: scryer sicstus scryer_tabs sicstus_tabs sysdata

scryer: scryer_tabs sysdata

sicstus: sicstus_tabs sysdata

sysdata : T2.tab T3.tab T4.tab T5.tab T6.tab T7.tab T8.tab make_sysdata_TUb.R
	Rscript make_sysdata_TUb.R

scryer_tabs : ../exec/prolog/scryer.pl ../exec/prolog/esc.pl
	scryer-prolog ../exec/prolog/scryer.pl ../exec/prolog/esc.pl

sicstus_tabs : ../exec/prolog/sicstus.pl ../exec/prolog/esc.pl
	sicstus -l ../exec/prolog/sicstus.pl -l ../exec/prolog/esc.pl
