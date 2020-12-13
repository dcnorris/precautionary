.PHONY: scryer sicstus scryer_tabs sicstus_tabs sysdata

scryer: scryer_tabs sysdata

sicstus: sicstus_tabs sysdata

sysdata : T2.tab T3.tab T4.tab T5.tab T6.tab T7.tab T8.tab make_sysdata_TUb.R
	Rscript make_sysdata_TUb.R

scryer_tabs : ../exec/scryer.pl ../exec/esc.pl
	scryer-prolog ../exec/scryer.pl ../exec/esc.pl

sicstus_tabs : ../exec/sicstus.pl ../exec/esc.pl
	sicstus -l ../exec/sicstus.pl -l ../exec/esc.pl
