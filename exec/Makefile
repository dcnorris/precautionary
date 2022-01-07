.PHONY: scryer sicstus scryer_tabs sicstus_tabs sysdata

scryer: scryer_tabs sysdata

sicstus: sicstus_tabs sysdata

sysdata : T2.tab T3.tab T4.tab T5.tab T6.tab T7.tab T8.tab make_sysdata_TUb.R
	Rscript make_sysdata_TUb.R

scryer_tabs : scryer.pl esc.pl
	scryer-prolog scryer.pl esc.pl

sicstus_tabs : sicstus.pl esc.pl
	sicstus -l sicstus.pl -l esc.pl

benchmark:
	scryer-prolog --version
	time scryer-prolog scryer.pl esc.pl
