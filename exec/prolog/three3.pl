% Attempt an enumeration of ALL 3+3 trials
:- use_module(library(clpfd)).

% Type constraints
istox(T) :- T in 0..3, indomain(T).
isdose(D) :- D in 0..3, indomain(D).

% Successive dose relation with constraints
dose_up1(D, D1) :- isdose(D), D1 #= D + 1, isdose(D1).

% Let esc(D) *describe* a list of 3+3 cohorts starting from dose D
esc(D) --> [D ^ T], { T in 0..3, indomain(T) },
	   (   { T #= 0 } -> (   { dose_up1(D,D1) } -> esc(D1)
			     ;	 [mtd_notfound(D)]
			     )
	   ;   { T #= 1 } -> sta(D)
	   ;   []  % TODO: Is such explicit termination needed?
	   ).
sta(D) --> [D ^ T], { T in 0..3, indomain(T) },
	   (   { T #= 0 } -> (   { dose_up1(D,D1) } -> esc(D1)
			     ;	 [declare_mtd(D)]
			     )
	   ;   { T #= 1 } -> [declare_mtd(D)]
	   ;   { dose_up1(D_1,D) } -> [declare_mtd(D_1)]
	   ).
