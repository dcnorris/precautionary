% Enumerate 3+3 trial outcomes for exact matrix calculations
:- use_module(library(clpfd)).
:- use_module(library(pio)).

% Prefix op * for 'generalizing away' goals (https://www.metalevel.at/prolog/debugging)
:- op(920, fy, *). *_.  % *Goal always succeeds

% DCG esc(D, Lo..Hi) describes a list of 3+3 cohorts FOLLOWING a dose D in Lo..Hi.
% One can read esc(D, Lo..Hi) as the DECISION to escalate from D to min(D+1,Hi).
% For example, esc(0, 0..5) *initiates* a trial that has 5 prespecified doses,
% and enrolls the first cohort at D=1.

tox(T) :- T in 0..3,
	  indomain(T).

% Mnemonic: * is ^ that 'splatted' on dose ceiling.
esc(Hi, Lo..Hi) --> { Lo #< Hi },
		    [Hi * T], { tox(T) },
		    (  {T #=< 1}, [mtd_notfound(Hi)]
		    ;  {T #>= 2}, des(Hi, Lo)
		    ).
esc(D, Lo..Hi) --> { D1 #= D + 1, D1 in Lo..Hi },
		   [D1 ^ T], { tox(T) },
		   (  {T #= 0}, esc(D1, Lo..Hi)
		   ;  {T #= 1}, sta(D1, Lo..Hi)
		   ;  {T #> 1}, des(D1, Lo)
		   ).

sta(D, _..D) --> [D - 0], [mtd_notfound(D)].
sta(D, Lo..Hi) --> { D #< Hi, D in Lo..Hi },
		   [D - 0],
		   esc(D, D..Hi).
sta(D, Lo.._) --> [D - T], { tox(T), T #> 0 },
		  des(D, Lo).

% As a mirror image of esc//2, des(D, Lo) moves
% downward FROM D, to max(D-1,Lo).
% NB: De-escalation to D-1 clamps Hi #= D - 1.
%% TODO: Does this mean des//2 could be written as des(Lo..D),
%%       somehow exploiting the impossibility of Y..X with Y > X?
des(D, Lo) --> { D_1 #= D - 1 },
	       (  {D_1 #= Lo}, [declare_mtd(Lo)]
	       ;  {D_1 #> Lo}, [D_1 : T], {tox(T)},
		  (  {T #=< 1}, [declare_mtd(D_1)]
		  ;  {T #>= 2}, des(D_1, Lo)
		  )
	       ).

% n_trials(+Drange, DN)
n_trials(Drange, DN) :-
    Dmax in Drange, indomain(Dmax),
    findall(Tr, phrase(esc(0, 0..Dmax), Tr), Trials),
    length(Trials, N),
    DN = (Dmax, N).

%?- n_trials(1..10, D_N).
%@ D_N =  (1, 10) ;
%@ D_N =  (2, 46) ;
%@ D_N =  (3, 154) ;
%@ D_N =  (4, 442) ;
%@ D_N =  (5, 1162) ;
%@ D_N =  (6, 2890) ;
%@ D_N =  (7, 6922) ;
%@ D_N =  (8, 16138) ;
%@ D_N =  (9, 36874) ;
%@ D_N =  (10, 82954).

?- phrase(esc(0, 0..2), Tr).
%@ Tr = [1^0, 2^0, 2*0, mtd_notfound(2)] ;
%@ Tr = [1^0, 2^0, 2*1, mtd_notfound(2)] ;
%@ Tr = [1^0, 2^0, 2*2, 1:0, declare_mtd(1)] ;
%@ Tr = [1^0, 2^0, 2*2, 1:1, declare_mtd(1)] ;
%@ Tr = [1^0, 2^0, 2*2, 1:2, declare_mtd(0)] ;
%@ Tr = [1^0, 2^0, 2*2, 1:3, declare_mtd(0)] ;
%@ Tr = [1^0, 2^0, 2*3, 1:0, declare_mtd(1)] ;
%@ Tr = [1^0, 2^0, 2*3, 1:1, declare_mtd(1)] ;
%@ Tr = [1^0, 2^0, 2*3, 1:2, declare_mtd(0)] ;
%@ Tr = [1^0, 2^0, 2*3, 1:3, declare_mtd(0)] ;
%@ Tr = [1^0, 2^1, 2-0, mtd_notfound(2)] ;
%@ Tr = [1^0, 2^1, 2-1, 1:0, declare_mtd(1)] ;
%@ Tr = [1^0, 2^1, 2-1, 1:1, declare_mtd(1)] ;
%@ Tr = [1^0, 2^1, 2-1, 1:2, declare_mtd(0)] ;
%@ Tr = [1^0, 2^1, 2-1, 1:3, declare_mtd(0)] ;
%@ Tr = [1^0, 2^1, 2-2, 1:0, declare_mtd(1)] ;
%@ Tr = [1^0, 2^1, 2-2, 1:1, declare_mtd(1)] ;
%@ Tr = [1^0, 2^1, 2-2, 1:2, declare_mtd(0)] ;
%@ Tr = [1^0, 2^1, 2-2, 1:3, declare_mtd(0)] ;
%@ Tr = [1^0, 2^1, 2-3, 1:0, declare_mtd(1)] ;
%@ Tr = [1^0, 2^1, 2-3, 1:1, declare_mtd(1)] ;
%@ Tr = [1^0, 2^1, 2-3, 1:2, declare_mtd(0)] ;
%@ Tr = [1^0, 2^1, 2-3, 1:3, declare_mtd(0)] ;
%@ Tr = [1^0, 2^2, 1:0, declare_mtd(1)] ;
%@ Tr = [1^0, 2^2, 1:1, declare_mtd(1)] ;
%@ Tr = [1^0, 2^2, 1:2, declare_mtd(0)] ;
%@ Tr = [1^0, 2^2, 1:3, declare_mtd(0)] ;
%@ Tr = [1^0, 2^3, 1:0, declare_mtd(1)] ;
%@ Tr = [1^0, 2^3, 1:1, declare_mtd(1)] ;
%@ Tr = [1^0, 2^3, 1:2, declare_mtd(0)] ;
%@ Tr = [1^0, 2^3, 1:3, declare_mtd(0)] ;
%@ Tr = [1^1, 1-0, 2^0, 2*0, mtd_notfound(2)] ;
%@ Tr = [1^1, 1-0, 2^0, 2*1, mtd_notfound(2)] ;
%@ Tr = [1^1, 1-0, 2^0, 2*2, declare_mtd(1)] ;
%@ Tr = [1^1, 1-0, 2^0, 2*3, declare_mtd(1)] ;
%@ Tr = [1^1, 1-0, 2^1, 2-0, mtd_notfound(2)] ;
%@ Tr = [1^1, 1-0, 2^1, 2-1, declare_mtd(1)] ;
%@ Tr = [1^1, 1-0, 2^1, 2-2, declare_mtd(1)] ;
%@ Tr = [1^1, 1-0, 2^1, 2-3, declare_mtd(1)] ;
%@ Tr = [1^1, 1-0, 2^2, declare_mtd(1)] ;
%@ Tr = [1^1, 1-0, 2^3, declare_mtd(1)] ;
%@ Tr = [1^1, 1-1, declare_mtd(0)] ;
%@ Tr = [1^1, 1-2, declare_mtd(0)] ;
%@ Tr = [1^1, 1-3, declare_mtd(0)] ;
%@ Tr = [1^2, declare_mtd(0)] ;
%@ Tr = [1^3, declare_mtd(0)] ;
%@ false.
%@ Tr = [1^0, 2^0, 2*0, mtd_notfound(2)] ;
%@ Tr = [1^0, 2^0, 2*1, mtd_notfound(2)] ;
%@ Tr = [1^0, 2^0, 2*2, 1:0, declare_mtd(1)] ;
%@ Tr = [1^0, 2^0, 2*2, 1:1, declare_mtd(1)] ;
%@ Tr = [1^0, 2^0, 2*2, 1:2, declare_mtd(0)] ;
%@ Tr = [1^0, 2^0, 2*2, 1:3, declare_mtd(0)] ;
%@ Tr = [1^0, 2^0, 2*3, 1:0, declare_mtd(1)] ;
%@ Tr = [1^0, 2^0, 2*3, 1:1, declare_mtd(1)] ;
%@ Tr = [1^0, 2^0, 2*3, 1:2, declare_mtd(0)] ;
%@ Tr = [1^0, 2^0, 2*3, 1:3, declare_mtd(0)] ;
%@ Tr = [1^0, 2^1, 2-0, mtd_notfound(2)] ;
%@ Tr = [1^0, 2^1, 2-1, 1:0, declare_mtd(1)] ;
%@ Tr = [1^0, 2^1, 2-1, 1:1, declare_mtd(1)] ;
%@ Tr = [1^0, 2^1, 2-1, 1:2, declare_mtd(0)] ;
%@ Tr = [1^0, 2^1, 2-1, 1:3, declare_mtd(0)] ;
%@ Tr = [1^0, 2^1, 2-2, 1:0, declare_mtd(1)] ;
%@ Tr = [1^0, 2^1, 2-2, 1:1, declare_mtd(1)] ;
%@ Tr = [1^0, 2^1, 2-2, 1:2, declare_mtd(0)] ;
%@ Tr = [1^0, 2^1, 2-2, 1:3, declare_mtd(0)] ;
%@ Tr = [1^0, 2^1, 2-3, 1:0, declare_mtd(1)] ;
%@ Tr = [1^0, 2^1, 2-3, 1:1, declare_mtd(1)] ;
%@ Tr = [1^0, 2^1, 2-3, 1:2, declare_mtd(0)] ;
%@ Tr = [1^0, 2^1, 2-3, 1:3, declare_mtd(0)] ;
%@ Tr = [1^0, 2^2, 1:0, declare_mtd(1)] ;
%@ Tr = [1^0, 2^2, 1:1, declare_mtd(1)] ;
%@ Tr = [1^0, 2^2, 1:2, declare_mtd(0)] ;
%@ Tr = [1^0, 2^2, 1:3, declare_mtd(0)] ;
%@ Tr = [1^0, 2^3, 1:0, declare_mtd(1)] ;
%@ Tr = [1^0, 2^3, 1:1, declare_mtd(1)] ;
%@ Tr = [1^0, 2^3, 1:2, declare_mtd(0)] ;
%@ Tr = [1^0, 2^3, 1:3, declare_mtd(0)] ;
%@ Tr = [1^1, 1-0, 2^0, 2*0, mtd_notfound(2)] ;
%@ Tr = [1^1, 1-0, 2^0, 2*1, mtd_notfound(2)] ;
%@ Tr = [1^1, 1-0, 2^0, 2*2, declare_mtd(1)] ;
%@ Tr = [1^1, 1-0, 2^0, 2*3, declare_mtd(1)] ;
%@ Tr = [1^1, 1-0, 2^1, 2-0, mtd_notfound(2)] ;
%@ Tr = [1^1, 1-0, 2^1, 2-1, declare_mtd(1)] ;
%@ Tr = [1^1, 1-0, 2^1, 2-2, declare_mtd(1)] ;
%@ Tr = [1^1, 1-0, 2^1, 2-3, declare_mtd(1)] ;
%@ Tr = [1^1, 1-0, 2^2, declare_mtd(1)] ;
%@ Tr = [1^1, 1-0, 2^3, declare_mtd(1)] ;
%@ Tr = [1^1, 1-1, declare_mtd(0)] ;
%@ Tr = [1^1, 1-2, declare_mtd(0)] ;
%@ Tr = [1^1, 1-3, declare_mtd(0)] ;
%@ Tr = [1^2, declare_mtd(0)] ;
%@ Tr = [1^3, declare_mtd(0)] ;
%@ false.

%% Transform dose-escalation path lists
%% to the arrays T(c,d,j).

path_matrix(P, D, M) :-
    phrase(pm_(M), P),
    M = C1-C2,
    length(C1,D),
    length(C2,D),
    maplist(ground_or_nil, C1),
    maplist(ground_or_nil, C2).

pm_(C1-C2) -->
    [D^T], % ascribe to D's 1st cohort
    { nth1(D, C1, T) },
    pm_(C1-C2).

pm_(C1-C2) -->
    (	[D-T]
    ;	[D*T]
    ;	[D:T]
    ), % ascribe to D's 2nd cohort
    { nth1(D, C2, T) },
    pm_(C1-C2).

pm_(C1-C2) -->
    (	[declare_mtd(_)]
    ;	[mtd_notfound(_)]
    ).

ground_or_nil(Term) :- ground(Term).
ground_or_nil('NA') :- true.

%?- phrase(pm_([A,B]-[C,D]), [1^0, 2^1, 2-0, mtd_notfound(2)]).
%@ A = D, D = 0,
%@ B = 1 ;
%@ false.

%?- path_matrix_([1^0, 2^1, 2-0, mtd_notfound(2)], M).
%@ M = [0, 1|_3816]-[_3826, 0|_3834].

%?- phrase(pm_(C1-C2), [1^0, 2^1, 2-0, mtd_notfound(2)]).
%@ C1 = [0, 1|_916],
%@ C2 = [_926, 0|_934] ;
%@ false.

pathmatrix(First-Second) -->
    row(First),
    row(Second).

row(Ls) -->
    { length(Ls,D),
      columns_format(D,Format) },
    format_(Format, Ls).

columns_format(1, '~w~n').
columns_format(N, F) :-
    N #> 1,
    N_1 #= N - 1,
    columns_format(N_1, F_1),
    atom_concat('~w\t', F_1, F).

%?- columns_format(6, Format).
%@ Format = '~w\t~w\t~w\t~w\t~w\t~w' ;
%@ false.

% Approximate format_//2 as provided by Scryer's library(format):
format_(Format, Ls) --> [ FLs ], { format(atom(FLs), Format, Ls) }.

format_matrix(Matrix) :- format('~s~n~s~n', Matrix).

%?- phrase(pathmatrix([1,'NA']-['NA','NA']), PM), format_matrix(PM).
%@ 1  NA 
%@ NA NA 
%@ PM = ['1  NA ', 'NA NA '].

%?- phrase(esc(0, 0..2), Tr), path_matrix(Tr, 2, M).
%@ Tr = [1^0, 2^0, 2*0, mtd_notfound(2)],
%@ M = [0, 0]-['NA', 0] ;
%@ Tr = [1^0, 2^0, 2*1, mtd_notfound(2)],
%@ M = [0, 0]-['NA', 1] ;
%@ Tr = [1^0, 2^0, 2*2, 1:0, declare_mtd(1)],
%@ M = [0, 0]-[0, 2] ;
%@ Tr = [1^0, 2^0, 2*2, 1:1, declare_mtd(1)],
%@ M = [0, 0]-[1, 2] ....

%?- phrase(esc(0, 0..2), Tr), path_matrix(Tr, 2, M), phrase(pathmatrix(M), PM), format_matrix(PM).
%@ 0  0  
%@ NA 0  
%@ Tr = [1^0, 2^0, 2*0, mtd_notfound(2)],
%@ M = [0, 0]-['NA', 0],
%@ PM = ['0  0  ', 'NA 0  '] ;
%@ 0  0  
%@ NA 1  
%@ Tr = [1^0, 2^0, 2*1, mtd_notfound(2)],
%@ M = [0, 0]-['NA', 1],
%@ PM = ['0  0  ', 'NA 1  '] ;
%@ 0  0  
%@ 0  2  
%@ Tr = [1^0, 2^0, 2*2, 1:0, declare_mtd(1)],
%@ M = [0, 0]-[0, 2],
%@ PM = ['0  0  ', '0  2  '] ;
%@ 0  0  
%@ 1  2  
%@ Tr = [1^0, 2^0, 2*2, 1:1, declare_mtd(1)],
%@ M = [0, 0]-[1, 2],
%@ PM = ['0  0  ', '1  2  '] ....

rep([], 0, _).
rep([E|Es], N, E) :-
    N #> 0,
    N_1 #= N - 1,
    rep(Es, N_1, E).

%?- rep([a,a,a,a,a], N, X).
%@ N = 5,
%@ X = a.

?- rep(Es, N, E).
%@ Es = [],
%@ N = 0 ;
%@ Es = [E],
%@ N = 1 ;
%@ Es = [E, E],
%@ N = 2 ;
%@ Es = [E, E, E],
%@ N = 3 ;
%@ Es = [E, E, E, E],
%@ N = 4 ...
% Dang!

%?- rep(a, 5, As).
%@ As = [a, a, a, a, a].

%?- repe(a, 5, As).
%@ As = [a, a, a, a, a] ;
%@ false.

%% Write out tab-delimited input files T<D>.tab
write_T(D) :-
    format(atom(Filename), 'T~d.tab', D),
    open(Filename, write, OS),
    findall(P, phrase(esc(0, 0..D), P), Paths),
    length(Paths, Len),
    format('~d ~4t ~d~n', [D, Len]), % feedback to console
    rep(Ds, Len, D),
    maplist(path_matrix, Paths, Ds, Ms),
    with_output_to(OS, see_esc_mtx(D, Ms)) -> close(OS).

see_esc_mtx(D, Ms) :-
    columns_format(D, Format),
    phrase(see_esc_mtx_(Format), Ms).

see_esc_mtx_(Format) -->
    [C1-C2],
    { format(Format, C1),
      format(Format, C2) },
    see_esc_mtx_(Format).
see_esc_mtx_(_) --> [].

%?- maplist(write_T, [2,3,4,5,6,7,8]).
%@ 2  46
%@ 3  154
%@ 4  442
%@ 5  1162
%@ 6  2890
%@ 7  6922
%@ 8  16138
%@ true.
%@ 2  46
%@ 3  154
%@ 4  442
%@ 5  1162
%@ 6  2890
%@ 7  6922
%@ 8  16138
%@ true.
%@ 2  46
%@ 3  154
%@ 4  442
%@ 5  1162
%@ 6  2890
%@ 7  6922
%@ 8  16138
%@ true.
