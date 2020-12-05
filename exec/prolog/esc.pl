% Enumerate 3+3 trial outcomes for exact matrix calculations
%%:- use_module(library(clpfd)). % SWI
/* Scryer */
:- use_module(library(format)).
:- use_module(library(lists)).
:- use_module(library(dcgs)).
:- use_module(library(clpz)).
:- use_module(library(time)).
/**/
:- set_prolog_flag(double_quotes, chars).

% Prefix op * for 'generalizing away' goals (https://www.metalevel.at/prolog/debugging)
:- op(920, fy, *). *_.  % *Goal always succeeds

% DCG esc(D, Lo..Hi) describes a list of 3+3 cohorts FOLLOWING a dose D in Lo..Hi.
% One can read esc(D, Lo..Hi) as the DECISION to escalate *from* D to min(D+1,Hi).
% For example, esc(0, 0..5) *initiates* a trial that has 5 prespecified doses,
% and enrolls the first cohort at D=1.

% - - - - - begin inset for paper - - - - -

tox(T) :- T in 0..3,
	  indomain(T).

% Mnemonic: * is ^ that 'splatted' on dose ceiling.
esc(Hi,Lo..Hi) --> [Hi * T], { tox(T) },
		   (  {T #=< 1}, [mtd_notfound(Hi)]
		   ;  {T #>= 2}, des(Hi, Lo)
		   ).
esc(D, Lo..Hi) --> { D #< Hi, D1 #= D + 1 },
		   [D1 ^ T], { tox(T) },
		   (  {T #= 0}, esc(D1, Lo..Hi)
		   ;  {T #= 1}, sta(D1, Lo..Hi)
		   ;  {T #> 1}, des(D1, Lo)
		   ).

sta(D,  _..D ) --> [D - 0], [mtd_notfound(D)].
sta(D, Lo..Hi) --> { D #< Hi, D in Lo..Hi },
		   [D - 0],
		   esc(D, D..Hi).
sta(D, Lo.._ ) --> [D - T], { tox(T), T #> 0 },
		   des(D, Lo).

% - - - - - end inset for paper - - - - -

% As a mirror image of esc//2, des(D, Lo) moves
% downward FROM D, to max(D-1,Lo).
% NB: De-escalation to D-1 clamps Hi #= D - 1.
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

%% SWI vs Scryer timings ...

%?- time(n_trials(1..10, D_N)).
%@ % 1,247 inferences, 0.000 CPU in 0.000 seconds (93% CPU, 5517699 Lips)
%@ D_N =  (1, 10) ;
%@ % 5,453 inferences, 0.001 CPU in 0.001 seconds (98% CPU, 6012128 Lips)
%@ D_N =  (2, 46) ;
%@ % 18,291 inferences, 0.002 CPU in 0.003 seconds (99% CPU, 7334002 Lips)
%@ D_N =  (3, 154) ;
%@ % 52,471 inferences, 0.008 CPU in 0.008 seconds (96% CPU, 6736552 Lips)
%@ D_N =  (4, 442) ;
%@ % 137,823 inferences, 0.017 CPU in 0.018 seconds (96% CPU, 7876500 Lips)
%@ D_N =  (5, 1162) ;
%@ % 342,511 inferences, 0.043 CPU in 0.044 seconds (98% CPU, 7929413 Lips)
%@ D_N =  (6, 2890) ;
%@ % 819,855 inferences, 0.102 CPU in 0.105 seconds (97% CPU, 8012421 Lips)
%@ D_N =  (7, 6922) ;
%@ % 1,910,479 inferences, 0.222 CPU in 0.230 seconds (97% CPU, 8601267 Lips)
%@ D_N =  (8, 16138) ;
%@ % 4,363,599 inferences, 0.496 CPU in 0.519 seconds (95% CPU, 8797313 Lips)
%@ D_N =  (9, 36874) ;
%@ % 9,813,571 inferences, 1.130 CPU in 1.162 seconds (97% CPU, 8686744 Lips)
%@ D_N =  (10, 82954).
%@    % CPU time: 0.025 seconds
%@    D_N = (1,10)
%@ ;  % CPU time: 0.128 seconds
%@    D_N = (2,46)
%@ ;  % CPU time: 0.454 seconds
%@    D_N = (3,154)
%@ ;  % CPU time: 1.358 seconds
%@    D_N = (4,442)
%@ ;  % CPU time: 3.697 seconds
%@    D_N = (5,1162)
%@ ;  % CPU time: 9.491 seconds
%@    D_N = (6,2890)
%@ ;  % CPU time: 23.205 seconds
%@    D_N = (7,6922)
%@ ;  % CPU time: 54.931 seconds
%@    D_N = (8,16138)
%@ ;  % CPU time: 127.262 seconds
%@    D_N = (9,36874)
%@ ;  % CPU time: 300.266 seconds
%@    D_N = (10,82954)
%@ ;  % CPU time: 300.346 seconds
%@    false.

%?- phrase(esc(0, 0..2), Tr).
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

%% Transform dose-escalation path lists to the arrays T(c,d,j).

path_matrix(P, D, M) :-
    phrase(pm_(M), P),
    M = C1-C2,
    length(C1,D),
    length(C2,D),
    maplist(ground_or_nil, C1),
    maplist(ground_or_nil, C2).

pm_(C1-C2) -->
    [D^T], % ascribe to D's 1st cohort
    { D0 #= D-1, nth0(D0, C1, T) },
    pm_(C1-C2).

pm_(C1-C2) -->
    (	[D-T]
    ;	[D*T]
    ;	[D:T]
    ), % ascribe to D's 2nd cohort
    { D0 #= D-1, nth0(D0, C2, T) },
    pm_(C1-C2).

pm_(_-_) -->
    (	[declare_mtd(_)]
    ;	[mtd_notfound(_)]
    ).

ground_or_nil(Term) :- ground(Term).
ground_or_nil('NA') :- true.

%?- phrase(pm_([A,B]-[C,D]), [1^0, 2^1, 2-0, mtd_notfound(2)]).
%@ caught: error(existence_error(procedure,phrase/2),phrase/2)
%@ A = D, D = 0,
%@ B = 1 ;
%@ false.

%?- path_matrix_([1^0, 2^1, 2-0, mtd_notfound(2)], M).
%@ M = [0, 1|_3816]-[_3826, 0|_3834].

%?- phrase(pm_(C1-C2), [1^0, 2^1, 2-0, mtd_notfound(2)]).
%@ C1 = [0, 1|_916],
%@ C2 = [_926, 0|_934] ;
%@ false.

columns_format(1, "~w~n").
columns_format(N, F) :-
    N #> 1,
    N_1 #= N - 1,
    columns_format(N_1, F_1),
    append("~w\t", F_1, F).

%% rep(?List, ?N, ?Elt) true if List is N repetitions of Elt.
rep([], 0, _).
rep([E|Es], N, E) :-
    N #> 0,
    N_1 #= N - 1,
    rep(Es, N_1, E).

%?- rep(Es, N, E).
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

%% Write out tab-delimited input files T<D>.tab
write_T(D) :-
    findall(P, phrase(esc(0, 0..D), P), Paths),
    length(Paths, Len),
    format("D = ~d ~t~8+J = ~d~n", [D, Len]), % feedback to console
    rep(Ds, Len, D),
    maplist(path_matrix, Paths, Ds, Ms) ->
	write_esc_array(D, Ms).

write_esc_array(D, Ms) :-
    %%format(atom(File), "T~d.tab", [D]), % SWI
    phrase(format_("T~d.tab", [D]), Filename), atom_chars(File, Filename), % Scryer
    open(File, write, OS),
    columns_format(D, Format),
    phrase(write_esc_array_(OS, Format), Ms) -> close(OS).

write_esc_array_(OS, Format) -->
    [C1-C2],
    { format(OS, Format, C1),
      format(OS, Format, C2) },
    write_esc_array_(OS, Format).
write_esc_array_(_, _) --> [].

%?- write_T(3).

%?- time(maplist(write_T, [2,3,4,5,6,7,8])).
%@ D = 2   J = 46
%@ D = 3   J = 154
%@ D = 4   J = 442
%@ D = 5   J = 1162
%@ D = 6   J = 2890
%@ D = 7   J = 6922
%@ D = 8   J = 16138
%@    % CPU time: 665.655 seconds
%@    true
%@ ;  ...
%@ % CPU time: 666.070 seconds
%@    
%@ D = 2   J = 46
%@ D = 3   J = 154
%@ D = 4   J = 442
%@ D = 5   J = 1162
%@ D = 6   J = 2890
%@ D = 7   J = 6922
%@ D = 8   J = 16138
%@ % 6,064,313 inferences, 1.132 CPU in 1.237 seconds (92% CPU, 5355090 Lips)
%@ true.

%%:- initialization maplist(write_T, [2,3,4,5,6,7,8]) -> halt.
