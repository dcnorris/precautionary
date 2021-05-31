%% Attempt to get a small repro of the if_ vs if-then performance issue

:- use_module(library(clpz)).
:- use_module(library(lists)).
:- use_module(library(dcgs)).
:- use_module(library(time)).
:- use_module(library(debug)).
:- use_module(library(dif)).
:- use_module(library(reif)).
:- use_module(library(format)).
:- use_module(library(lambda)).


qcompare(=<, T1/N1, T2/N2) :-
    T1 + max(0, N2 - N1) #=< T2.
	
qcompare(>=, T1/N1, T2/N2) :-
    T1 #>= T2 + max(0, N1 - N2).

%% Reified versions of the above, as done at bottom of clpz.pl
qcompare(=<, T1/N1, T2/N2, Truth) :-
    /*
    %% Let's try using fast arithmetic, when possible
    (	ground(T1/N1 - T2/N2) ->
	(   DN is N2 - N1,
	    (	DN >= 0 ->
		(   T1plusDN is T1 + DN,
		    T1plusDN =< T2 -> Truth = true
		;   Truth = false
		)
	    ;	% DN < 0, so a simpler condition applies
		(   T1 =< T2 -> Truth = true
		;   Truth = false
		)
	    )
	)
    ;	% general case (non-ground args 2 or 3) is handled by CLP(Z) ...
    */(	
	T1 + max(0, N2 - N1) #=< T2 #<==> B,
	zo_t(B, Truth)
    ).

qcompare(>=, T1/N1, T2/N2, Truth) :-
    /*
    %% Let's try using fast arithmetic, when possible
    (	ground(T1/N1 - T2/N2) ->
	(   DN is N1 - N2,
	    (	DN >= 0 ->
		(   T2plusDN is T2 + DN,
		    T2plusDN =< T1 -> Truth = true
		;   Truth = false
		)
	    ;	% DN < 0, so a simpler condition applies
		(   T1 >= T2 -> Truth = true
		;   Truth = false
		)
	    )
	)
    ;	% general case (non-ground args 2 or 3) is handled by CLP(Z) ...
    */(	
	T1 #>= T2 + max(0, N1 - N2) #<==> B,
	zo_t(B, Truth)
    ).

zo_t(0, false).
zo_t(1, true).


:- op(900, xfx, &=<).
&=<(Q1, Q2) :-
    qcompare(=<, Q1, Q2).

&=<(Q1, Q2, Truth) :- % reified
    qcompare(=<, Q1, Q2, Truth).

:- op(900, xfx, &>=).
&>=(Q1, Q2) :-
    qcompare(>=, Q1, Q2).

&>=(Q1, Q2, Truth) :- % reified
    qcompare(>=, Q1, Q2, Truth).


hit_ceiling_t(_, [], false).
hit_ceiling_t(Q, [C|Cs], Truth) :-
    if_(Q &>= C
	, Truth = true
	, hit_ceiling_t(Q, Cs, Truth)
       ).
/*
    (	Q &>= C -> Truth = true
    ;	hit_ceiling_t(Q, Cs, Truth)
    ).
*/

hit_floor_t(_, [], false).
hit_floor_t(Q, [F|Fs], Truth) :-
    if_(Q &=< F
	, Truth = true
	, hit_floor_t(Q, Fs, Truth)
       ).
/*
    (	Q &=< F -> Truth = true
    ;	hit_floor_t(Q, Fs, Truth)
    ).
*/

%% tally_decision_ccd(?Q, ?Decision, +CCD) relates tallies Q to Decisions,
%% for a GIVEN ground cumulative-cohort design (CCD) which takes the form
%% of a triplet of boundaries, followed by max enrollment per cohort.
%% TODO: I believe the if-then cascade, with its default final 'escape clause',
%%       together with the determinism of hit_ceiling_t/3 and hit_floor_t/3,
%%       proves the determinism of tally_decision/2 so long as the defining
%%       boundaries are lists from ℚ.
tally_decision_ccd(Q, Decision, ccd(RemovalBdy, DeescBdy, EscBdy, FullCoh)) :-
    Q = T/N,
    N in 0..FullCoh, indomain(N),
    T in 0..N,
    if_(hit_ceiling_t(Q, RemovalBdy)
	, Decision = remove
	, if_(hit_ceiling_t(Q, DeescBdy)
	      , Decision = deescalate
	      , if_(hit_floor_t(Q, EscBdy)
		    , Decision = escalate
		    , Decision = stay
		   )
	     )
       ).
/*
    (	hit_ceiling_t(Q, RemovalBdy, true) -> Decision = remove
    ;	hit_ceiling_t(Q, DeescBdy, true) -> Decision = deescalate
    ;	hit_floor_t(Q, EscBdy, true) -> Decision = escalate
    ;	Decision = stay
    ).
*/

%% For testing purposes, we hard-code the CCD to obtain tally_decision/2
tally_decision(Q, Decision) :-
    tally_decision_ccd(Q, Decision, ccd([3/5, 4/8, 5/10, 6/12],
					[1/1, 2/4, 3/6, 4/8, 5/11],
					[0/1, 1/8],
					12)).


%% This enroll/3 goal, with its reified 'success' arg #3, creates fine
%% opportunities to bring various CCD-adapted STOPPING CRITERIA to bear.
%% Presently, we are simply checking whether cohort N0 is already 'full'.
enroll(T0/N0, T1/N1, Truth) :-
    N0 #>= 6 #<==> CohortFull, % suggestive of a line from a trial 'config file'?
    if_(CohortFull #= 1
       , Truth = false
       , ( N1 #= N0 + 1,
           T in 0..1, % T is the pending tox assessment of newly-enrolled patient
           indomain(T), % TODO: How to 'parametrize' this? Use OPTIONS?
           T1 #= T0 + T,
           Truth = true
         )
       ).
/*
    (	N0 #>= 6 -> Truth = false
    ;	N1 #= N0 + 1,
	T in 0..1, % T is the pending tox assessment of newly-enrolled patient
	indomain(T), % TODO: How to 'parametrize' this? Use OPTIONS?
	T1 #= T0 + T,
	Truth = true
    ).
*/

length_plus_1(Ls, MTD) :-
    length(Ls, MTDminus1),
    MTD #= MTDminus1 + 1.

stay(Ls ^ [R | Rs] ^ Es, State) :-
    if_(enroll(R, R1)
	, State = Ls ^ [R1 | Rs] ^ Es
	, ( length_plus_1(Ls, MTD),
	    State = declare_mtd(MTD)
	  )
       ).
/*
    enroll(R, R1, Truth),
    (	Truth == true -> State = Ls ^ [R1 | Rs] ^ Es
    ;	length_plus_1(Ls, MTD),
	State = declare_mtd(MTD)
    ).
*/

escalate(Ls ^ [R] ^ Es, State) :- % NB: this is a 'clamped' situation
    stay(Ls ^ [R] ^ Es, State).

escalate(Ls ^ [Q, R | Rs] ^ Es, State) :-
    if_(enroll(R, R1)
	, State = [Q | Ls] ^ [R1 | Rs] ^ Es
	, ( length_plus_1(Ls, MTD),
	    State = declare_mtd(MTD)
	  )
       ).
/*
    enroll(R, R1, Truth),
    (	Truth == true -> State = [Q | Ls] ^ [R1 | Rs] ^ Es
    ;	%% If the next dose up (R) cannot be enrolled, that's because
	%% it's already full. What's more, it must have recommended
	%% de-escalation---which is how we got to the current dose Q!
	%% Accordingly, we declare the current dose to be MTD:
	( length_plus_1(Ls, MTD),
	  State = declare_mtd(MTD)
	)
    ).
*/

deescalate([] ^ _ ^ _, declare_mtd(0)). % deescalate from already-lowest dose

deescalate([L | Ls] ^ Rs ^ Es, State) :-
    if_(enroll(L, L1)
	, State = Ls ^ [L1 | Rs] ^ Es
	, ( length_plus_1(Ls, MTD),
	    State = declare_mtd(MTD)
	  )
       ).
/*
    enroll(L, L1, Truth),
    (	Truth == true -> State = Ls ^ [L1 | Rs] ^ Es
    ;	( length_plus_1(Ls, MTD),
	  State = declare_mtd(MTD)
	)
    ).
*/

remove(Ls ^ Rs ^ Es, State) :-
    append(Rs, Es, RsEs),
    deescalate(Ls ^ [] ^ RsEs, State).

ccd_state0_action_state(CCD, Ls ^ [R | Rs] ^ Es, Action, State) :-
    %% TODO: Check whether total enrollment has been reached, and stop?
    tally_decision_ccd(R, Action, CCD),
    call(Action, Ls ^ [R | Rs] ^ Es, State).

%% Note how this state-machine naturally starts up from a blank slate:
%?- (Action-State)+\(default_ccd(CCD), ccd_state0_action_state(CCD, [] ^ [0/0, 0/0, 0/0] ^ [], Action, State)).
%@    Action = stay, State = []^[0/1,0/0,0/0]^[]
%@ ;  Action = stay, State = []^[1/1,0/0,0/0]^[].

ccd_actions(_, declare_mtd(_)) --> [].
ccd_actions(CCD, S0) --> [(S0->A->S)],
			 { ccd_state0_action_state(CCD, S0, A, S) },
			 ccd_actions(CCD, S).

:- op(900, xfx, ~>).

%% First, extract the path matrix:
path_matrix, [S ~> MTD] --> [ (S->_->declare_mtd(MTD)) ].
path_matrix, [] --> [(_->_->_^_^_)], path_matrix. % 'skip to end'

default_ccd(ccd([3/5, 4/8, 5/10, 6/12],		
		[1/1, 2/4, 3/6, 4/8, 5/11],
		[0/1, 1/8],
		12)).

%% I've implemented this predicate to more clearly exhibit the
%% sharing of arguments with ccd_state0_action_state/4.
ccd_state0_matrix(CCD, State0, Matrix) :-
    phrase(ccd_actions(CCD, State0), Path),
    phrase(path_matrix, Path, Matrix).

%% This predicate implements the common special case
%% where a CCD trial starts from the lowest dose.
ccd_d_matrix(CCD, D, Matrix) :-
    length(Tallies, D), maplist(=(0/0), Tallies),
    ccd_state0_matrix(CCD, []^Tallies^[], Matrix).

ccd_d_path(CCD, D, Path) :-
    length(Tallies, D), maplist(=(0/0), Tallies),
    phrase(ccd_actions(CCD, []^Tallies^[]), Path).

%?- Matrix+\(default_ccd(CCD), ccd_d_matrix(CCD, 2, Matrix)).
%@    Matrix = [[0/1]^[0/6]^[]~>2]
%@ ;  Matrix = [[0/1]^[1/6]^[]~>2]
%@ ;  Matrix = [[0/1]^[1/6]^[]~>2]
%@ ;  Matrix = [[0/1]^[2/6]^[]~>2]
%@ ;  ...

%?- J+\(default_ccd(CCD), D=1, time(findall(M, ccd_d_path(CCD, D, P), Ps)), length(Ps, J)).
%@    % CPU time: 10.140 seconds
%@    % CPU time: 10.144 seconds
%@    J = 20. % ^ removing fast arithmetic branches from qcompare/4
%@    % CPU time: 1.216 seconds
%@    % CPU time: 1.220 seconds
%@    J = 20. % ^ ... now in enroll/3
%@    % CPU time: 0.405 seconds
%@    % CPU time: 0.410 seconds
%@    J = 20. % ^ now with if_ in escalate/2 too
%@    % CPU time: 0.414 seconds
%@    % CPU time: 0.418 seconds
%@    J = 20. % ^ now with if_ in deescalate/2 as well
%@    % CPU time: 0.396 seconds
%@    % CPU time: 0.401 seconds
%@    J = 20. % ^ swapping if_ into stay/2
%@    % CPU time: 0.388 seconds
%@    % CPU time: 0.392 seconds
%@    J = 20. % ^ swapping if_ into hit_ceilingfloor_t/3 predicates
%@    % CPU time: 0.273 seconds
%@    % CPU time: 0.277 seconds
%@    J = 20. % ^ swapping if_ into tally_decision_ccd/3
%@    % CPU time: 0.185 seconds
%@    % CPU time: 0.189 seconds
%@    J = 20.

%?- J+\(default_ccd(CCD), D=2, time(findall(M, ccd_d_path(CCD, D, P), Ps)), length(Ps, J)).
%@    % CPU time: 108.094 seconds
%@    % CPU time: 108.098 seconds
%@    J = 212. % ^ removing fast arithmetic branches from qcompare/4
%@    % CPU time: 13.608 seconds
%@    % CPU time: 13.613 seconds
%@    J = 212. % ^ ... now in enroll/3
%@    % CPU time: 4.272 seconds
%@    % CPU time: 4.277 seconds
%@    J = 212. % ^ now with if_ in escalate/2 too
%@    % CPU time: 4.269 seconds
%@    % CPU time: 4.274 seconds
%@    J = 212. % ^ now with if_ in deescalate/2 as well
%@    % CPU time: 4.240 seconds
%@    % CPU time: 4.244 seconds
%@    J = 212. % ^ swapping if_ into stay/2
%@    % CPU time: 4.049 seconds
%@    % CPU time: 4.053 seconds
%@    J = 212. % ^ swapping if_ into hit_ceilingfloor_t/3 predicates
%@    % CPU time: 2.727 seconds
%@    % CPU time: 2.731 seconds
%@    J = 212. % ^ swapping if_ into tally_decision_ccd/3
%@    % CPU time: 1.887 seconds
%@    % CPU time: 1.891 seconds
%@    J = 212.

%?- J+\(default_ccd(CCD), D=3, time(findall(M, ccd_d_path(CCD, D, P), Ps)), length(Ps, J)).
%@    % CPU time: 606.836 seconds
%@    % CPU time: 606.840 seconds
%@    J = 1151. % ^ removing fast arithmetic branches from qcompare/4
%@    % CPU time: 73.445 seconds
%@    % CPU time: 73.449 seconds
%@    J = 1151. % ^ ... now in enroll/3
%@    % CPU time: 24.449 seconds
%@    % CPU time: 24.453 seconds
%@    J = 1151. % ^ now with if_ in escalate/2 too
%@    % CPU time: 23.212 seconds
%@    % CPU time: 23.217 seconds
%@    J = 1151. % ^ now with if_ in deescalate/2 as well
%@    % CPU time: 23.018 seconds
%@    % CPU time: 23.022 seconds
%@    J = 1151. % ^ swapping if_ into stay/2
%@    % CPU time: 22.220 seconds
%@    % CPU time: 22.224 seconds
%@    J = 1151. % ^ swapping if_ into hit_ceilingfloor_t/3 predicates
%@    % CPU time: 14.681 seconds
%@    % CPU time: 14.685 seconds
%@    J = 1151. % ^ swapping if_ into tally_decision_ccd/3
%@    % CPU time: 10.430 seconds
%@    % CPU time: 10.434 seconds
%@    J = 1151.

%?- J+\(default_ccd(CCD), D=4, time(findall(M, ccd_d_path(CCD, D, P), Ps)), length(Ps, J)).
%@    % CPU time: 132.071 seconds
%@    % CPU time: 132.076 seconds
%@    J = 6718. % ^ swapping if_ into hit_ceilingfloor_t/3 predicates
%@    % CPU time: 85.553 seconds
%@    % CPU time: 85.558 seconds
%@    J = 6718. % ^ swapping if_ into tally_decision_ccd/3
%@    % CPU time: 61.984 seconds
%@    % CPU time: 61.989 seconds
%@    J = 6718.

%?- X is 3.14 * 2.
%@    X = 6.28.
