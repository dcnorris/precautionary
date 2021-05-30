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
	T1 + max(0, N2 - N1) #=< T2 #<==> B,
	zo_t(B, Truth)
    ).

qcompare(>=, T1/N1, T2/N2, Truth) :-
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
    (	Q &>= C -> Truth = true
    ;	hit_ceiling_t(Q, Cs, Truth)
    ).

hit_floor_t(_, [], false).
hit_floor_t(Q, [F|Fs], Truth) :-
    (	Q &=< F -> Truth = true
    ;	hit_floor_t(Q, Fs, Truth)
    ).


/*

Our 'spatial' understanding of floors/ceilings defined with respect to
(&=) and (&>=) shows that their minimal representation is obtained by
listing only their (outer) vertices.

*/

ceiling_vertex_t(Qs, T/N, Truth) :-
    memberd_t(T/N, Qs, true),
    N1 #= N + 1,
    if_(hit_ceiling_t(T/N1, Qs)
	, Truth = false
	, ( Tminus1 #= T - 1,
	    Nminus1 #= N - 1,
	    if_(hit_ceiling_t(Tminus1/Nminus1, Qs)
		, Truth = false
		, Truth = true
	       )
	  )
       ).

%?- time(ceiling_vertex_t([3/3,3/4,3/5,4/6,4/7,4/8,5/9,5/10,6/11,6/12], V, true)).
%@    % CPU time: 0.017 seconds
%@    V = 3/5
%@ ;  % CPU time: 0.061 seconds
%@    V = 4/8
%@ ;  % CPU time: 0.105 seconds
%@    V = 5/10
%@ ;  % CPU time: 0.154 seconds
%@    V = 6/12
%@ ;  % CPU time: 0.174 seconds
%@    false. %    ^^^^^ using memberd_t(T/N, Qs, true)
%@    % CPU time: 0.011 seconds
%@    V = 3/5
%@ ;  % CPU time: 0.031 seconds
%@    V = 4/8
%@ ;  % CPU time: 0.052 seconds
%@    V = 5/10
%@ ;  % CPU time: 0.074 seconds
%@    V = 6/12
%@ ;  % CPU time: 0.086 seconds
%@    false. %    ^^^^^ using member(T/N, Qs)

%?- ceiling_vertex_t([3/3,3/4,3/5,4/6,5/7,6/8,6/11,6/12], V, true).
%@    V = 3/5
%@ ;  V = 6/12
%@ ;  false.

floor_vertex_t(Qs, T/N, Truth) :-
    memberd_t(T/N, Qs, true),
    Nminus1 #= N - 1,
    if_(hit_floor_t(T/Nminus1, Qs)
	, Truth = false
	, ( T1 #= T + 1,
	    N1 #= N + 1,
	    if_(hit_floor_t(T1/N1, Qs)
		, Truth = false
		, Truth = true
	       )
	  )
       ).

%?- floor_vertex_t([0/3,1/4,2/5,4/8,5/12], V, true).
%@    V = 2/5
%@ ;  V = 4/8
%@ ;  V = 5/12
%@ ;  false.

%?- floor_vertex_t([0/3,1/5,4/8,5/12], V, true).
%@    V = 0/3
%@ ;  V = 4/8
%@ ;  V = 5/12
%@ ;  false.

%?- floor_vertex_t([0/3,1/5,4/8,5/12], V, false).
%@    V = 1/5
%@ ;  false.

%?- tfilter(floor_vertex_t([0/3,1/5,4/8,5/12]), [0/3,1/5,4/8,5/12], Vs).
%@    Vs = [0/3,4/8,5/12]
%@ ;  false.

%% It will help to have unique, minimal ('canonical') representations
%% for floor- and ceiling-type boundaries.
%% TODO: Are these concepts better represented through SETS than lists?
%%       Shouldn't I require that each tally in a canonical floor/ceiling
%%       contributes something? Does this requirement ensure uniqueness?
ceiling_canonical(Qs, Ks) :-
    tfilter(ceiling_vertex_t(Qs), Qs, Ks_),
    sort(Ks_, Ks).

%?- ceiling_canonical([3/3,3/4,3/5,4/6,4/7,4/8,5/9,5/10,6/11,6/12], K).
%@    K = [3/5,4/8,5/10,6/12]
%@ ;  false.

%?- ceiling_canonical([6/12,4/8,3/3,3/4,3/5,4/6,4/7,5/9,5/10,6/11], K).
%@    K = [3/5,4/8,5/10,6/12]
%@ ;  false.

floor_canonical(Qs, Ks) :-
    tfilter(floor_vertex_t(Qs), Qs, Ks_),
    sort(Ks_, Ks).

%?- floor_canonical([0/1, 0/2, 0/3, 0/4, 0/5, 0/6, 0/7, 1/8, 1/9, 1/10, 1/11, 1/12], K).
%@    K = [0/1,1/8]
%@ ;  false.

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
    (	hit_ceiling_t(Q, RemovalBdy, true) -> Decision = remove
    ;	hit_ceiling_t(Q, DeescBdy, true) -> Decision = deescalate
    ;	hit_floor_t(Q, EscBdy, true) -> Decision = escalate
    ;	Decision = stay
    ).

%% For testing purposes, we hard-code the CCD to obtain tally_decision/2
tally_decision(Q, Decision) :-
    tally_decision_ccd(Q, Decision, ccd([3/5, 4/8, 5/10, 6/12],
					[1/1, 2/4, 3/6, 4/8, 5/11],
					[0/1, 1/8],
					12)).

%?- tally_decision(T/5, Decision).
%@    Decision = remove, clpz:(T in 3..5).
%@    Decision = stay, clpz:(T in 1..2)
%@ ;  Decision = escalate, T = 0
%@ ;  Decision = remove, clpz:(T in 3..5).

%?- tally_decision(Q, Decision).
%@    Q = 0/0, Decision = stay
%@ ;  Q = 1/1, Decision = deescalate
%@ ;  Q = 2/2, Decision = deescalate
%@ ;  Q = 3/3, Decision = remove
%@ ;  Q = _A/4, Decision = remove, clpz:(_A in 3..4)
%@ ;  Q = _A/5, Decision = remove, clpz:(_A in 3..5)
%@ ;  Q = _A/6, Decision = remove, clpz:(_A in 4..6)
%@ ;  Q = _A/7, Decision = remove, clpz:(_A in 5..7)
%@ ;  Q = _A/8, Decision = remove, clpz:(_A in 6..8)
%@ ;  Q = _A/9, Decision = remove, clpz:(_A in 7..9)
%@ ;  Q = _A/10, Decision = remove, clpz:(_A in 8..10)
%@ ;  Q = _A/11, Decision = remove, clpz:(_A in 9..11)
%@ ;  Q = _A/12, Decision = remove, clpz:(_A in 10..12).
%@    Q = 0/0, Decision = stay
%@ ;  Q = 0/1, Decision = escalate
%@ ;  Q = 1/1, Decision = deescalate
%@ ;  Q = 1/2, Decision = stay
%@ ;  Q = 0/2, Decision = escalate
%@ ;  Q = 2/2, Decision = deescalate
%@ ;  Q = 1/3, Decision = stay
%@ ;  Q = 0/3, Decision = escalate
%@ ;  Q = 2/3, Decision = deescalate
%@ ;  Q = 3/3, Decision = remove
%@ ;  Q = 1/4, Decision = stay
%@ ;  Q = 0/4, Decision = escalate
%@ ;  Q = 2/4, Decision = deescalate
%@ ;  Q = _A/4, Decision = remove, clpz:(_A in 3..4)
%@ ;  Q = _A/5, Decision = stay, clpz:(_A in 1..2)
%@ ;  Q = 0/5, Decision = escalate
%@ ;  ...


%% This enroll/3 goal, with its reified 'success' arg #3, creates fine
%% opportunities to bring various CCD-adapted STOPPING CRITERIA to bear.
%% Presently, we are simply checking whether cohort N0 is already 'full'.
enroll(T0/N0, T1/N1, Truth) :-
    (	N0 #>= 6 -> Truth = false
    ;	N1 #= N0 + 1,
	T in 0..1, % T is the pending tox assessment of newly-enrolled patient
	indomain(T), % TODO: How to 'parametrize' this? Use OPTIONS?
	T1 #= T0 + T,
	Truth = true
    ).

length_plus_1(Ls, MTD) :-
    length(Ls, MTDminus1),
    MTD #= MTDminus1 + 1.

stay(Ls ^ [R | Rs] ^ Es, State) :-
    enroll(R, R1, Truth),
    (	Truth == true -> State = Ls ^ [R1 | Rs] ^ Es
    ;	length_plus_1(Ls, MTD),
	State = declare_mtd(MTD)
    ).

escalate(Ls ^ [R] ^ Es, State) :- % NB: this is a 'clamped' situation
    stay(Ls ^ [R] ^ Es, State).

escalate(Ls ^ [Q, R | Rs] ^ Es, State) :-
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

deescalate([] ^ _ ^ _, declare_mtd(0)). % deescalate from already-lowest dose

deescalate([L | Ls] ^ Rs ^ Es, State) :-
    enroll(L, L1, Truth),
    (	Truth == true -> State = Ls ^ [L1 | Rs] ^ Es
    ;	( length_plus_1(Ls, MTD),
	  State = declare_mtd(MTD)
	)
    ).

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

%% Examine the smallest possible trial -- a trial with just 1 dose!
%?- Trial+\(default_ccd(CCD), phrase(ccd_actions(CCD, []^[0/0]^[]), Trial)).
%@    Trial = [([]^[0/0]^[]->stay->[]^[0/1]^[]),([]^[0/1]^[]->escalate->[]^[0/2]^[]),([]^[0/2]^[]->escalate->[]^[0/3]^[]),([]^[0/3]^[]->escalate->[]^[0/4]^[]),([]^[0/4]^[]->escalate->[]^[0/5]^[]),([]^[0/5]^[]->escalate->[]^[0/6]^[]),([]^[0/6]^[]->escalate->declare_mtd(1))]
%@ ;  Trial = [([]^[0/0]^[]->stay->[]^[0/1]^[]),([]^[0/1]^[]->escalate->[]^[0/2]^[]),([]^[0/2]^[]->escalate->[]^[0/3]^[]),([]^[0/3]^[]->escalate->[]^[0/4]^[]),([]^[0/4]^[]->escalate->[]^[0/5]^[]),([]^[0/5]^[]->escalate->[]^[1/6]^[]),([]^[1/6]^[]->stay->declare_mtd(1))]
%@ ;  Trial = [([]^[0/0]^[]->stay->[]^[0/1]^[]),([]^[0/1]^[]->escalate->[]^[0/2]^[]),([]^[0/2]^[]->escalate->[]^[0/3]^[]),([]^[0/3]^[]->escalate->[]^[0/4]^[]),([]^[0/4]^[]->escalate->[]^[1/5]^[]),([]^[1/5]^[]->stay->[]^[1/6]^[]),([]^[1/6]^[]->stay->declare_mtd(1))]
%@ ;  ...

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

%?- Matrix+\(default_ccd(CCD), ccd_d_matrix(CCD, 2, Matrix)).
%@    Matrix = [[0/1]^[0/6]^[]~>2]
%@ ;  Matrix = [[0/1]^[1/6]^[]~>2]
%@ ;  Matrix = [[0/1]^[1/6]^[]~>2]
%@ ;  ...

%?- J+\(default_ccd(CCD), D=1, time(findall(M, ccd_d_matrix(CCD, D, M), Ms)), length(Ms, J)).
%@    % CPU time: 0.181 seconds
%@    % CPU time: 0.185 seconds
%@    J = 20.

%?- J+\(default_ccd(CCD), D=2, time(findall(M, ccd_d_matrix(CCD, D, M), Ms)), length(Ms, J)).
%@    % CPU time: 1.982 seconds
%@    % CPU time: 1.986 seconds
%@    J = 212.

%?- J+\(default_ccd(CCD), D=3, time(findall(M, ccd_d_matrix(CCD, D, M), Ms)), length(Ms, J)).
%@    % CPU time: 10.613 seconds
%@    % CPU time: 10.617 seconds
%@    J = 1151.

%?- J+\(default_ccd(CCD), D=4, time(findall(M, ccd_d_matrix(CCD, D, M), Ms)), length(Ms, J)).
%@    % CPU time: 62.156 seconds
%@    % CPU time: 62.161 seconds
%@    J = 6718.

%?- J+\(default_ccd(CCD), D=5, time(findall(M, ccd_d_matrix(CCD, D, M), Ms)), length(Ms, J)).
%@    % CPU time: 370.651 seconds
%@    % CPU time: 370.655 seconds
%@    J = 39289.
