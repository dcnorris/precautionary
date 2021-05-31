%% Attempt to get a small repro of the if_ vs if-then performance issue

:- use_module(library(clpz)).
:- use_module(library(lists)).
:- use_module(library(dcgs)).
:- use_module(library(time)).
:- use_module(library(debug)).
:- use_module(library(dif)).
:- use_module(library(reif)).
:- use_module(library(lambda)).

qcompare(=<, T1/N1, T2/N2) :- T1 + max(0, N2 - N1) #=< T2.
qcompare(>=, T1/N1, T2/N2) :- T1 #>= T2 + max(0, N1 - N2).

%% Reified versions of the above, as done at bottom of clpz.pl
qcompare(=<, T1/N1, T2/N2, Truth) :-
    (	%% Use fast arithmetic, when possible
    	ground(T1/N1 - T2/N2) ->
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

qcompare(>=, Q1, Q2, Truth) :- qcompare(=<, Q2, Q1, Truth).

zo_t(0, false).
zo_t(1, true).


:- op(900, xfx, &=<).
&=<(Q1, Q2) :- qcompare(=<, Q1, Q2).
&=<(Q1, Q2, Truth) :- % reified
    qcompare(=<, Q1, Q2, Truth).

:- op(900, xfx, &>=).
&>=(Q1, Q2) :- qcompare(>=, Q1, Q2).
&>=(Q1, Q2, Truth) :- % reified
    qcompare(>=, Q1, Q2, Truth).


hit_ceiling_t(_, [], false).
hit_ceiling_t(Q, [C|Cs], Truth) :-
    if_(Q &>= C
	, Truth = true
	, hit_ceiling_t(Q, Cs, Truth)
       ).
/* Alternate formulation yielding 25% speedup
hit_ceiling_t(Q, [C|Cs], Truth) :-
    (	&>=(Q, C, true) -> Truth = true
    ;	hit_ceiling_t(Q, Cs, Truth)
    ).
*/

hit_floor_t(_, [], false).
hit_floor_t(Q, [F|Fs], Truth) :-
    if_(Q &=< F
	, Truth = true
	, hit_floor_t(Q, Fs, Truth)
       ).
/* Alternate formulation yielding 25% speedup
hit_floor_t(Q, [F|Fs], Truth) :-
    (	&=<(Q, F, true) -> Truth = true
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
/* I see about 10-15% savings from this substitution ... (but at what cost?)
    (	hit_ceiling_t(Q, RemovalBdy, true) -> Decision = remove
    ;	hit_ceiling_t(Q, DeescBdy, true) -> Decision = deescalate
    ;	hit_floor_t(Q, EscBdy, true) -> Decision = escalate
    ;	Decision = stay
    ).
*/

cohort_full(N, true) :- N >= 6.
cohort_full(N, false) :- N < 6.

%?- cohort_full(N, false).
%@ caught: error(instantiation_error,(is)/2)
%@ caught: error(instantiation_error,(is)/2)

%% This enroll/3 goal, with its reified 'success' arg #3, creates fine
%% opportunities to bring various CCD-adapted STOPPING CRITERIA to bear.
%% Presently, we are simply checking whether cohort N0 is already 'full'.
enroll(T0/N0, T1/N1, Truth) :-
    %%N0 #>= 6 #<==> CohortFull, % suggestive of a line from a trial 'config file'?
    %%if_(CohortFull #= 1
    if_(cohort_full(N0) % yields ~2.2x speedup vs reified CohortFull above
       , Truth = false
       , ( N1 #= N0 + 1,
           T in 0..1, % T is the pending tox assessment of newly-enrolled patient
           indomain(T), % TODO: How to 'parametrize' this? Use OPTIONS?
           T1 #= T0 + T,
           Truth = true
         )
       ).

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

escalate(Ls ^ [R] ^ Es, State) :- % NB: this is a 'clamped' situation
    stay(Ls ^ [R] ^ Es, State).

escalate(Ls ^ [Q, R | Rs] ^ Es, State) :-
    if_(enroll(R, R1)
	, State = [Q | Ls] ^ [R1 | Rs] ^ Es
	, ( length_plus_1(Ls, MTD),
	    State = declare_mtd(MTD)
	  )
       ).

deescalate([] ^ _ ^ _, declare_mtd(0)). % deescalate from already-lowest dose

deescalate([L | Ls] ^ Rs ^ Es, State) :-
    if_(enroll(L, L1)
	, State = Ls ^ [L1 | Rs] ^ Es
	, ( length_plus_1(Ls, MTD),
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

ccd_actions(_, declare_mtd(_)) --> [].
ccd_actions(CCD, S0) --> [(S0->A->S)],
			 { ccd_state0_action_state(CCD, S0, A, S) },
			 ccd_actions(CCD, S).

default_ccd(ccd([3/5, 4/8, 5/10, 6/12],		
		[1/1, 2/4, 3/6, 4/8, 5/11],
		[0/1, 1/8],
		12)).

ccd_d_path(CCD, D, Path) :-
    length(Tallies, D), maplist(=(0/0), Tallies),
    phrase(ccd_actions(CCD, []^Tallies^[]), Path).

%?- J+\(default_ccd(CCD), D=1, time(findall(M, ccd_d_path(CCD, D, P), Ps)), length(Ps, J)).
%@    % CPU time: 9.260 seconds
%@    % CPU time: 9.264 seconds
%@    J = 20. % ^ Showing 16x slowdown from dropping fast arithmetic branch of qcompare/4!
%@    % CPU time: 1.194 seconds
%@    % CPU time: 1.198 seconds
%@    J = 20. % ^ Showing the > 2x slowdown from reified CohortFull in enroll/3
%@    % CPU time: 0.565 seconds
%@    % CPU time: 0.569 seconds
%@    J = 20. % ^ PURE BASELINE
%@    % CPU time: 0.185 seconds
%@    % CPU time: 0.189 seconds
%@    J = 20.

%?- J+\(default_ccd(CCD), D=2, time(findall(M, ccd_d_path(CCD, D, P), Ps)), length(Ps, J)).
%@    % CPU time: 98.820 seconds
%@    % CPU time: 98.824 seconds
%@    J = 212. % ^ Showing 17x slowdown from dropping fast arithmetic branch of qcompare/4!
%@    % CPU time: 13.242 seconds
%@    % CPU time: 13.246 seconds
%@    J = 212. % ^ Showing 2.2x slowdown from reified CohortFull in enroll/3
%@    % CPU time: 5.777 seconds
%@    % CPU time: 5.781 seconds
%@    J = 212. % ^ PURE BASELINE
%@    % CPU time: 1.887 seconds
%@    % CPU time: 1.891 seconds
%@    J = 212.

%?- J+\(default_ccd(CCD), D=3, time(findall(M, ccd_d_path(CCD, D, P), Ps)), length(Ps, J)).
%@    % CPU time: 554.185 seconds
%@    % CPU time: 554.189 seconds
%@    J = 1151. % ^ Showing 17x slowdown from dropping fast arithmetic branch of qcompare/4!
%@    % CPU time: 31.437 seconds
%@    % CPU time: 31.441 seconds
%@    J = 1151.
%@    % CPU time: 71.946 seconds
%@    % CPU time: 71.950 seconds
%@    J = 1151. % ^ Showing 2.3x slowdown from reified CohortFull in enroll/3
%@    % CPU time: 31.445 seconds
%@    % CPU time: 31.449 seconds
%@    J = 1151. % ^ PURE BASELINE
%@    % CPU time: 10.430 seconds
%@    % CPU time: 10.434 seconds
%@    J = 1151.

%?- J+\(default_ccd(CCD), D=4, time(findall(M, ccd_d_path(CCD, D, P), Ps)), length(Ps, J)).
%@    % CPU time: 184.071 seconds
%@    % CPU time: 184.076 seconds
%@    J = 6718. % ^ PURE BASELINE
%@    % CPU time: 61.984 seconds
%@    % CPU time: 61.989 seconds
%@    J = 6718.
