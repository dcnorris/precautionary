% Generalizing the operational semantics of 'aliquots' for rolling enrollment
:- use_module(library(clpz)).
:- use_module(library(pio)).
:- use_module(library(lists)).
:- use_module(library(dcgs)).

% Prefix op * for 'generalizing away' goals (https://www.metalevel.at/prolog/debugging)
:- op(920, fy, *).
*_.  % *Goal always succeeds

/*

Under a scheme of rolling enrollment, we no longer have the alternating sequence
of enrollment and assessment 'coroutines'. Rather, the GENERAL state is one where
EITHER further enrollment or resolution of pending DLT assessments may occur.

What if I plunged in, and explicitly used Prolog's distinction between ground and
uninstantiated variables? The sequence of enrolled patients is then a list of vars,
each in domain 0..1. The state transitions consist of either the instantiation of
non-ground vars, or else the extension of the list with new enrollment. One happy
consequence of this design is that it looks open to generalization to ordinal tox.

*/

/*

That looks promising! So the state representation now becomes the previous
left-right lists, except that the element of these lists are themselves lists!

*/

%% A 'dose cohort' (or 'cohort' for short) is a list of toxicity assessments
%% which are var upon enrollment, and become ground when completed.
%% (The association with a dose will be effected through the L^R stack pair
%% as previously implemented in the 'aliquots' design.)

cohort --> [].
cohort --> [T], { T in 0..1 }, cohort.

%% We will judge a cohort according to its *tally*, a proper fraction:
cohort_tally(C, T/N) :-
    length(C, N),
    phrase(cohort, C), %% TODO: Just say C ins 0..1
    sum(C, #=, T). % clpz:sum/3 avoids awful "maplist(indomain, C), sum_list(C, T)"!

%?- phrase(cohort, C).
%@    C = []
%@ ;  C = [_A], clpz:(_A in 0..1)
%@ ;  C = [_A,_B], clpz:(_A in 0..1), clpz:(_B in 0..1)
%@ ;  C = [_A,_B,_C], clpz:(_A in 0..1), clpz:(_B in 0..1), clpz:(_C in 0..1)
%@ ;  C = [_A,_B,_C,_D], clpz:(_A in 0..1), clpz:(_B in 0..1), clpz:(_C in 0..1), clpz:(_D in 0..1)
%@ ;  ...

%?- cohort_tally(C, T/N).
%@    C = [], T = 0, N = 0
%@ ;  C = [T], N = 1, clpz:(T in 0..1)
%@ ;  C = [_A,_B], N = 2, clpz:(_A+_B#=T), clpz:(T in 0..2), clpz:(_A in 0..1), clpz:(_B in 0..1)
%@ ;  C = [_A,_B,_C], N = 3, clpz:(_A+_B+_C#=T), clpz:(_A in 0..1), clpz:(_B in 0..1), clpz:(_C in 0..1), clpz:(T in 0..3)
%@ ;  C = [_A,_B,_C,_D], N = 4, clpz:(_A+_B+_C+_D#=T), clpz:(_A in 0..1), clpz:(_B in 0..1), clpz:(_C in 0..1), clpz:(_D in 0..1), clpz:(T in 0..4)
%@ ;  C = [_A,_B,_C,_D,_E], N = 5, clpz:(_A+_B+_C+_D+_E#=T), clpz:(_A in 0..1), clpz:(_B in 0..1), clpz:(_C in 0..1), clpz:(_D in 0..1), clpz:(_E in 0..1), clpz:(T in 0..5)
%@ ;  C = [_A,_B,_C,_D,_E,_F], N = 6, clpz:(_A+_B+_C+_D+_E+_F#=T), clpz:(_A in 0..1), clpz:(_B in 0..1), clpz:(_C in 0..1), clpz:(_D in 0..1), clpz:(_E in 0..1), clpz:(_F in 0..1), clpz:(T in 0..6)
%@ ;  ...

%% Let's also be able to obtain best- and worst-case tallies
%% from a partially complete cohort of DLT assessments.
cohort_bestcase([], []).
cohort_bestcase([T|Ts], [B|Bs]) :-
    (	var(T),
	B #= 0
    ;	ground(T),
	B = T
    ),
    cohort_bestcase(Ts, Bs).

cohort_worstcase([], []).
cohort_worstcase([T|Ts], [W|Ws]) :-
    (	var(T),
	W #= 1
    ;	ground(T),
	W = T
    ),
    cohort_worstcase(Ts, Ws).

%?- cohort_bestcase([0,1,T], B).
%@    B = [0,1,0]
%@ ;  false.

%?- cohort_worstcase([0,1,T], W).
%@    W = [0,1,1]
%@ ;  false.

cohort_mintally(C, T/N) :-
    cohort_bestcase(C, B),
    cohort_tally(B, T/N).

cohort_maxtally(C, T/N) :-
    cohort_worstcase(C, W),
    cohort_tally(W, T/N).

%?- cohort_mintally([0,1,T], Q).
%@    Q = 1/3
%@ ;  false.

%?- cohort_maxtally([0,1,T], Q).
%@    Q = 2/3
%@ ;  false.

%% A fair enumeration of tallies will be helpful for developing
%% and testing comparison relations between tallies.

tally(DLTs/Enrolled) :-
    Enrolled in 0..6, % suffices for 3+3 and 'rolling 6' designs
    indomain(Enrolled),
    DLTs in 0..Enrolled,
    %%0 #=< DLTs, DLTs #=< Enrolled, 
    indomain(DLTs).

%% Proof that tally/1 terminates!
%?- tally(_), false.
%@ false.

%?- tally(T/N).
%@    T = 0, N = 0
%@ ;  T = 0, N = 1
%% ... as expected ...
%@ ;  T = 4, N = 6
%@ ;  T = 5, N = 6
%@ ;  T = 6, N = 6.

%% How do tallies COMPARE?
safer_than(T0/N0, T1/N1) :-
    tally(T0/N0), N0 #> 0, % NB: These N>0 conditions are implicit in the #< below,
    tally(T1/N1), N1 #> 0, %     but are stated explicitly here to expose the logic.
    N0 #>= N1,
    T0*N1 #< N0*T1. % 'cross-multiply'

:- op(900, xfx, &<).
&<(Q1, Q2) :- safer_than(Q1, Q2).
%?- Q &< 1/3.
%@    Q = 0/3
%@ ;  Q = 0/4
%@ ;  Q = 1/4
%@ ;  Q = 0/5
%@ ;  Q = 1/5
%@ ;  Q = 0/6
%@ ;  Q = 1/6
%@ ;  false.

noworse_than(T0/N0, T1/N1) :-
    tally(T0/N0), N0 #> 0, % NB: By contrast with the N>0 conditions in safer_than/2,
    tally(T1/N1), N1 #> 0, %     these are strictly necessary for correctness here.
    (	N0 #>= N1,
	T0*N1 #=< N0*T1 % 'cross-multiply'
    ;	N0 #< N1,
	T0 + (N1 - N0) #=< T1 % extend T0/N0 under worst-case scenario
    ).

:- op(900, xfx, &=<).
&=<(Q1, Q2) :- noworse_than(Q1, Q2).
%?- 1/3 &=< Q.
%@    Q = 1/1
%@ ;  Q = 1/2
%@ ;  Q = 2/2
%@ ;  Q = 1/3
%@ ;  Q = 2/3
%@ ;  Q = 3/3
%@ ;  Q = 2/4
%@ ;  Q = 3/4
%@ ;  Q = 4/4
%@ ;  Q = 3/5
%@ ;  Q = 4/5
%@ ;  Q = 5/5
%@ ;  Q = 4/6
%@ ;  Q = 5/6
%@ ;  Q = 6/6.

on_par(T0/N0, T1/N1) :-
    tally(T0/N0), N0 #> 0,
    tally(T1/N1), N1 #> 0,
    T0*N1 #= N0*T1.

:- op(900, xfx, &=).
&=(Q1, Q2) :- on_par(Q1, Q2).

% Strict tally-safer (&<) excludes tally-equals (&=)
%?- Q &< R, Q &= R.
%@ false.

% (&=<) and '&>' are MUTUALLY EXCLUSIVE:
%?- Q &=< R, R &< Q.
%@ false.

%% Here was my motivation for writing down on_par/2 in the first place:
%?- Q &= 1/3. %% 1/3 is (roughly) the 3+3 design's 'target toxicity rate'
%@    Q = 1/3
%@ ;  Q = 2/6
%@ ;  false.

%% TODO: Consider deriving all tally comparisons from cross-product zcompare,
%%       taking care to exclude 'incomparables' somehow.

%?- cohort_tally([0,0,T], Q), safer_than(Q, 2/3).
%@    Q = 0/3, T = 0
%@ ;  Q = 1/3, T = 1
%@ ;  false.

%% Uh-oh!
%?- Vs = [A,B,C], Vs ins 0..1, sum(Vs, #=, Sum), labeling([min(Sum)], Vs).
%@    Vs = [0,A,A], A = 0, B = 0, Sum = 0, C = 0
%@ ;  Vs = [0,A,1], A = 0, B = 0, Sum = 1, C = 1
%@ ;  Vs = [0,1,A], A = 0, B = 1, Sum = 1, C = 0
%@ ;  Vs = [1,0,B], A = 1, B = 0, Sum = 1, C = 0
%@ ;  Vs = [0,1,B], A = 0, B = 1, Sum = 2, C = 1
%@ ;  Vs = [1,0,A], A = 1, B = 0, Sum = 2, C = 1
%@ ;  Vs = [1,A,0], A = 1, B = 1, Sum = 2, C = 0
%@ ;  ...

%?- X = a, nonvar(X).
%@    X = a. 

%?- nonvar(X).
%@ false.

%% Therefore var/nonvar makes no sense declaratively!
%% var/1 and nonvar/1 are metalogical

%% So instead we can do this:
% Vs = [A,B,C], Vs ins 0..1, sum(Vs, #=, Sum), labeling([min(Sum)], Vs).
%% NB: If I want the first minimum, wrap this in once/1.


%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Pick up from here...

%% Opportunity to replace the hard-coding here with meta-interpretation.
%% You might interpret 'more fairly', e.g.

%% ENROLL new patient 'T' into the lowest right-hand dose
state0_action_state(Ls ^ [R | Rs], enroll, Ls ^ [[T|R] | Rs]) :-
    enrollable_cohort(R),
    T in 0..1, % var T is pending tox assessment of the newly-enrolled patient
    %% TODO: impose restriction on acceptable tally
    true.

enrollable_cohort(R) :-
    length(R, N), N #< 6, % limit max cohort size
    true.


% Which tallies do not rule out further enrollment at a dose?
enrollable_tally(T0/N0) :-
    tally(T0/N0),
    (	N0 #= 0
    ;	N0 #= 3, T0 #=< 1
    ).
%?- enrollable_tally(Q).
%@ Q = 0/0 ;
%@ Q = _3234/3,
%@ _3234 in 0..1 ;
%@ false.

% Which tallies RULE OUT further enrollment at a dose?
unenrollable_tally(T/3) :- tally(T/3), T #> 1. % 'too toxic'
unenrollable_tally(T/6) :- tally(T/6). % 'enough already'

%?- state0_action_state(S0, enroll, S).
%@ S0 = _7610:[0/0|_7618],        % We enroll a new dose 0/0,
%@ S = _7610^[_7652/3|_7618],
%@ _7652 in 0..3 ;
%@ S0 = _11382:[_11394/3|_11390], % or an already-tried dose ..
%@ S = _11382^[_11424/6|_11390],
%@ _11394 in 0..1,                % of 0/3 or 1/3
%@ _11394+_11474#=_11424,         % and add correctly!
%@ _11474 in 0..3,                % This is # of new toxicities,
%@ _11424 in 0..4 ;               % which can sum to 4/6 at most.
%@ false.

/*
Instead of regarding merely TALLIES as presumably safe/toxic,
we apply these JUDGMENTS rather (and more properly) to the DOSES
themselves -- as represented through lists.
*/
presumably_safe([]). % [] corresponds to zero dose
presumably_safe([Q|_]) :- % "no need to enroll further at this dose to say it's okay"
    tally(Q),
    Q = T/6,   % TODO: Might I obtain the no-deescalation variant of 3+3
    T in 0..1. %       simply by modifying these 2 goals?

presumably_toxic([]) :- false.
presumably_toxic([Q|_]) :-
    tally(Q),
    Q = T/_,
    T #> 1.

%?- presumably_toxic([Q|_]).
%@ Q = _5352/3,
%@ _5352 in 2..3 ;
%@ Q = _6862/6,
%@ _6862 in 2..6.

%% PROOF: Presumable tallies are NOT enrollable
%?- (presumably_safe([T]); presumably_toxic([T])), enrollable_tally(T).
%@ false.

%% PROOF: Enrollable and unenrollable tallies are disjoint
%?- enrollable_tally(Q), unenrollable_tally(Q).
%@ false.

%?- enrollable_tally(T/N).
%@ T = N, N = 0 ;
%@ N = 3,
%@ T in 0..1 ;
%@ false.

% Handling the (exceptional) MTD-not-found case in the ENROLLMENT
% phase is formally quite natural, even if it appears odd from a
% perspective that ':' is the state upon which the 'routine' trial
% operations act. (The literal interpretation of ':' as the phase
% where trial administrative staff are at work is thus misleading.)
%% AHA! This demands modification along with the 'conspicuous'
%% clause below.
%% state0_action_state(Ls : [], stop, mtd_notfound(D)) :-
%%     presumably_safe(Ls),
%%     length(Ls, D). % RP2D is highest prespecified dose, D.
% This is by far the MOST CONSPICUOUS clause of state0_action_state/3,
% since it clearly catches an edge case _:[] and changes the stacks.
% There is a real sense in which this exceptional case is handled as
% an 'oopsie' here, rather than straightforwardly.
% Perhaps we see here a tension between 'elegance' and 'bluntness',
% BOTH of which seem like attrbutes of good Prolog programming.
% ---
% Now that I find this clause causing problems during DCG translation
% to the compact D{^,-,*}T notation, I believe a special term is
% warranted to distinguish this case.
%% state0_action_state([Q0|Ls] : [], enroll, Ls ^ [Q]) :-
%%     tally0_cohort_tally(Q0, _, Q).


%% JUDGE the toxicities, and WEIGH the dose-escalation decision.
%% NB: We require a dose R that we can escalate TO!
state0_action_state(Ls ^ [Q,R|Rs], escalate, [Q|Ls] : [R|Rs]) :-
    Q &< 1/3.
% In the case where we WOULD ESCALATE, BUT CAN'T, we record this
% as a 'clamp' action, by analogy with voltage clamping.
state0_action_state(Ls ^ [Q], clamp, Ls : [Q]) :-
    Q &< 1/3,
    enrollable_tally(Q).
state0_action_state(Ls^[Q|Rs], stay, Ls:[Q|Rs]) :-
    Q &= 1/3, % <-- Roughly a statement about meeting 'target toxicity rate'
    enrollable_tally(Q).
state0_action_state([L|Ls] ^ Rs, deescalate, Ls : [L]) :-
    enrollable_tally(L),
    presumably_toxic(Rs).
state0_action_state(Ls ^ [R|Rs], stop, declare_mtd(MTD)) :-
    presumably_safe(Ls),
    presumably_toxic([R|Rs]),
    length(Ls, MTD).
% Here again we treat specially the case where 1st arg is Ls ^ [Q] ...
state0_action_state(Ls ^ [Q], stop, declare_mtd(MTD)) :-
    presumably_safe([Q|Ls]),
    length([Q|Ls], MTD).


% What does a  cohort look like?
cohort(DLTs/N) :-
    N #= 3, % initially, we permit only cohorts of 3
    DLTs in 0..N,
    *indomain(DLTs).
%?- cohort(C).
%@ C = _44170/3,
%@ _44170 in 0..3.


%?- unenrollable_tally(Q).
%@ Q = _5464/3,
%@ _5464 in 2..sup ;
%@ Q = _6286/6.

% How do cohorts ACCUMULATE into tallies?
tally0_cohort_tally(T0/N0, T_/N_, T/N) :-
    enrollable_tally(T0/N0),
    cohort(T_/N_),
    tally(T/N),
    T #= T0 + T_,
    N #= N0 + N_.
%?- tally0_cohort_tally(T0, C, T).
%@ T0 = 0/0,
%@ C = T, T = _10336/3,
%@ _10336 in 0..3 ;
%@ T0 = _14210/3,
%@ C = _14228/3,
%@ T = _14246/6,
%@ _14210 in 0..1,
%@ _14210+_14228#=_14246,
%@ _14228 in 0..3,
%@ _14246 in 0..4 ;
%@ false.


/* A truly ESSENTIAL core pharmacologic concept undergirding dose-finding
 * trials is the monotonicity of the dose-toxicity function. To manifest
 * this notion in our code, we model the DOSE as a DESCENDING LIST.
 * (Cf. the treatment of the hyperreals as convergent sequences of reals.)
 * In this treatment, it is not necessary to carry along dose labels in
 * these lists. The mere STRUCTURE of the list suffices to enable integer
 * doses to be 'read off' via length/2. Accordingly, we are free to use
 * the list elements to maintain tallies of toxicities observed at the
 * several doses of the trial.
*/

%% NB: You can do this instead, since tally/1 terminates:
%% tallies(Qs) :- maplist(tally, Qs).

% Describe a LIST of TALLIES for consecutive prespecified doses:
tallies([]).
tallies([Q|Qs]) :-
    length(Qs, _), % enumerate solutions fairly
    tally(Q),
    tallies(Qs).
%?- length(C, 1), tallies(C).
%@ C = [0/0] ;
%@ C = [_7170/3],
%@ _7170 in 0..3 ;
%@ C = [_8490/6],
%@ _8490 in 0..6. % NB: Only 0..3 could occur in a 3+3 trial

%?- length(C, 2), tallies(C), false. % Does it terminate?
%@ false.                            % Yes.

/* An assumption we make about the dose-escalation trials modeled here
 * is that AT ANY GIVEN MOMENT they always partition the available doses
 * into 2 groups:
 * 'L' - A descending list of doses considered 'too Low' to enroll RIGHT NOW
 * 'R' - An ascending list of doses of which the head looks Right to enroll.
 * We can imagine these lists as 2 stacks, on the 'left' and 'right'.
 * 
 * Since it proves useful to describe the trial in two alternating phases,
 * we designate the state of the trial via L and R lists conjoined with
 * alternating infix functors:
 * '^' - Imagine: 'scales of Justice' where we make judgements about tallies
 * ':' - Imagine: an ellipsis (or dripping faucet) where we enroll patients
 *       (who 'trickle in'), wait for & assess toxicities.
 */



/****
The following general queries seem to PROVE
key properties of 3+3 have been attained.
 ****/
%?- RP2D in 0..sup, state0_action_state(S0, stop, mtd_notfound(RP2D)).
%@ RP2D = 0,
%@ S0 = []:[] ;
%@ RP2D = 1,
%@ S0 = [_6146/6]:[],
%@ _6146 in 0..1 ;
%@ RP2D = 2,
%@ S0 = [_7652/6, _7658]:[],
%@ _7652 in 0..1 ;
%@ RP2D = 3,
%@ S0 = [_9170/6, _9176, _9182]:[],
%@ _9170 in 0..1 ; ...

%?- MTD in 0..sup, state0_action_state(S0, stop, declare_mtd(MTD)).
%@ MTD = 0,
%@ S0 = []^[_2400/3|_2396],
%@ _2400 in 2..3 ;
%@ MTD = 0,
%@ S0 = []^[_4346/6|_4342],
%@ _4346 in 2..6 ;
%@ MTD = 1,
%@ S0 = [_7710/6]^[_7722/3|_7718],
%@ _7710 in 0..1,
%@ _7722 in 2..3 ;
%@ MTD = 2,
%@ S0 = [_9476/6, _9482]^[_9494/3|_9490],
%@ _9476 in 0..1,
%@ _9494 in 2..3 ;
%@ MTD = 3,
%@ S0 = [_11254/6, _11260, _11266]^[_11278/3|_11274],
%@ _11254 in 0..1,
%@ _11278 in 2..3 ; ...

%?- MTD in 0..3, indomain(MTD), state0_action_state(S0, stop, declare_mtd(MTD)).
%@ MTD = 0,
%@ S0 = []^[_2842/3|_2838],
%@ _2842 in 2..3 ;
%@ MTD = 0,
%@ S0 = []^[_4762/6|_4758],
%@ _4762 in 2..6 ;
%@ MTD = 1,
%@ S0 = [_8222/6]^[_8234/3|_8230],
%@ _8222 in 0..1,
%@ _8234 in 2..3 ;
%@ MTD = 1,
%@ S0 = [_10396/6]^[_10408/6|_10404],
%@ _10396 in 0..1,
%@ _10408 in 2..6 ;
%@ MTD = 2,
%@ S0 = [_14060/6, _14066]^[_14078/3|_14074],
%@ _14060 in 0..1,
%@ _14078 in 2..3 ;
%@ MTD = 2,
%@ S0 = [_16252/6, _16258]^[_16270/6|_16266],
%@ _16252 in 0..1,
%@ _16270 in 2..6 ;
%@ MTD = 3,
%@ S0 = [_19894/6, _19900, _19906]^[_19918/3|_19914],
%@ _19894 in 0..1,
%@ _19918 in 2..3 ;
%@ MTD = 3,
%@ S0 = [_22104/6, _22110, _22116]^[_22128/6|_22124],
%@ _22104 in 0..1,
%@ _22128 in 2..6.
/* ---
 This makes a very nice characterization of the 'typical' or 'motivating'
 case of successful conclusion of the trial. In this variant of 3+3, the
 most desirable conclusion is to have observed {0,1}/6 toxities at the
 declared MTD, and excessive toxicity (2 or more) at the next-higher dose.
 */

%?- state0_action_state(S0, stop, mtd_notfound(MTD)).
%@ S0 = []:[],
%@ MTD = 0 ;
%@ S0 = [_1510]:[],
%@ MTD = 1 ;
%@ S0 = [_1510, _2558]:[],
%@ MTD = 2 ;
%@ S0 = [_1510, _2558, _3606]:[],
%@ MTD = 3 .

%% DISCUSS: Non-terminating. From a verification perspective,
%% this is not so troubling, since it merely asks Prolog to
%% confirm something we can already verify BY INSPECTION.
%?- state0_action_state(_:[], A, _), dif(A, stop).
%@ Action (h for help) ? abort
%@ % Execution Aborted

%?- state0_action_state(S0, deescalate, S).
%@ S0 = [0/3|_9574]^[[_9596/3|_9592]|_9586],
%@ S = _9574:[0/3],
%@ _9596 in 2..3 ;
%@ S0 = [0/3|_11542]^[[_11564/6|_11560]|_11554],
%@ S = _11542:[0/3],
%@ _11564 in 2..6.
/* ---
 What does this show about de-escalation?
 1. It only ever occurs back to a dose where tally is 0/3
 2. It never occurs after any of {0,1}/3, {0,1}/6
 3. Resulting state S is of 'pending' type ":"
 4. This 0/3 tally is the 'focus dose' at head of right-hand list
 5. The dose we de-escalated FROM looked like T/6 with T in 3..6
 6. The dose we de-escalated FROM gets TRUNCATED from right-hand list
 */

%?- state0_action_state(S0, escalate, S).
%@ S0 = _504^[0/3|_512],
%@ S = [0/3|_504]:_512 ;
%@ S0 = _4778^[_4790/6|_4786],
%@ S = [_4790/6|_4778]:_4786,
%@ _4790 in 0..1 ;
%@ false.
/* ---
 What does this say about escalation?
 1. It may occur ONLY after a 0/3, 0/6 or 1/6 observation
    (NB: The 0/6 observation is not excluded BY THIS PREDICATE ALONE.)
 2. Escalation leaves us in a ':' state
 3. The RHS of the ':' state MAY be empty list [].
 */

%?- state0_action_state(S0, stay, S).
%@ S0 = _7196^[1/3|_7204],
%@ S = _7196:[1/3|_7204] ;
%@ false.
/* ---
 That result speaks for itself!
 */


%?- X is 4 rdiv 7.
%@ X = 4r7.

%?- X is log(3).
%@ X = 1.0986122886681098.
%@ ERROR: Syntax error: Operator expected
%@ ERROR: X is log
%@ ERROR: ** here **
%@ ERROR:  3 . 

%% About which (indeterminate) tallies do we presume nothing?
%?- tally(I), \+ presumably_safe([I]), \+ presumably_toxic([I]).
%@ I = 0/0 ;
%@ I = 0/3 ;
%@ I = 1/3 ;
%@ false.
%% NOTE: These are precisely the enrollable doses!

%?- enrollable_tally(T).
%@ T = 0/0 ;
%@ T = _8562/3,
%@ _8562 in 0..1 ;
%@ false.

