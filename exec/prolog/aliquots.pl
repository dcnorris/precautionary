% Attempt an enumeration of ALL 3+3 trials
:- use_module(library(clpfd)).
:- use_module(library(pio)).

% Prefix op * for 'generalizing away' goals (https://www.metalevel.at/prolog/debugging)
:- op(920, fy, *). *_.  % *Goal always succeeds

% What does a  cohort look like?
cohort(DLTs/N) :-
    N #= 3, % initially, we permit only cohorts of 3
    DLTs in 0..N,
    *indomain(DLTs).
%?- cohort(C).
%@ C = _44170/3,
%@ _44170 in 0..3.

% What does a DLT tally look like?
tally(DLTs/Enrolled) :-
    Enrolled in 0..6,
    *indomain(Enrolled),
    % -------------------------------------------------------------------------------
    Ncohorts in 0..2,         % TODO: Release this constraint to permit finer-grained
    indomain(Ncohorts),       %       enrollment 0..6 as above, which will be needed
    Enrolled #= Ncohorts * 3, %       to model Simon &al 1997 accelerated titration.
    % -------------------------------------------------------------------------------
    DLTs in 0..Enrolled,
    *indomain(DLTs).
%?- tally(C).
%@ C = 0/0 ;
%@ C = _9656/3,
%@ _9656 in 0..3 ;
%@ C = _10970/6,    % NB: This allows certain tallies that cannot actually
%@ _10970 in 0..6.  %     arise during a 3+3 trial: 3/6, 4/6, 5/6, 6/6.

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

% How do tallies COMPARE?
safer_than(T0/N0, T1/N1) :-
    tally(T0/N0), N0 #> 0, % NB: These N>0 conditions are implicit in the #< below,
    tally(T1/N1), N1 #> 0, %     but are stated explicitly here to expose the logic.
    T0*N1 #< N0*T1. % 'cross-multiply'

:- op(900, xfx, user:(&<)).
&<(Q1, Q2) :- safer_than(Q1, Q2).
%?- Q &< 1/3.
%@ Q = 0/3 ;
%@ Q = _9866/6,
%@ _9866 in 0..1 ;
%@ false.

noworse_than(T0/N0, T1/N1) :-
    tally(T0/N0), N0 #> 0, % NB: By contrast with the N>0 conditions in safer_than/2,
    tally(T1/N1), N1 #> 0, %     these are strictly necessary for correctness here
    T0*N1 #=< N0*T1.       % <-- because the weaker comparison '#=<' is used.

:- op(900, xfx, user:(&=<)).
&=<(Q1, Q2) :- noworse_than(Q1, Q2).
%?- 1/3 &=< Q.
%@ Q = _2456/3,
%@ _2456 in 1..3 ;
%@ Q = _4352/6,
%@ _4352 in 2..6 ;
%@ false.

on_par(T0/N0, T1/N1) :-
    tally(T0/N0), N0 #> 0,
    tally(T1/N1), N1 #> 0,
    T0*N1 #= N0*T1.

:- op(900, xfx, user:(&=)).
&=(Q1, Q2) :- on_par(Q1, Q2).

%?- Q1 &= Q2. %% DISCUSS: Duplicated equivalent solutions:
%@ Q1 = Q2, Q2 = _10596/3,
%@ _10596 in 0..3 ;
% ...
%@ Q1 = Q2, Q2 = _25876/3,
%@ _25876 in 0..3 ;
% ...

% Strict tally-safer (&<) excludes tally-equals (&=)
%?- Q &< R, Q &= R.
%@ false.

% (&=<) and '&>' are MUTUALLY EXCLUSIVE:
%?- Q &=< R, R &< Q.
%@ false.

%% Here was my motivation for writing down on_par/2 in the first place:
%?- Q &= 1/3. %% 1/3 is (roughly) the 3+3 design's 'target toxicity rate'
%@ Q = 1/3 ;
%@ Q = 2/6 ;
%@ false.

%% TODO: Consider deriving all tally comparisons from cross-product zcompare,
%%       taking care to exclude 'incomparables' somehow.

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


:- discontiguous state0_action_state/3.
% NOTE: The currently quite trivial ENROLLMENT PHASE may well grow
%       more complex with introduction of an accelerated titration
%       phase or other such variations on standard 3+3.

%% ENROLL into the ':' state
state0_action_state(Ls : [Q0 | Rs], enroll, Ls ^ [Q | Rs]) :-
    tally0_cohort_tally(Q0, _, Q).
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

%%actions(_) --> [].
actions(mtd_notfound(_)) --> [].
actions(declare_mtd(_)) --> [].
actions(S0) --> [S0->A->S],
        { state0_action_state(S0, A, S) },
        actions(S).

%% TODO: Separate logical properties of what I write out
%%       from the actual writing-out of that output.
%%       Ideal description is always formulated by pure means.

%% Let me BEGIN with a simple CONDENSATION of the S0->A-S terms.
%% NB: This will help to DEFINE the condensed notation *itself*,
%%     in terms of the domain-level concepts elaborated through
%%     state0_action_state/3 and actions//1.
%%
%% Might I do this best using SEMICONTEXT NOTATION? I wonder if
%% such a 'translation' effectively DEMONSTRATES "how to read"
%% or "how to *interpret*" the events and course of a trial.
%% On that view, you might even say "the more DCGs, the better",
%% since DIVERSE 'interpretations' lend DEPTH and FALSIFIABILITY
%% to the CONTENT of the program. Every new interpretation of
%% the dose-escalation design generates new ways to express our
%% expectations and understandings OBJECTIVELY.
%%

%% TODO: Make this DCG translator more abstract
%%       by implementing a 'state-difference'
%%       operator that yields enrolled cohort.
%% NB: This will also reduce duplicated code!

state0_cohort_state_dose(L:[Q0|R], T/N, L^[Q1|R], D) :-
    tally0_cohort_tally(Q0, T/N, Q1),
    length([Q1|L], D),
    indomain(T), indomain(N),
    *indomain(D).

condensed, [D^T] -->
    [ _->escalate->_, S0->enroll->S ],
    { state0_cohort_state_dose(S0, T/_, S, D) },
    condensed.
condensed, [D^T] --> % handle special case of *start* of trial
    [ S0->enroll->S ],
    { state0_cohort_state_dose(S0, T/_, S, D) },
    condensed.
condensed, [D-T] -->
    [ _->stay->_, S0->enroll->S ],
    { state0_cohort_state_dose(S0, T/_, S, D) },
    condensed.
condensed, [D:T] -->
    [ _->deescalate->_, S0->enroll->S ],
    { state0_cohort_state_dose(S0, T/_, S, D) },
    condensed.
condensed, [D*T] -->
    [ _->clamp->_, S0->enroll->S ],
    { state0_cohort_state_dose(S0, T/_, S, D) },
    condensed.
condensed, [declare_mtd(MTD)] --> [ _->stop->declare_mtd(MTD) ].
condensed, [mtd_notfound(MTD)] --> [ _->stop->mtd_notfound(MTD) ].
%condensed --> []. %% Uncomment this for 'catch-all' permitting partial translations.

%% Example of set membership constraint.
%?- tuples_in([[A,B,C],[D,E,F]], [[1,2,3],[4,5,6]]), indomain(A).
%@ A = 1,
%@ B = 2,
%@ C = 3,
%@ D in 1\/4,
%@ tuples_in([[D, E, F]], [[1, 2, 3], [4, 5, 6]]),
%@ E in 2\/5,
%@ F in 3\/6 ;
%@ A = 4,
%@ B = 5,
%@ C = 6,
%@ D in 1\/4,
%@ tuples_in([[D, E, F]], [[1, 2, 3], [4, 5, 6]]),
%@ E in 2\/5,
%@ F in 3\/6.

%?- clause(condensed(Ls0,Ls), Body).
%@ Body =  (Ls0=[(_47356->escalate->_47364),  (_47374->enroll->_47382)|_47370], state0_cohort_state_dose(_47374, _47392/_47394, _47382, _47404), _47408=_47370, condensed(_47408, _47416), Ls=[_47404^_47392|_47416]) ;
%@ Body =  (Ls0=[(_48280->enroll->_48288)|_48276], state0_cohort_state_dose(_48280, _48298/_48300, _48288, _48310), _48314=_48276, condensed(_48314, _48322), Ls=[_48310^_48298|_48322]) ;
%@ Body =  (Ls0=[(_49186->stay->_49194),  (_49204->enroll->_49212)|_49200], state0_cohort_state_dose(_49204, _49222/_49224, _49212, _49234), _49238=_49200, condensed(_49238, _49246), Ls=[_49234-_49222|_49246]) ;
%@ Body =  (Ls0=[(_50110->deescalate->_50118),  (_50128->enroll->_50136)|_50124], state0_cohort_state_dose(_50128, _50146/_50148, _50136, _50158), _50162=_50124, condensed(_50162, _50170), Ls=[_50158:_50146|_50170]) ;
%@ Body =  (Ls0=[(_51034->clamp->_51042),  (_51052->enroll->_51060)|_51048], state0_cohort_state_dose(_51052, _51070/_51072, _51060, _51082), _51086=_51048, condensed(_51086, _51094), Ls=[_51082*_51070|_51094]) ;
%@ Ls0 = [(_51952->stop->declare_mtd(_51964))|_51948],
%@ Ls = [declare_mtd(_51964)|_51948],
%@ Body = true ;
%@ Ls0 = [(_53394->stop->mtd_notfound(_53406))|_53390],
%@ Ls = [mtd_notfound(_53406)|_53390],
%@ Body = true.

%% Examine the smallest possible trial -- a trial with just 1 dose!
%?- phrase(actions([]:[0/0]), Trial), phrase(condensed, Trial, Translation).
%@ Trial = [([]:[0/0]->enroll->[]^[0/3]),  ([]^[0/3]->clamp->[]:[0/3]),  ([]:[0/3]->enroll->[]^[_16320/6]),  ([]^[_16320/6]->stop->declare_mtd(0))],
%@ Translation = [1^0, 1*_16320, declare_mtd(0)],
%@ _16320 in 2..3 ;
%@ Trial = [([]:[0/0]->enroll->[]^[0/3]),  ([]^[0/3]->clamp->[]:[0/3]),  ([]:[0/3]->enroll->[]^[_20634/6]),  ([]^[_20634/6]->stop->declare_mtd(1))],
%@ Translation = [1^0, 1*_20634, declare_mtd(1)],
%@ _20634 in 0..1 ;
%@ Trial = [([]:[0/0]->enroll->[]^[1/3]),  ([]^[1/3]->stay->[]:[1/3]),  ([]:[1/3]->enroll->[]^[_28502/6]),  ([]^[_28502/6]->stop->declare_mtd(0))],
%@ Translation = [1^1, 1-_28562, declare_mtd(0)],
%@ _28502 in 2..4,
%@ 1+_28562#=_28502,
%@ 1+_28644#=_28502,
%@ _28562 in 1..3,
%@ _28644 in 1..3 ;
%@ Trial = [([]:[0/0]->enroll->[]^[1/3]),  ([]^[1/3]->stay->[]:[1/3]),  ([]:[1/3]->enroll->[]^[1/6]),  ([]^[1/6]->stop->declare_mtd(1))],
%@ Translation = [1^1, 1-0, declare_mtd(1)] ;
%@ Trial = [([]:[0/0]->enroll->[]^[_5564/3]),  ([]^[_5564/3]->stop->declare_mtd(0))],
%@ Translation = [1^_5564, declare_mtd(0)],
%@ _5564 in 2..3 ;
%@ false.

%% NEXT: Translate paths from the 'condensed' notation into MATRICES.
%% Transform the A path representations to the arrays T(c,d,j).

path_matrix(D, P, M) :-
    length(C1,D),
    length(C2,D),
    M = (C1,C2),
    path_matrix_(P, M),
    maplist(ground_or_nil, C1),
    maplist(ground_or_nil, C2).

path_matrix_([D^T|Rest], (C1,C2)) :-
    nth1(D, C1, T),
    path_matrix_(Rest, (C1,C2)).

path_matrix_([D-T|Rest], (C1,C2)) :-
    nth1(D, C2, T),
    path_matrix_(Rest, (C1,C2)).

path_matrix_([D*T|Rest], (C1,C2)) :- path_matrix_([D-T|Rest], (C1,C2)).
path_matrix_([D:T|Rest], (C1,C2)) :- path_matrix_([D-T|Rest], (C1,C2)).

path_matrix_([declare_mtd(_)], _).
path_matrix_([mtd_notfound(_)], _).

path_matrix_([], ([],[])).

ground_or_nil(Term) :- ground(Term).
ground_or_nil('NA') :- true.

pathmatrix((First,Second)) -->
    row(First),
    row(Second).

% Christmas-tree predicate!
columns_format(1, '~w~t~3+').
columns_format(2, '~w~t~3+~w~t~3+').
columns_format(3, '~w~t~3+~w~t~3+~w~t~3+').
columns_format(4, '~w~t~3+~w~t~3+~w~t~3+~w~t~3+').

% Approximate format_//2 as provided by Scryer's library(format):
format_(Format, Ls) --> [ FLs ], { format(atom(FLs), Format, Ls) }.

row(Ls) -->
    { length(Ls,D),
      columns_format(D,Format) },
    format_(Format, Ls).

%?- phrase(row(X), Row).
%@ X = [_28598],
%@ Row = ['_28598'] .
%@ X = [_44750],
%@ Row = ['_44750'] ;
%@ X = [_44750, _45792],
%@ Row = ['_44750_45792'] ;
%@ X = [_44750, _45792, _46834],
%@ Row = ['_44750_45792_46834'] ;
%@ X = [_44750, _45792, _46834, _47876],
%@ Row = ['_44750_45792_46834_47876'] .

format_matrix(Matrix) :- format('~s~n~s~n', Matrix).

%?- phrase(pathmatrix(([1,'NA'],['NA','NA'])), PM), format_matrix(PM).
%@ 1   NA  
%@ NA  NA  
%@ PM = ['1   NA  ', 'NA  NA  '].

%?- phrase(row([this,is,a]), Row), format('~s', Row).
%@ thisis a  
%@ Row = ['thisis a  '].
%@ this is a

%?- phrase(actions([] : [0/0,0/0]), Trial), phrase(condensed, Trial, Translation), path_matrix(2, Translation, M), phrase(pathmatrix(M), Matrix), format_matrix(Matrix).
%@ 0  0  
%@ 2  2  
%@ Trial = [([]:[0/0, 0/0]->enroll->[]^[0/3, 0/0]),  ([]^[0/3, 0/0]->escalate->[0/3]:[0/0]),  ([0/3]:[0/0]->enroll->[0/3]^[0/3]),  ([0/3]^[0/3]->clamp->[0/3]:[0/3]),  ([0/3]:[0/3]->enroll->[... / ...]^[... / ...]),  ([... / ...]^[... / ...]->deescalate->[]:[...]),  ([]:[...]->enroll-> ... ^ ...),  (... ^ ... -> ... -> ...)],
%@ Translation = [1^0, 2^0, 2*2, 1:2, declare_mtd(0)],
%@ M =  ([0, 0], [2, 2]),
%@ Matrix = ['0  0  ', '2  2  '] ;
%@ 0  0  
%@ 2  3  
%@ Trial = [([]:[0/0, 0/0]->enroll->[]^[0/3, 0/0]),  ([]^[0/3, 0/0]->escalate->[0/3]:[0/0]),  ([0/3]:[0/0]->enroll->[0/3]^[0/3]),  ([0/3]^[0/3]->clamp->[0/3]:[0/3]),  ([0/3]:[0/3]->enroll->[... / ...]^[... / ...]),  ([... / ...]^[... / ...]->deescalate->[]:[...]),  ([]:[...]->enroll-> ... ^ ...),  (... ^ ... -> ... -> ...)],
%@ Translation = [1^0, 2^0, 2*3, 1:2, declare_mtd(0)],
%@ M =  ([0, 0], [2, 3]),
%@ Matrix = ['0  0  ', '2  3  '] ;
%@ 0  0  
%@ 3  2  
%@ Trial = [([]:[0/0, 0/0]->enroll->[]^[0/3, 0/0]),  ([]^[0/3, 0/0]->escalate->[0/3]:[0/0]),  ([0/3]:[0/0]->enroll->[0/3]^[0/3]),  ([0/3]^[0/3]->clamp->[0/3]:[0/3]),  ([0/3]:[0/3]->enroll->[... / ...]^[... / ...]),  ([... / ...]^[... / ...]->deescalate->[]:[...]),  ([]:[...]->enroll-> ... ^ ...),  (... ^ ... -> ... -> ...)],
%@ Translation = [1^0, 2^0, 2*2, 1:3, declare_mtd(0)],
%@ M =  ([0, 0], [3, 2]),
%@ Matrix = ['0  0  ', '3  2  '] ;
%@ 0  0  
%@ 3  3  
%@ Trial = [([]:[0/0, 0/0]->enroll->[]^[0/3, 0/0]),  ([]^[0/3, 0/0]->escalate->[0/3]:[0/0]),  ([0/3]:[0/0]->enroll->[0/3]^[0/3]),  ([0/3]^[0/3]->clamp->[0/3]:[0/3]),  ([0/3]:[0/3]->enroll->[... / ...]^[... / ...]),  ([... / ...]^[... / ...]->deescalate->[]:[...]),  ([]:[...]->enroll-> ... ^ ...),  (... ^ ... -> ... -> ...)],
%@ Translation = [1^0, 2^0, 2*3, 1:3, declare_mtd(0)],
%@ M =  ([0, 0], [3, 3]),
%@ Matrix = ['0  0  ', '3  3  '] ;
%@ 0  0  
%@ 0  2  
%@ Trial = [([]:[0/0, 0/0]->enroll->[]^[0/3, 0/0]),  ([]^[0/3, 0/0]->escalate->[0/3]:[0/0]),  ([0/3]:[0/0]->enroll->[0/3]^[0/3]),  ([0/3]^[0/3]->clamp->[0/3]:[0/3]),  ([0/3]:[0/3]->enroll->[... / ...]^[... / ...]),  ([... / ...]^[... / ...]->deescalate->[]:[...]),  ([]:[...]->enroll-> ... ^ ...),  (... ^ ... -> ... -> ...)],
%@ Translation = [1^0, 2^0, 2*2, 1:0, declare_mtd(1)],
%@ M =  ([0, 0], [0, 2]),
%@ Matrix = ['0  0  ', '0  2  '] .

%?- phrase(actions([] : [0/0]), Trial), phrase(condensed, Trial, Translation).
%@ Trial = [([]:[0/0]->enroll->[]^[0/3]),  ([]^[0/3]->clamp->[]:[0/3]),  ([]:[0/3]->enroll->[]^[_12274/6]),  ([]^[_12274/6]->stop->declare_mtd(0))],
%@ Translation = [1^0, 1*_12274, declare_mtd(0)],
%@ _12274 in 2..3 ;
%@ Trial = [([]:[0/0]->enroll->[]^[0/3]),  ([]^[0/3]->clamp->[]:[0/3]),  ([]:[0/3]->enroll->[]^[_17516/6]),  ([]^[_17516/6]->stop->declare_mtd(1))],
%@ Translation = [1^0, 1*_17516, declare_mtd(1)],
%@ _17516 in 0..1 ;
%@ Trial = [([]:[0/0]->enroll->[]^[1/3]),  ([]^[1/3]->stay->[]:[1/3]),  ([]:[1/3]->enroll->[]^[_27436/6]),  ([]^[_27436/6]->stop->declare_mtd(0))],
%@ Translation = [1^1, 1-_27496, declare_mtd(0)],
%@ _27436 in 2..4,
%@ 1+_27496#=_27436,
%@ 1+_27578#=_27436,
%@ _27496 in 1..3,
%@ _27578 in 1..3 ;
%@ Trial = [([]:[0/0]->enroll->[]^[1/3]),  ([]^[1/3]->stay->[]:[1/3]),  ([]:[1/3]->enroll->[]^[1/6]),  ([]^[1/6]->stop->declare_mtd(1))],
%@ Translation = [1^1, 1-0, declare_mtd(1)] ;
%@ Trial = [([]:[0/0]->enroll->[]^[_5876/3]),  ([]^[_5876/3]->stop->declare_mtd(0))],
%@ Translation = [1^_5876, declare_mtd(0)],
%@ _5876 in 2..3 ;
%@ false.
%@ ERROR: '$clause_term_position'/3: Cannot represent due to `int'
%@ ERROR: '$clause_term_position'/3: Cannot represent due to `int'


%% NB: The trivial trial with no doses is 'impossible' in this formulation.
%?- phrase(actions([] : []), Trial), phrase(condensed, Trial, Translation).
%@ false.

ndoses_state0(1, []:[0/0]).
ndoses_state0(N, []:[0/0|Init_1]) :-
    N #> 1,
    N_1 #= N - 1,
    ndoses_state0(N_1, []:Init_1).

%?- ndoses_state0(2,I).
%@ I = []:[0/0, 0/0] ;
%@ false.

% n_trials(+Drange, XY)
ndoses_ntrials(Drange, XY) :-
    Dmax in Drange, indomain(Dmax),
    ndoses_state0(Dmax, State0),
    findall(Trial,
	    (	phrase(actions(State0), Trial) % delays labeling
		,phrase(condensed, Trial, _)  % forces labeling
	    ),
	    Paths),
    length(Paths, N),
    XY = (Dmax, N).

%?- ndoses_ntrials(1..8, XY). 
%@ XY =  (1, 10) ;
%@ XY =  (2, 46) ;
%@ XY =  (3, 154) ;
%@ XY =  (4, 442) ;
%@ XY =  (5, 1162) ;
%@ XY =  (6, 2890) ;
%@ XY =  (7, 6922) ;
%@ XY =  (8, 16138) ; %% These all match counts from esc//2
%@ false.

%% TODO: Use sets to ensure exact matches between esc//2
%%       and condensed//0 translation of actions//1.


/* --------------------------------------------------

 Today, my efforts have to be focused on achieving easily vettable
 displays of TRIAL PATHS.

1. Introduce a format/4 predicate that maps (S0,Act,S) --> Term.
2. 

 */

%% Let's examine conduct of a trial with JUST 1 (nonzero) dose level.
%% (Trials always start off in state []:[0/0, ..., 0/0].)
%?- phrase(actions([] : [0/0]), Trial).
%@ Trial = [([]:[0/0]->enroll->[]^[_9666/3]),  ([]^[_9666/3]->clamp->[]:[_9666/3]),  ([]:[_9666/3]->enroll->[]^[_9732/6]),  ([]^[_9732/6]->stop->declare_mtd(0))],
%@ _9666 in 0..1,
%@ _9666+_9804#=_9732,
%@ _9804 in 1..3,
%@ _9732 in 2..4 ;
%@ Trial = [([]:[0/0]->enroll->[]^[1/3]),  ([]^[1/3]->stay->[]:[1/3]),  ([]:[1/3]->enroll->[]^[_16640/6]),  ([]^[_16640/6]->stop->declare_mtd(0))],
%@ _16640 in 2..4,
%@ 1+_16712#=_16640,
%@ _16712 in 1..3 ;
%@ Trial = [([]:[0/0]->enroll->[]^[_18968/3]),  ([]^[_18968/3]->stop->declare_mtd(0))],
%@ _18968 in 2..3 ;
%@ false.

%?- X #< Y, Y #< X.
%@ Y#=<X+ -1,
%@ X#=<Y+ -1.

% This models constraints as a GRAPH. "R-consistency"
?- [X,Y,Z] ins 0..2, all_distinct([X,Y,Z]), Z in 0..1, X in 0..1.
%@ Y = 2,
%@ X in 0..1,
%@ all_distinct([X, 2, Z]),
%@ Z in 0..1.
%@ X in 0..2,
%@ all_distinct([X, Y, Z]),
%@ Y in 0..2,
%@ Z in 0..1.


?- [X,Y,Z] ins 0..1, all_different([X,Y,Z]).
%@ X in 0..1,
%@ all_different([X, Y, Z]),
%@ Y in 0..1,
%@ Z in 0..1.

%% FASCINATING! I had never thought about this case, but the decision of the program
%% in retrospect looks reasonable. I do rather suspect that the spirit of the '=< 1/6'
%% definition would require that any dose 'reported out' of the trial (e.g. as RP2D)
%% ought to have been tested in 6 patients.

/* TODO ...

1. Let the 3+3 'rules' be DEMONSTRATED by complete *singleton* solutions
   to a sequence of Prolog QUERIES. This is the 'proof of correctness'
   that I was supposing I had 'approximately' obtained with the matching
   solution-set sizes for trial//1 and esc//2.

2. The RATIONALITY of any 3+3 program will be tested by how readily it
   supports generalization to other variants and extensions.

a) Demonstrate that the alternate 3+3 design (without de-escalation) can
   readily be obtained by a small (ideally, *parametrized*) modification.

b) Similarly introduce the accelerated titration. Ideally here, we might
   demonstrate COMPOSITIONALITY such that the 3+3 with an A.T. phase can
   be expressed by sequencing within a DCG.

*/
