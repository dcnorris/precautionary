% [Adapted from previous 'rolling.pl', on the assumption
% that it represents my latest thinking on this subject,
% although it dealt specifically with extending 'aliquots'
% to the case of rolling enrollment.]

% Implementing BOIN within a larger class of designs.

/*

Having at last read [1] with an eye toward implementation, I see that BOIN
lends a formal status to a whole class of dose-escalation designs which may
prove extremely easy to implement and explore using Prolog.

Table 2 in the paper is especially helpful, suggesting that a general class
of transition rules may be relatively easy to specify in CLP terms, without
recourse to whatever Real-analytic computations might be used as heuristics
for *finding* them:

             cumulative patients treated at current dose (n_j):
              1    2    3    4    5    6    7    8    9    10    11    12

lambda_{1,j} 0/1  0/2  0/3  0/4  0/5  0/6  0/7  1/8  1/9  1/10  1/11  1/12

lambda_{2,j} 1/1  2/2  2/3  2/4  3/5  3/6  4/7  4/8  5/9  5/10  5/11  6/12

elimination   -    -   3/3  3/4  3/5  4/6  4/7  4/8  5/9  5/10  6/11  6/12


The dose-elimination concept also readily invites computations over *lists*,
which are subject to truncation in the same manner.

Furthermore, with the 'local BOIN' design having so few free parameters (indeed,
under the default choice delta = +/-40%, the target toxicity rate is all that we
have left!), we obtain the real possibility of caching precomputed T arrays for
a finite number of possible designs! Thus, we might cache BOIN T arrays on a grid
of (TTR,D) in {0.1, 0.15, 0.2, 0.25, 0.3} x {3, 4, 5, 6, 7, 8} which amounts to
only 30 distinct arrays. (An interesting question is whether it might suffice to
cache only the D=8 arrays, and then collapse these to the D<8 arrays by a quick
array computation. In this case, it might well be reasonable to cache on a finer
TTR grid.)

Also helpful for situating this type of design within a 'taxonomy' is [2], with
its notion of a 'cumulative cohort design'. This is a design where the escalation
decisions are made solely on the tally *at the current dose*, without regard to
the tallies at other doses. This is the basic pattern to be implemented here,
albeit with elimination of overly-toxic doses which departs from strict CCD form.


1. Liu S, Yuan Y. Bayesian optimal interval designs for phase I clinical trials.
   J R Stat Soc C. 2015;64(3):507-523. doi:10.1111/rssc.12089

2. Ivanova A, Flournoy N, Chung Y. Cumulative cohort design for dose-finding.
   Journal of Statistical Planning and Inference. 2007;137(7):2316-2327.
   doi:10.1016/j.jspi.2006.07.009

--------------------------------------------------------------------------------

*/
:- use_module(library(clpz)).
:- use_module(library(pio)).
:- use_module(library(lists)).
:- use_module(library(dcgs)).
:- use_module(library(lambda)).
:- use_module(library(time)).
%@    true.

% Prefix op * for 'generalizing away' goals (https://www.metalevel.at/prolog/debugging)
:- op(920, fy, *).
*_.  % *Goal always succeeds

%% -----------------------------------------------------------------------------

% My initial emphasis is on generating all possible paths (CPE) for the BOIN
% design set forth in the table above. Although the BOIN design of [1] lacks
% any terminating principle except elimination of all doses, we do need such
% a rule here. The most natural formtype of rule, in view of the Table above,
% might be a 'stop-for-consensus' type of rule as found in package 'dtpcrm'.
% This is specified as a maximum number of patients to enroll at any 1 dose.

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

%% A fair enumeration of tallies will be helpful for developing
%% and testing comparison relations between tallies.

%% In general, I may require a max-enrollment *parameter*,
%% however unsightly it may be tagging along like this ...
tally_maxn(DLTs/Enrolled, MaxN) :-
    ground(MaxN),
    Enrolled in 0..MaxN,
    indomain(Enrolled),
    DLTs in 0..Enrolled,
    indomain(DLTs).

%?- tally_maxn(T/N, 6).
%@    T = 0, N = 0
%@ ;  T = 0, N = 1
%@ ;  T = 1, N = 1
%% ... as expected ...
%@ ;  T = 4, N = 6
%@ ;  T = 5, N = 6
%@ ;  T = 6, N = 6.

%% Having a 'default' max cohort size will make dev & test a bit easier:
tally(DLTs/Enrolled) :- tally_maxn(DLTs/Enrolled, 12).

%% tally/1 terminates:
%?- tally(_), false.
%@ false.

%?- tally(T/N).
%@    T = 0, N = 0
%@ ;  T = 0, N = 1
%@ ;  T = 1, N = 1
%@ ;  T = 0, N = 2
%@ ;  T = 1, N = 2
%% ... as expected ...
%@ ;  T = 8, N = 12
%@ ;  T = 9, N = 12
%@ ;  T = 10, N = 12
%@ ;  T = 11, N = 12
%@ ;  T = 12, N = 12.

%% How do tallies COMPARE?

/*
The essence of tally comparisons is one of REACHABILITY, i.e. EXISTENCE OF PATHS.
Generally, in fact, it is NON-REACHABILITY that supports any definite statement.
This might be stated correctly in terms of the DCG cohort//0. But this would not
be efficient. (Still, I could use such a statement as a correctness check!)

Strict inequalities like safer_than/2 mean that whichever argument is *earlier*
(i.e., has the smaller denominator) cannot under any conceivable path cross over
the other argument. Thus, the earlier is safer than the later if even a path of
100% toxicity that brings its denominator to the later's leaves its numerator
strictly less. Conversely, the later is safer if even a 100% non-toxic path from
the earlier does not make it look better than the later.
*/

safer_than(T0/N0, T1/N1) :-
    tally(T0/N0),
    tally(T1/N1),
    (	N0 #>= N1, T0 #< T1
    ;	N0 #< N1,
	MaxTox = N1 - N0,
	T0 + MaxTox #< T1
    ).

%?- safer_than(Q, 2/3).
%@    Q = 0/2
%@ ;  Q = 0/3
%@ ;  Q = 1/3
%@ ;  Q = 0/4
%@ ;  Q = 1/4
%@ ;  Q = 0/5
%@ ;  Q = 1/5
%@ ;  Q = 0/6
%@ ;  Q = 1/6
%@ ;  Q = 0/7
%@ ;  Q = 1/7
%@ ;  Q = 0/8
%@ ;  Q = 1/8
%@ ;  Q = 0/9
%@ ;  Q = 1/9
%@ ;  Q = 0/10
%@ ;  Q = 1/10
%@ ;  Q = 0/11
%@ ;  Q = 1/11
%@ ;  Q = 0/12
%@ ;  Q = 1/12
%@ ;  false.

%% Note e.g. that we cannot say safer_than(1/2, 2/3)
%% because we might enroll 1 with 1/2 ~~> 2/3, which
%% is not *strictly* less than 2/3.

:- op(900, xfx, &<).
&<(Q1, Q2) :- safer_than(Q1, Q2).
%?- Q &< 1/3.
%@    Q = 0/3
%@ ;  Q = 0/4
%@ ;  Q = 0/5
%@ ;  Q = 0/6
%@ ;  Q = 0/7
%@ ;  Q = 0/8
%@ ;  Q = 0/9
%@ ;  Q = 0/10
%@ ;  Q = 0/11
%@ ;  Q = 0/12
%@ ;  false.

%% Note that, under the new REACHABILITY semantics of tally comparisons,
%% we no longer make a claim such as 1/6 &< 1/3! What remains true is
%% only that 1/6 &=< 1/3.
%% NOTE: Reasoning about reachability amounts to reasoning about as-yet
%%       uninstantiated toxicity assessments, and therefore accords
%%       perfectly with this essential strength of Prolog.
%% TODO: It seems likely that reachability semantics are more stringent
%%       than my previous formalism, so that we now can say fewer things.
%%       I ought to prove this---or let Prolog prove it for me.
%% TODO: What is the relationship (if any) between meta-interpretation
%%       and this change of tally-comparison semantics? Does MI affect
%%       only the IMPLEMENTATION---and not the meaning---of a program?
%% TODO: Given how far this reachability semantics departs from normal
%%       pharmacologic intuitions, perhaps even words like 'safer_than'
%%       should be reconsidered. I may well have abandoned all hope of
%%       capturing such intuitions in these predicates!

noworse_than(T0/N0, T1/N1) :-
    tally(T0/N0),
    tally(T1/N1),
    (	N0 #>= N1, T0 #=< T1
    ;	N0 #< N1,
	MaxTox = N1 - N0,
	T0 + MaxTox #=< T1
    ).

%?- safer_than(1/3, 2/3).
%@    true
%@ ;  false.

%?- noworse_than(1/3, 2/3).
%@    true
%@ ;  false.

%?- tally(1/3), tally(2/3), 3 #>= 3, 1 #=< 2.
%@    true.


%?- noworse_than(Q, 2/3).
%@    Q = 0/1
%@ ;  Q = 0/2
%@ ;  Q = 1/2
%@ ;  Q = 0/3
%@ ;  Q = 1/3
%@ ;  Q = 2/3
%@ ;  Q = 0/4
%@ ;  Q = 1/4
%@ ;  Q = 2/4
%@ ;  Q = 0/5
%@ ;  Q = 1/5
%@ ;  Q = 2/5
%@ ;  Q = 0/6
%@ ;  Q = 1/6
%@ ;  Q = 2/6
%@ ;  Q = 0/7
%@ ;  Q = 1/7
%@ ;  Q = 2/7
%@ ;  Q = 0/8
%@ ;  Q = 1/8
%@ ;  Q = 2/8
%@ ;  Q = 0/9
%@ ;  Q = 1/9
%@ ;  Q = 2/9
%@ ;  Q = 0/10
%@ ;  Q = 1/10
%@ ;  Q = 2/10
%@ ;  Q = 0/11
%@ ;  Q = 1/11
%@ ;  Q = 2/11
%@ ;  Q = 0/12
%@ ;  Q = 1/12
%@ ;  Q = 2/12
%@ ;  false.

%?- safer_than(T/N, 1/2).
%@    T = 0, N = 2
%@ ;  T = 0, N = 3
%@ ;  T = 0, N = 4
%@ ;  T = 0, N = 5
%@ ;  T = 0, N = 6
%@ ;  T = 0, N = 7
%@ ;  T = 0, N = 8
%@ ;  T = 0, N = 9
%@ ;  T = 0, N = 10
%@ ;  T = 0, N = 11
%@ ;  T = 0, N = 12
%@ ;  false.

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

%?- Q &= R.
%@    Q = 0/1, R = 0/1
%@ ;  Q = 0/1, R = 0/2
%@ ;  Q = 0/1, R = 0/3
%@ ;  Q = 0/1, R = 0/4
%@ ;  Q = 0/1, R = 0/5
%@ ;  Q = 0/1, R = 0/6
%@ ;  Q = 0/1, R = 0/7
%@ ;  Q = 0/1, R = 0/8
%@ ;  Q = 0/1, R = 0/9
%@ ;  Q = 0/1, R = 0/10
%@ ;  Q = 0/1, R = 0/11
%@ ;  Q = 0/1, R = 0/12
%@ ;  Q = 1/1, R = 1/1
%@ ;  Q = 1/1, R = 2/2
%@ ;  Q = 1/1, R = 3/3
%@ ;  Q = 1/1, R = 4/4
%@ ;  Q = 1/1, R = 5/5
%@ ;  Q = 1/1, R = 6/6
%@ ;  ...

%?- Q &< R.
%@    Q = 0/1, R = 1/1
%@ ;  Q = 0/1, R = 2/2
%@ ;  Q = 0/1, R = 3/3
%@ ;  Q = 0/1, R = 4/4
%@ ;  Q = 0/1, R = 5/5
%@ ;  Q = 0/1, R = 6/6
%@ ;  ...

% Strict tally-safer (&<) excludes tally-equals (&=)
%?- time((Q &< R, Q &= R)).
%@    % CPU time: 56.976 seconds
%@ false.

% (&=<) and '&>' are MUTUALLY EXCLUSIVE:
%?- time((Q &< R, R &=< Q)).
%@    % CPU time: 63.436 seconds
%@ false.


%% Here was my motivation for writing down on_par/2 in the first place:
%?- Q &= 1/3. %% 1/3 is (roughly) the 3+3 design's 'target toxicity rate'
%@    Q = 1/3
%@ ;  Q = 2/6
%@ ;  Q = 3/9
%@ ;  Q = 4/12
%@ ;  false.

%?- cohort_tally([0,0,T], Q), safer_than(Q, 2/3).
%@    Q = 0/3, T = 0
%@ ;  Q = 1/3, T = 1
%@ ;  false.

%% -----------------------------------------------------------------------------

% The BOIN rules suggest that comparisons of tallies might be generalized to
% comparisons between a tally and a *list* of tallies. Thus, e.g. the rule:
%
% lambda_{1,j} 0/1  0/2  0/3  0/4  0/5  0/6  0/7  1/8  1/9  1/10  1/11  1/12
%
% generates interest in comparisons such as T/N &=< [0/1, 0/2, 0/3, 0/4, ...].
%
% Does this comparison have an unambiguous or obvious meaning that might allow
% us to condense the list? It would be a shame to have to completely list all
% possible fractions in a sequence like this!
%
% With &=< at least, perhaps we would like to write instead:
%   T/N &=< [0/1, 1/8, 2/15],
% or possibly even abbreviate this as:
%   T/N &=< [1, 8, 15],
% although I am inclined to reject the latter as more 'cute' than *clear*.
%
% The comparison would then hold if ANY of the mapped comparisons holds.
%
% What about the flip-side?
%
% lambda_{2,j} 1/1  2/2  2/3  2/4  3/5  3/6  4/7  4/8  5/9  5/10  5/11  6/12
%
% Here, we require comparisons like T/N &>= [1/1, 2/2, 2/3, ...] or possibly
% T/N &>= [1/1, 2/4, 3/6, ...] or even T/N &>= [1, 4, 6, ...].
%
% Again, this complex comparison holds if ANY component comparison holds.

% Define what a BOIN-type threshold is ...

zip([], [], []).
zip([X|Xs], [Y|Ys], [X/Y | XsYs]) :-
    zip(Xs, Ys, XsYs).

%?- zip([1,2,3], [4,5,6], Z).
%@    Z = [1/4,2/5,3/6].

%?- zip(T, N, [1/4,2/5,3/6]).
%@    T = [1,2,3], N = [4,5,6]
%@ ;  false.

tox_boundary(Bdy) :-
    zip(Ts, Ns, Bdy),
    sort(Ns, Ns), % TODO: Do I gain anything by thus removing degeneracy?
    maplist(\X^(X #> 0), Ns), % TODO: Is there a clpz idiom for this?
    %% Note there are further constraints on what might be considered
    %% reasonable tox boundaries. There shouldn't be duplicated Ns, e.g.,
    %% and all the Ts should be non-negative. But perhaps these are best
    %% postponed until some kind of *search* over the space of possible
    %% tox boundaries becomes desirable.
    maplist(#=<, Ts, Ns).

%?- tox_boundary([0/2, 1/4]).
%@    true
%@ ;  false.

%?- tox_boundary([0/5, 1/4]).
%@ false.

%% TODO: Why does sortedness of Ns not appear in constraints below?
%?- tox_boundary(Bdy).
%@    Bdy = []
%@ ;  Bdy = [_B/_A], clpz:(_A#>=_B), clpz:(_A in 1..sup)
%@ ;  Bdy = [_B/_A,_D/_C], clpz:(_A#>=_B), clpz:(_C#>=_D), clpz:(_A in 1..sup), clpz:(_C in 1..sup)
%@ ;  Bdy = [_B/_A,_D/_C,_F/_E], clpz:(_A#>=_B), clpz:(_C#>=_D), clpz:(_E#>=_F), clpz:(_A in 1..sup), clpz:(_C in 1..sup), clpz:(_E in 1..sup)
%@ ;  ...

/*
The question about a cumulative-cohort tally vis-à-vis a given threshold,
is whether it has touched or crossed the threshold toward 'extreme' values.
The past participle 'hit' seems reasonable, since it covers both the case
of touching and of having 'broken through'.
*/

%% TODO: Consider not leaning on (&<)-type tally comparisons, but rather
%%       going directly to CLP(Z) constraints. This will be faster, and
%%       also enable me to avoid superfluous choice points, etc.
%% OOH! Maybe I'm discovering that all the above efforts were unnecessary!
%%      Could it be that this file can just about start HERE?

/*
If I do 'cut out the middleman' in this way, what is the right representation
for the toxicity boundaries?

lambda_{1,j} 0/1  0/2  0/3  0/4  0/5  0/6  0/7  1/8  1/9  1/10  1/11  1/12
lambda_{2,j} 1/1  2/2  2/3  2/4  3/5  3/6  4/7  4/8  5/9  5/10  5/11  6/12
elimination   -    -   3/3  3/4  3/5  4/6  4/7  4/8  5/9  5/10  6/11  6/12

RELATIVE to a maximum cohort size of 12, the above boundaries gain sufficient
and minimal representation as the following lists:

lambda_1 [0/1, 1/8]
lambda_2 [1/1, 2/4, 3/6, 4/8, 5/10, 6/12]
elim     [3/5, 4/8, 5/10, 6/12]

This is because lambda_1 is a lower bounds, and so the removed limits were
non-binding. Likewise, since lambda_2 and elim were upper bounds the removed
limits were non-binding.

In general, minimal representations of tox boundaries will consist of lists
of tallies BdyRep satisfying:

minimal_maxn(BdyRep, MaxN) :-
    zip(Ts, Ns, BdyRep),
    T1 #=< Tn,
    Ts in T1..Tn,
    sort(Ns, Ns),
    maplist(\X^(X #=< MaxN), Ns).


*/

%% I unapologetically assume that the Bdy list is a minimal representation.
hit_floor(_/_, []) :- false.
hit_floor(T/N, [Tb/Nb | Bdys]) :-
    %minimal_maxn([Tb/Nb|Bdys], 12), %% TODO: Try activating this assertion
    %% TODO: The correctness of this predicate may depend only on the
    %%       sortedness of Ts and Ns. Work out the details!
    (	T #= Tb, N #>= Nb
    ;	T #> Tb, hit_floor(T/N, Bdys)
    ).
    
%?- hit_floor(0/3, [0/1, 0/2, 0/3, 0/4, 0/5, 0/6, 0/7, 1/8, 1/9, 1/10, 1/11, 1/12]).
%@    true
%@ ;  false.

%?- hit_floor(0/3, [0/1, 1/8]).
%@    true
%@ ;  false.

hit_ceiling(_/_, []) :- false.
hit_ceiling(T/N, [Tb/Nb | Bdys]) :-
    %minimal_maxn([Tb/Nb|Bdys], 12), %% TODO: Try activating this assertion
    (	T #= Tb, N #=< Nb
    ;	T #> Tb, hit_ceiling(T/N, Bdys)
    ).

%?- hit_ceiling(3/5, [1/1, 2/4, 3/6, 4/8, 5/10, 6/12]).
%@    true
%@ ;  false.

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Pick up from here...

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

/*
Above text is included for comparison with current outlook.
*/

%% Instead of carrying the tox boundaries along as parameters,
%% let's simply assert them into the database ...
%% Note that this naturally invites consideration of a DSL
%% in which such design rules could be expressed generally!

escalate(Q) :- hit_floor(Q, [0/1, 1/8]).
deescalate(Q) :- hit_ceiling(Q, [1/1, 2/4, 3/6, 4/8, 5/10, 6/12]).
remove(Q) :- hit_ceiling(Q, [3/5, 4/8, 5/10, 6/12]).

%% TODO: Note the duplication between remove/2 and deescalate/2.
%%       Consider whether these should be refined to avoid this.
%% TODO: Might a zcompare/3-like idiom further the cause of purity?

/*
Even as I set out to write the first lines of this predicate,
I can see the value of having complementary &** predicates
that enable me to state definitively whether a boundary has
been breached, or not.
In a way, what I need is a kind of zcompare/3 for tox boundaries!
Let me not impose burdens on state-machine code, that more properly
belong to basic comparison operators!
I don't think I have yet fully worked out the MEANING of toxicity
boundaries as used in BOIN.
Keep in mind that I want the code ultimately to be verifiable
on inspection by Ying Yuan and others -- at least where it is not
self-verifying by proofs executed in Prolog itself.
*/

%% NB: The left-hand list of Ls ^ Rs is sorted in descending order.
%%     So heads L and R in [L|Ls] ^ [R|Rs] are *adjacent* doses.

cohort_full(N, yes) :- N #>= 12.
cohort_full(N, no) :- N in 0..11.
cohort_full(_/N, YesNo) :- cohort_full(N, YesNo).

enroll(T0/N0, T1/N1) :-
    cohort_full(N0, no),
    N1 #= N0 + 1,
    T in 0..1, % T is the pending tox assessment of newly-enrolled patient
    T1 #= T0 + T.

%% Overload enroll/2 on lists of tallies, enrolling head tally
enroll([], _) :- false. % (**)
enroll([T0/N0 | Qs], [T1/N1 | Qs]) :-
    enroll(T0/N0, T1/N1).

length_plus_1(Ls, MTD) :-
    length(Ls, MTD_1),
    MTD #= MTD_1 + 1.

stay(Ls ^ [R | Rs], Ls ^ [R1 | Rs]) :- enroll(R, R1).
stay(Ls ^ [_/N | _], declare_mtd(MTD)) :- cohort_full(N, yes),
					  length_plus_1(Ls, MTD).

escalate(Ls ^ [T/N], State) :- % NB: this is a 'clamped' situation
    stay(Ls ^ [T/N], State).

escalate(Ls ^ [Q, R | _], [Q | Ls] ^ [R1 | _]) :- enroll(R, R1).
escalate(Ls ^ [Q, R | _], declare_mtd(MTD)) :- cohort_full(R, yes),
						length_plus_1([Q|Ls], MTD).

deescalate([] ^ _, declare_mtd(0)). % deescalate from already-lowest dose

%% TODO: Use 'guard' to parse this into separate clauses.
deescalate([L | Ls] ^ Rs, Ls ^ [L1 | Rs]) :- enroll(L, L1).
deescalate([L | Ls] ^ _, declare_mtd(MTD)) :- cohort_full(L, yes),
					      length(Ls, MTD).

%% TODO: Note how TRIVIAL state0_action_state/3 has become.
%%       Does this demand refactoring, or does it represent other
%%       types of opportunity -- even meta-interpretation?

state0_action_state(Ls ^ [R | Rs], escalate, State) :-
    escalate(R),
    escalate(Ls ^ [R | Rs], State).

state0_action_state(Ls ^ [R | Rs], deescalate, State) :-
    deescalate(R),
    deescalate(Ls ^ [R | Rs], State).

state0_action_state(Ls ^ [R | Rs], stay, Ls ^ [R1 | Rs]) :-
    (	escalate(R) -> false
    ;	deescalate(R) -> false
    ;	%% Until I gain greater (e.g., zcompare/3-style) control over the
	%% boundary-hitting concept, I must treat 'stay' as a DEFAULT action:
	enroll(R, R1)
    ).

%?- state0_action_state([] ^ [0/0, 0/0, 0/0], Action, State).
%@ caught: error(existence_error(procedure,between/3),between/3)
