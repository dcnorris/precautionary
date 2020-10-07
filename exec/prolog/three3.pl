% Attempt an enumeration of ALL 3+3 trials
:- use_module(library(clpfd)).

% Prefix op * for 'generalizing away' goals (https://www.metalevel.at/prolog/debugging)
:- op(920, fy, *). *_.  % *Goal always succeeds

/*

Aim initially to describe the 'Design 1' of Simon &al (1997).
The main difficulty in describing this design seems to be the
arbitrarily distant lookback to earlier cohorts, as required
to ensure that at most 1/6 toxicities occur at the declared
MTD.

What may suffice, however, is to define simple (Markov?) rules,
then PROVE that desired properties (such as '< 1/6') result.

ALTERNATIVELY (and surely *better*!), I might try to state the
desiderata, and allow the dose-escalation process to EMERGE as a
*result* of those specifications. For example, Skolnik &al (2008)
offer this definition:

> The MTD is defined as the dose level at which none or one
> of six participants (0% to 17%) experience a DLT, when at least
> two of three to six participants (33% to 67%) experience a DLT
> at the next highest dose.

To support such a declarative specification, I would have to
implement a scan of the full list generated thus far, perhaps
using semicontext notation as in the tree leaf-counting example
by Triska. The goals at any point would then be:
1. Can I declare an MTD?
2. If not, what should the next cohort be?
3. If I can't create another cohort, then declare "no MTD found".

NOTE HOW WONDERFULLY GOAL-DIRECTED THIS WOULD BE!!

Hmmm!! I think I begin to see implicit concepts REVEALED as I try
to think about the 3+3 in this manner:
a) A dose considered 'too toxic', initially this is (Hi+1)
b) A dose known to be 'safe', initially this is 0

Might I even get away with calling the 'safe' dose MTD?

REFERENCES

Skolnik, Jeffrey M., Jeffrey S. Barrett, Bhuvana Jayaraman, Dimple Patel,
and Peter C. Adamson. “Shortening the Timeline of Pediatric Phase I Trials:
The Rolling Six Design.” JCO 26, no. 2 (Jan 10, 2008): 190–95.
https://doi.org/10.1200/JCO.2007.12.7712.

=====

In all honesty, though, perhaps I should acknowledge that a
'running range' is carried forward implicitly as a state variable,
providing context for (de-)escalation decisions.

esc(D, Lo..Hi) --> ...

*/

% Let esc(D, Lo..Hi) *describe* a list of 3+3 cohorts starting from dose D in Lo..Hi.
% Read esc(D, Lo..Hi) as escalating FROM D to min(D+1,Hi).

%% TODO: WLOG set Lo == 1 by convention.

esc(Hi, Lo..Hi) --> [Hi * T], { Lo #< Hi, T in 0..3, indomain(T) }, % * is a ^ that went *splat*
		    (   { T in 0..1 } -> [mtd_notfound(Hi)]
		    ;   des(Hi, Lo)
		    ).
esc(D, Lo..Hi) --> [D1 ^ T], { D1 #= D + 1, D1 in Lo..Hi, T in 0..3, indomain(T) },
		   (   { T #= 0 } -> esc(D1, Lo..Hi)
		   ;   { T #= 1 } -> sta(D1, Lo..Hi)
		   ;   des(D1, Lo)
		   ).

% TODO: Do I really need the special case of sta(D, _..D)?
sta(D, _..D) --> [D - 0], [mtd_notfound(D)].
sta(D, Lo..Hi) --> [D - 0], { D #< Hi, D in Lo..Hi, indomain(D) }, esc(D, D..Hi).
sta(D, Lo.._) --> [D - T], { T in 1..3, indomain(T) }, des(D, Lo).

% As a mirror image of esc, des(D, Lo) moves downward FROM D to max(D-1,Lo).
% NB: Implicit in de-escalation to D-1 is that Hi #= D - 1.
des(D, Lo) --> { D_1 #= D - 1 },
	       (   { D_1 #= Lo } -> [declare_mtd(Lo)]
	       ;   [D_1 : T], { T in 0..3, indomain(T) },
		   (   { T in 0..1 } -> [declare_mtd(D_1)]
		   ;   des(D_1, Lo)
		   )
	       ).

%% TODO: Write cohorts as (Dose ^ Tox / Enrolled).
%%       This will enable cohorts of size other than 3,
%%       as needed e.g. to describe an accelerated titration phase.

% n_trials(+Drange, XY)
n_trials(Drange, XY) :-
    Dmax in Drange, indomain(Dmax),
    findall(Tr, phrase(esc(0, 0..Dmax), Tr), Trials),
    length(Trials, N),
    XY = (Dmax, N).
