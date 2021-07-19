%% Exploration of non-monotonicity in CLP(ℤ)

:- use_module(library(clpz)).

%% I have extracted here the most fundamental relation between tallies,
%% implemented in 'tally.pl' -- namely, REACHABILITY. For example,
%% could a tally of 3/6 (3 toxicities out of 6 patients dosed) ever be
%% ensue during a trial where we've already dosed 4 patients, with zero
%% toxicity? That is the question whether 0/4 ~~ 3/6, which is FALSE.
%% Conversely, 3/6 is a reachable tally from a state of the trial where
%% a 1/4 tally has been observed.
:- op(900, xfx, ~~).
~~(T1/N1, T2/N2) :- % REACHABILITY
    %% The following conditions translate as follows, under 3 scenarios:
    %% (=) If N1==N2, then we get T2 =< T1 =< T2 which holds iff T1==T2.
    %% (<) If N1 < N2, these conds translate to T1 =< T2 =< T1 + MaxTox.
    %% (>) If N1 > N2, these conds translate to T2 =< T1 =< T2 + MaxTox.
    T2 #=< T1 + max(0, N2 - N1),
    T1 #=< T2 + max(0, N1 - N2).


%% As noted above ...

%?- 0/4 ~~ 3/6.
%@ false. % even if patients 5 & 6 have tox, worst tally would be 2/6

%?- 1/4 ~~ 3/6.
%@    true. % if patients 5 & 6 both experience toxicities, we get 1+2=3

%% To get to the heart of the matter, I need to ask for all tallies
%% from which a given tally could be reached, then show 'we missed one'.

%?- T/N ~~ 2/6, N in 0..5, indomain(N), T in 0..N.
%@    T = 0, N = 0
%@ ;  N = 1, clpz:(T+5#=_A), clpz:(_A in 5..6), clpz:(T in 0..1)
%@ ;  N = 2, clpz:(T+4#=_A), clpz:(_A in 4..6), clpz:(T in 0..2)
%@ ;  N = 3, clpz:(T+3#=_A), clpz:(_A in 3..5), clpz:(T in 0..2)
%@ ;  N = 4, clpz:(T+2#=_A), clpz:(_A in 2..4), clpz:(T in 0..2)
%@ ;  N = 5, clpz:(T+1#=_A), clpz:(_A in 2..3), clpz:(T in 1..2).

%?- T/N ~~ 2/6, N in 0..5, indomain(N), T in 0..N, indomain(T).
%@    T = 0, N = 0
%@ ;  T = 0, N = 1
%@ ;  T = 1, N = 1
%@ ;  T = 0, N = 2
%@ ;  T = 1, N = 2
%@ ;  T = 2, N = 2
%@ ;  T = 0, N = 3
%@ ;  T = 1, N = 3
%@ ;  T = 2, N = 3
%@ ;  T = 0, N = 4
%@ ;  T = 1, N = 4
%@ ;  T = 2, N = 4
%@ ;  T = 1, N = 5
%@ ;  T = 2, N = 5.

%% In investigating reachability, could we ever get tripped up by the
%% non-monotonicity of clpz?

% Could 1/6 ever be reached from an 'earlier' tally T/N with N<6 but T>1?

%?- T/N ~~ 1/6, T #> 1, N #>= T, N #< 6.
%@ false.

% What if I take steps toward indicating my view that N is composed of
% the sum of X's and O's (toxicities and non-toxicities, as commonly
% depicted in dose-escalation diagrams).

%   Added goal
%   vvvvvvvvv
%?- N = X + O, T/N ~~ 1/6, T #> 1, N #>= T, N #< 6.
%@    N = X+O, clpz:(_A#=max(0,_B)), clpz:(_C#=max(0,_D)), clpz:(_B+_E#=6), clpz:(_D+6#=_F), clpz:(T+_A#=_G), clpz:(X+O#=_E), clpz:(X+O#=_F), clpz:(X+O#=_H), clpz:(1+_C#=_I), clpz:(T#=<X+O), clpz:(_A#>=_B), clpz:(_C#>=_D), clpz:(_I#>=T), clpz:(_A in 0..sup), clpz:(_C in 0..sup), clpz:(_I in 2..sup), clpz:(_G in 2..sup), clpz:(T in 2..sup), clpz:(_H in inf..5).

% INTERESTING! I told clpz 'more stuff' (even if I didn't really work
% it out fully, by explaining that O and X are non-negative integers),
% and that left clpz unable to 'tidy up loose ends' during constraint
% propagation?

% I appreciate now the distinction between an ANSWER and SOLUTIONS.
% What clpz gives now is a more 'diffuse' answer, but whether that
% answer contains any solutions that could be produced upon labeling
% is a different matter.

% It's hard for me to believe that clpz could possibly yield such a
% false 'solution'. It seems that most of the nonmonotonicity here
% arises out of the unhelpful/distracting REPRESENTATION (N = X + O)
% which I have suggested to clpz. Surely, no actual LABELING of the
% variables here could sneak past the constraints!

% Let me at least try to demonstrate (or refute!) this by labeling.

% The first step in labeling is restricting the variables to finite
% domains:

%?- N = X + O, T/N ~~ 1/6, T #> 1, N #>= T, N #< 6, N in 1..24.
%@ caught: error(type_error(integer,_1384253+_1384255),unknown(_1384253+_1384255)-1)

%?- N = X + O, T/N ~~ 1/6, T #> 1, N #>= T, N #< 6, N #> 0.
%@    N = X+O, clpz:(_A#=max(0,_B)), clpz:(_C#=max(0,_D)), clpz:(_B+_E#=6), clpz:(_D+6#=_F), clpz:(T+_A#=_G), clpz:(X+O#=_E), clpz:(X+O#=_F), clpz:(X+O#=_H), clpz:(X+O#=_I), clpz:(1+_C#=_J), clpz:(T#=<X+O), clpz:(_A#>=_B), clpz:(_C#>=_D), clpz:(_J#>=T), clpz:(_A in 0..sup), clpz:(_C in 0..sup), clpz:(_J in 2..sup), clpz:(_G in 2..sup), clpz:(T in 2..sup), clpz:(_H in inf..5), clpz:(_I in 1..sup).

% ANOTHER SURPRISE! Adding superfluous constraints N in 1..24 or N #> 0
% had distinct effects. Is (N in 1..24) NOT same as (N #> 0, N #=< 24)? 

%?- N = X + O, T/N ~~ 1/6, T #> 1, N #>= T, N #< 6, #N #> 0.
%@ caught: error(type_error(integer,_1384251+_1384253),must_be/2)

% How do I reflect on this, AS A USER of Prolog? It looks to me
% as if (N = X + O) represented a 'typo', which clpz was able to
% detect under the more stringent syntactical demands it imposes.

%     ,---- oops! forgot to type the '#'
%     v
%?- N = X + O, T/N ~~ 1/6, T #> 1, N #>= T, N #< 6.
%@    N = X+O, clpz:(_A#=max(0,_B)), clpz:(_C#=max(0,_D)), clpz:(_B+_E#=6), clpz:(_D+6#=_F), clpz:(T+_A#=_G), clpz:(X+O#=_E), clpz:(X+O#=_F), clpz:(X+O#=_H), clpz:(1+_C#=_I), clpz:(T#=<X+O), clpz:(_A#>=_B), clpz:(_C#>=_D), clpz:(_I#>=T), clpz:(_A in 0..sup), clpz:(_C in 0..sup), clpz:(_I in 2..sup), clpz:(_G in 2..sup), clpz:(T in 2..sup), clpz:(_H in inf..5).

%     ,---- phew! I remembered it this time
%     v
%?- N #= X + O, T/N ~~ 1/6, T #> 1, N #>= T, N #< 6.
%@ false.

%?- N #= X + O.
%@    clpz:(X+O#=N).


% My 'conjecture' is that the PRACTICAL effect of clpz:monotonic
% is to impose unpleasant syntactical demands on the programmer,
% merely in order to detect typos.

% One way to operationalize this conjecture would be to demonstrate
% that syntactical analysis of such programs would suffice to flag
% all such typos. Why can't clpz recognize that N appears on LHS of
% a #=/2, and 'put a # on it' automatically? What *legitimate* goal
% would thereby be made impossible to write?

%?- (X+1) in 0..10.
%@ caught: error(type_error(integer,_1384249+1),unknown(_1384249+1)-1)

%@    true.

%?- must_be(integer, X).
%@ caught: error(instantiation_error,must_be/2)

%?- N = X + O, T/N ~~ 1/6, T #> 1, N #>= T, N #< 6, X in 0..12, O in 0..24.
%@    N = X+O, clpz:(_A#=max(0,_B)), clpz:(_C#=max(0,_D)), clpz:(_D+6#=_E), clpz:(_B+_F#=6), clpz:(X+O#=_G), clpz:(X+O#=_E), clpz:(X+O#=_F), clpz:(T+_A#=_H), clpz:(1+_C#=_I), clpz:(T#=<X+O), clpz:(_A#>=_B), clpz:(_C#>=_D), clpz:(_I#>=T), clpz:(_G in 0..5), clpz:(_E in 0..10), clpz:(_F in 0..10), clpz:(_D in-6..4), clpz:(_B in-4..6), clpz:(_A in 0..6), clpz:(_C in 1..4), clpz:(_H in 2..11), clpz:(X in 0..5), clpz:(_I in 2..5), clpz:(T in 2..5), clpz:(O in 0..5).
%@    N = X+O, clpz:(_A#=max(0,_B)), clpz:(_C#=max(0,_D)), clpz:(_D+6#=_E), clpz:(_B+_F#=6), clpz:(T+_A#=_G), clpz:(X+O#=_H), clpz:(X+O#=_E), clpz:(X+O#=_F), clpz:(1+_C#=_I), clpz:(T#=<X+O), clpz:(_A#>=_B), clpz:(_C#>=_D), clpz:(_I#>=T), clpz:(_H in-10..5), clpz:(_E in-10..17), clpz:(_F in-10..17), clpz:(_D in-16..11), clpz:(_B in-11..16), clpz:(_A in 0..16), clpz:(_C in 1..11), clpz:(_G in 2..28), clpz:(O in-10..5), clpz:(_I in 2..12), clpz:(T in 2..12), clpz:(X in 0..12).
%@    N = X+O, clpz:(_A#=max(0,_B)), clpz:(_C#=max(0,_D)), clpz:(_B+_E#=6), clpz:(_D+6#=_F), clpz:(T+_A#=_G), clpz:(X+O#=_E), clpz:(X+O#=_F), clpz:(X+O#=_H), clpz:(1+_C#=_I), clpz:(T#=<X+O), clpz:(_A#>=_B), clpz:(_C#>=_D), clpz:(_I#>=T), clpz:(_A in 0..sup), clpz:(_C in 0..sup), clpz:(_I in 2..sup), clpz:(_G in 2..sup), clpz:(T in 2..sup), clpz:(_H in inf..5), clpz:(O in inf..4), clpz:(X in 1..sup).
%@    N = X+O, clpz:(_A#=max(0,_B)), clpz:(_C#=max(0,_D)), clpz:(_B+_E#=6), clpz:(_D+6#=_F), clpz:(T+_A#=_G), clpz:(X+O#=_E), clpz:(X+O#=_F), clpz:(X+O#=_H), clpz:(1+_C#=_I), clpz:(T#=<X+O), clpz:(_A#>=_B), clpz:(_C#>=_D), clpz:(_I#>=T), clpz:(_A in 0..sup), clpz:(_C in 0..sup), clpz:(_I in 2..sup), clpz:(_G in 2..sup), clpz:(T in 2..sup), clpz:(_H in inf..5).



%?- assertz(clpz:monotonic).
%@    true.

% What if I try that same query again, AFTER making the above assertion?

%?- N = X + O, T/N ~~ 1/6, T #> 1, N #>= T, N #< 6.
%@ caught: error(instantiation_error,instantiation_error(unknown(_1384246),1))

%?- #N #= #X + #O, #T/ #N ~~ 1/6, #T #> 1, #N #>= #T, #N #< 6.
%@ false.

% But what about those first queries?
%?- N = X + O.
%@    N = X+O.

%?- #N #= #X + #O.
%@    clpz:(#X+ #O#= #N).

/*

"You meant to exclude these solutions, but you didn't." -- MT

"By using clpz:monotonic, you enforce that the clean representation is used."

"If you have doubts about the completeness of the answer, then [you're done for.]"

"If something is lost within the system... then I have to reimplement the whole system."

"If I get solutions that I meant to exclude, then I can exclude them. I can write a verifier..."

"But if the system doesn't REPORT solutions that it does accept..." then how could we trust it?

"The key selling point is that we get everything that there is."

%% TODO: Stress generalization/specialization techniques. Must demonstrate
%%       we have a program that we can generalize readily.

%% TODO: The take-away from the call and the 'factorial.pl' exploration is (I think)
%%       that we don't want the system to invite 'magical thinking', allowing us to
%%       freely discuss solutions with representations that WERE NEVER INTENDED.

*/
