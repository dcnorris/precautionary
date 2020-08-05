% Attempt an enumeration of ALL 3+3 trials


% Predefined dose levels
dose(D) :- between(1, 5, D).

% 0 - 3 DLTs per cohort
tox(N) :- between(0, 3, N).

% TODO: Distinguish an MTD-not-found terminal case.

% Trial terminates with declared MTD:
declare_mtd(D, [D ^ 1 | Prev]) :-
    stay(D, Prev).
declare_mtd(D_1, [D ^ N | Prev]) :-
    tox(N), N > 1,
    dose(D),
    D_1 is D - 1,
    escalate(D, Prev).

stay(D, [D ^ 1 | Prev]) :-
    dose(D),
    escalate(D, Prev).

escalate(1, []). % base case for recursion
escalate(D1, [D ^ 0 | Prev]) :-
    dose(D),
    D1 is D + 1,
    dose(D1),
    (   escalate(D, Prev)
    ;   stay(D, Prev)
    ).


% ---- Attempt a more concise program ----

declare(D, [D ^ 1, D ^ 1 | E]) :-
    escalation(D, E).
declare(D_1, [D ^ T | E]) :-
    escalation(D, E),
    tox(T), T > 1,
    D_1 is D - 1.

escalation(1, []).

escalation(D1, [D ^ 0 | E]) :-
    dose(D1),
    D is D1 - 1,
    D > 0,
    escalation(D, E).

% In these goals, should I be unifying E with
% a more detailed representation?

% Do I need 2 different predicates, next/2 and trial/2?
next(D, E) :-
    (	up(D, E)
    ;	sta(D, E)
    ).
up(1, []).
up(U, [D^0 | E]) :- dose(D), U is D+1, dose(U), next(D,E).
sta(D, [D^1 | E]) :- up(D, E).
mtd(D, [D^1 | E]) :- sta(D, E).
mtd(D_, [D^T | E]) :- tox(T), T>1, up(D, E), D_ is D-1.
result(R, E) :- % TODO: Make #doses a parameter
    (	mtd(D, E), R = mtd(D)
    ;	up(6, E), R = mtd_not_found
    ).

% TODO: Try a formulation via PRODUCTION RULES!

/* I think the way this works is, I have terms on the left of '-->',
 * each describing an event that may happen during dose escalation.
*/

% Aha! Interesting, that I can appreciate the grammar of 3+3 trials
% as being non-context-free!

istox(T) :- member([3,2,1,0], T). % count down from high to low tox
trial --> [1^T, mtd(0)], { istox(T), T>0 }.
trial --> [D^0], { dose(D) }, esc(D).
esc(D) --> [D1^T], { dose(D1), D1 is D+1, tox(T) }. % at most a 2nd cohort allowed
%% trial --> 
%% trial --> [esc].
%% trial, [sty] --> trial.
%% trial, [esc] --> trial.
%% esc --> [1 ^ T], { tox(T) }.
%% esc --> esc | sty.
%% sty --> esc.
%% %% esc --> [1^T], { tox(T) }.
%% esc, [2^T, 1^0] --> [1^0], { tox(T) }, esc.
%% esc, [D^T, D_^0] --> [D_^0], { dose(D_), D is D_ + 1, tox(T) }. %%, (esc | sty).
%% sty, [D^T, D^1] --> [D^1], { dose(D), tox(T) }, esc.
%% %% esc --> [D^T], [ , [1^T] --> [].

%% ntox(T) :- between(0, 3, T).
%% esc, [D ^ T], 
%% esc --> [1^T], ntox(T).
%% esc --> notox

%% Try the Fibonacci exercise from
%% https://web.archive.org/web/20090503181751/http://www.cotilliongroup.com/arts/DCG.html

% So fibonacci//2 is actually fibonacci/4, with the 1st argument X threading the
% desired Nth number in the form fib(N,Fn), and the 2nd argument being a alternately
% a 2-element list that contains the first 2 elements at the head of the results list,
% or else a whole list with as-yet undetermined tail.
% The 'hidden' 3rd & 4th arguments are the lists managed implicitly by DCG notation.

% If the sought X=fib(N,Fn) can be unified with head of results list, then stop.
fibonacci(X, [X|_]) --> []. % hmm! this *terminates* the recursion

% In the following clauses, X is merely threaded thru.
% The substance is the stuff after the "X,".

% Using a time-series-inspired lag-type notation, I'll write 4 successive Fibonacci numbers
% as F_2, F_1, F0, F1.
% This clause is about taking the last 2 memoized elements of sequence (F_1, F_2),
% and generating a new pair F0 := F_2 + F_1, then F1 := F_1 + F).
fibonacci(X, [F_2, F_1]) -->
    next_fib(X, F_2, F_1, F0),
    next_fib(X, F_1, F0, F1),
    fibonacci(X, [F0, F1]).

% Alternately, the 
fibonacci(X, [A, B | Rest]) -->
    succ_push(X, A),
    fibonacci(X, [B | Rest]). % recur to 1st clause to check if X unifies with B

next_fib(X, fib(Idx1, Num1), fib(Idx2, Num2), fib(Idx3, Num3)) -->
    succ_push(X, fib(Idx1, Num1)),
    {
	Idx3 is Idx2 + 1,
	Num3 is Num2 + Num1
    }.

% I suppose that we could also check for NumX > NumA when NumX is instantiated.
% (I think in this program NumA will always be instantiated.)
gt(fib(IndX,NumX), fib(IndA,NumA)) :-
    IndX > IndA.

succ_push(X, A, Old, New) :-
    gt(X, A), % this fails if A is the desired answer
    concat([Old|T1] - T1, [B|T2] - T2, New).

concat(A1 - Z1, Z1 - Z2, A1 - Z2).

% simple_list//1 is actually simple_list/3
simple_list([]) --> [].
simple_list([H|Diff] - Diff, L0, List) :-
    
simple_list([H, Datum|Diff] - Diff, L0, List) :-
    simple_list(H, [Datum|L0], List).
