%% Exploration of non-monotonicity in CLP(ℤ)

:- use_module(library(clpz)).
:- use_module(library(error)).

n_factorial(0, 1).
n_factorial(N, F) :-
    N #> 0,
    N1 #= N - 1,
    F #= N * F1,
    n_factorial(N1, F1).

%% A quick note about termination:

%?- n_factorial(10, F).
%@    F = 3628800
%@ ;  false.

%?- time(n_factorial(N, 3628800)).
%@    % CPU time: 2.179 seconds
%@    N = 10
%@ ;  % CPU time: 2.203 seconds
%@    false. % <-- NB: having found 1 solution, the program knows immediately there are none further

%% Even if we ask in cases where no solution exists, the program DOES terminate:
%?- time((F #= 3628801, n_factorial(N, F))).
%@    % CPU time: 79.779 seconds
%@ false.

%% But we could offer some help to termination...
%?- n_factorial(N, F), N^4 #> F.
%@    N = 2, F = 2
%@ ;  N = 3, F = 6
%@ ;  N = 4, F = 24
%@ ;  N = 5, F = 120
%@ ;  N = 6, F = 720
%@ ;  <hangs>
%@ caught: error('$interrupt_thrown',repl)

%% Our implementation doesn't 'know' (and can't find out)
%% that N! > N^4 ∀ N > 6. But by supplying this extra info
%% we can speed up termination:

%?- time((N #< 7 #\/ N^4 #< F, F #= 3628801, n_factorial(N, F))).
%@    % CPU time: 0.473 seconds
%@ false.

%% That said, let's exhibit what n_factorial/2 knows...

%% Although the ordinary programmer might think this is 'pretty cool',
%% for the logic programmer this is perhaps closer to the bare minimum:
%% "Tell me everything that holds."

%?- n_factorial(N, F).
%@    N = 0, F = 1
%@ ;  N = 1, F = 1
%@ ;  N = 2, F = 2
%@ ;  N = 3, F = 6
%@ ;  N = 4, F = 24
%@ ;  N = 5, F = 120
%@ ;  N = 6, F = 720
%@ ;  ...

%% But check THIS out! Does this not seem slightly magical?

%?- n_factorial(N+2, F-1).
%@    N = -1, F = 2
%@ ;  N = 0, F = 3
%@ ;  N = 1, F = 7
%@ ;  N = 2, F = 25
%@ ;  N = 3, F = 121
%@ ;  N = 4, F = 721
%@ ;  N = 5, F = 5041
%@ ;  ...

%% But WAIT A MINUTE! Where's my N = -2, F = 2 answer?
%?- N #= -2, F #= 2, n_factorial(N+2, F-1).
%@ false.

%?- n_factorial(-2+2, 2-1). % <-- very strong 'hints' indeed!
%@ false.

%% Even very strong hints don't help elicit this solution! What went wrong here?
%% Should I have implemented the base case in some more general fashion?

n_fact(Zero, One) :-
    Zero #= 0,
    One #= 1.
n_fact(N, F) :- % identical to general clause of n_factorial/2 above
    N #> 0,
    N1 #= N - 1,
    F #= N * F1,
    n_fact(N1, F1).

%?- n_fact(N+2, F-1).
%@    N = -2, F = 2 % PHEW! We got it this time.
%@ ;  N = -1, F = 2
%@ ;  N = 0, F = 3
%@ ;  N = 1, F = 7
%@ ;  N = 2, F = 25
%@ ;  N = 3, F = 121
%@ ;  N = 4, F = 721
%@ ;  ...

%% All the above was done with non-monotonic clpz:
%?- clpz:monotonic.
%@ false.

%% So what happens if I now demand monotonicity?
%?- assertz(clpz:monotonic).
%@    true.

%?- n_factorial(N, F).
%@    N = 0, F = 1
%@ ;  caught: error(instantiation_error,instantiation_error(unknown(_1384246),1))

%% Oops! I forgot I would have to rewrite to avoid defaulty representation...

m_factorial(0, 1).
m_factorial(N, F) :-
    #N #> 0,
    #N1 #= #N - 1,
    #F #= #N * #F1,
    m_factorial(N1, F1).

%?- m_factorial(N, F).
%@    N = 0, F = 1
%@ ;  N = 1, F = 1
%@ ;  N = 2, F = 2
%@ ;  N = 3, F = 6
%@ ;  N = 4, F = 24
%@ ;  ...

%% So monotonicity doesn't ruin our 'bare minimum' requirements of course.
%% But do we still get the magic?

%?- m_factorial(N-1, F).
%@ caught: error(type_error(integer,_1384249-1),must_be/2)

%?- m_factorial(1, F).
%@    F = 1
%@ ;  false.

%?- m_factorial(10, F).
%@    F = 3628800
%@ ;  false.

%?- m_factorial(N, 6).
%@    N = 3
%@ ;  false.

%% Do I dare ...?
%?- m_factorial(N, 3628800).
%@    N = 10
%@ ;  false.

%?- m_factorial(N, 3628801).
%@ false.

%% CONCLUSION: Yes, monotonicity does ruin the 'magic' of Prolog.
%% COROLLARY: _BUT_ 'magic' might not be such a great thing!
