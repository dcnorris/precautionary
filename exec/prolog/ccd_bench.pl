%% Benchmarking for ccd module

:- use_module(ccd).

%% Because ccd:regression FAILS to indicate that regression has not occurred,
%% we employ DISJUNCTION (;) in order to run this test as a benchmark, under
%% variously toggled goal expansions.
%%
%% Calls to consult/1 have CUMULATIVE side-effects, despite the disjunction,
%% so only one selected 'elegance' may be restored before restarting session
%% as e.g. with C-u C-u [F10] in ediprolog.
benchmark(Elegance) :-
    (	format("Pure code that will go in the paper:~n", []),
	ccd:regression
    ;	format("Faster code, exploiting several inelegances:~n", []),
	consult(inelegance),
	use_module(ccd), % recompile
	ccd:regression
    ;	format("Partial restoration of elegance ... [~q].~n", [Elegance]),
	consult(Elegance), % restore just this elegance
	consult(inelegance),
	use_module(ccd), % recompile
	ccd:regression
    ).

version(Version) :-
    '$scryer_prolog_version'(Version).

%?- version(Version).
%@    Version = "v0.8.123-595-g176da ...".

%?- benchmark(elegant_if_). % remember to use C-u C-u F10
%@ Pure code that will go in the paper:
%@  D = 1 ...   % CPU time: 0.866 seconds
%@    % CPU time: 0.870 seconds
%@  J(1) = 20.
%@  D = 2 ...   % CPU time: 10.524 seconds
%@    % CPU time: 10.528 seconds
%@  J(2) = 212.
%@ Faster code, exploiting several inelegances:
%@ Warning: overwriting goal_expansion/2
%@  D = 1 ...   % CPU time: 1.151 seconds
%@    % CPU time: 1.155 seconds
%@  J(1) = 20.
%@  D = 2 ...   % CPU time: 13.564 seconds
%@    % CPU time: 13.568 seconds
%@  J(2) = 212.
%@ Partial restoration of elegance ... [elegant_if_].
%@ Warning: overwriting goal_expansion/2
%@  D = 1 ...   % CPU time: 1.207 seconds
%@    % CPU time: 1.211 seconds
%@  J(1) = 20.
%@  D = 2 ...   % CPU time: 14.263 seconds
%@    % CPU time: 14.267 seconds
%@  J(2) = 212.
%@ false.

%?- benchmark(elegant_qcompare). % remember to use C-u C-u F10
%@ Pure code that will go in the paper:
%@  D = 1 ...   % CPU time: 0.924 seconds
%@    % CPU time: 0.928 seconds
%@  J(1) = 20.
%@  D = 2 ...   % CPU time: 11.387 seconds
%@    % CPU time: 11.391 seconds
%@  J(2) = 212.
%@ Faster code, exploiting several inelegances:
%@ Warning: overwriting goal_expansion/2
%@  D = 1 ...   % CPU time: 1.217 seconds
%@    % CPU time: 1.222 seconds
%@  J(1) = 20.
%@  D = 2 ...   % CPU time: 14.498 seconds
%@    % CPU time: 14.503 seconds
%@  J(2) = 212.
%@ Partial restoration of elegance ... [elegant_qcompare].
%@ Warning: overwriting goal_expansion/2
%@  D = 1 ...   % CPU time: 0.870 seconds
%@    % CPU time: 0.874 seconds
%@  J(1) = 20.
%@  D = 2 ...   % CPU time: 10.730 seconds
%@    % CPU time: 10.734 seconds
%@  J(2) = 212.
%@ false.

%% ------------------------------------------------------------
%%
%% The following benchmarks were done with older version (my bad!)

%% With these definitions, I get:
%?- version(Version).
%@    Version = "v0.8.123-587-g87ef3 ...".

%?- benchmark(elegant_if_). % remember to use C-u C-u F10
%@ Pure code that will go in the paper:
%@  D = 1 ...   % CPU time: 11.179 seconds
%@    % CPU time: 11.183 seconds
%@  J(1) = 20.
%@  D = 2 ...   % CPU time: 120.674 seconds
%@    % CPU time: 120.678 seconds
%@  J(2) = 212.
%@ Faster code, exploiting several inelegances:
%@ Warning: overwriting goal_expansion/2
%@  D = 1 ...   % CPU time: 1.034 seconds
%@    % CPU time: 1.039 seconds
%@  J(1) = 20.
%@  D = 2 ...   % CPU time: 12.301 seconds
%@    % CPU time: 12.305 seconds
%@  J(2) = 212.
%@ Partial restoration of elegance ... [elegant_if_].
%@ Warning: overwriting goal_expansion/2
%@  D = 1 ...   % CPU time: 2.233 seconds
%@    % CPU time: 2.239 seconds
%@  J(1) = 20.
%@  D = 2 ...   % CPU time: 24.493 seconds
%@    % CPU time: 24.498 seconds
%@  J(2) = 212.
%@ false.
%@ Pure code that will go in the paper:
%@  D = 1 ...   % CPU time: 11.794 seconds
%@    % CPU time: 11.800 seconds
%@  J(1) = 20.
%@ Faster code, exploiting several inelegances:
%@ Warning: overwriting goal_expansion/2
%@  D = 1 ...   % CPU time: 1.128 seconds
%@    % CPU time: 1.132 seconds
%@  J(1) = 20.
%@ Partial restoration of elegance ... [elegant_if_].
%@ Warning: overwriting goal_expansion/2
%@  D = 1 ...   % CPU time: 2.260 seconds
%@    % CPU time: 2.264 seconds
%@  J(1) = 20.
%@ false.
%@ Pure code that will go in the paper:
%@  D = 1 ...   % CPU time: 10.995 seconds
%@    % CPU time: 11.000 seconds
%@  J(1) = 20.
%@ Faster code, exploiting several inelegances:
%@ Warning: overwriting goal_expansion/2
%@  D = 1 ...   % CPU time: 1.035 seconds
%@    % CPU time: 1.039 seconds
%@  J(1) = 20.
%@ Partial restoration of elegance ... [elegant_if_].
%@ Warning: overwriting goal_expansion/2
%@  D = 1 ...   % CPU time: 2.083 seconds
%@    % CPU time: 2.087 seconds
%@  J(1) = 20.
%@ false.

%?- benchmark(elegant_qcompare). % remember to use C-u C-u F10
%@ Pure code that will go in the paper:
%@  D = 1 ...   % CPU time: 11.188 seconds
%@    % CPU time: 11.194 seconds
%@  J(1) = 20.
%@  D = 2 ...   % CPU time: 120.373 seconds
%@    % CPU time: 120.378 seconds
%@  J(2) = 212.
%@ Faster code, exploiting several inelegances:
%@ Warning: overwriting goal_expansion/2
%@  D = 1 ...   % CPU time: 0.999 seconds
%@    % CPU time: 1.003 seconds
%@  J(1) = 20.
%@  D = 2 ...   % CPU time: 12.227 seconds
%@    % CPU time: 12.231 seconds
%@  J(2) = 212.
%@ Partial restoration of elegance ... [elegant_qcompare].
%@ Warning: overwriting goal_expansion/2
%@  D = 1 ...   % CPU time: 9.589 seconds
%@    % CPU time: 9.594 seconds
%@  J(1) = 20.
%@  D = 2 ...   % CPU time: 107.605 seconds
%@    % CPU time: 107.610 seconds
%@  J(2) = 212.
%@ false.
%@ Pure code that will go in the paper:
%@  D = 1 ...   % CPU time: 11.280 seconds
%@    % CPU time: 11.285 seconds
%@  J(1) = 20.
%@ Faster code, exploiting several inelegances:
%@ Warning: overwriting goal_expansion/2
%@  D = 1 ...   % CPU time: 1.084 seconds
%@    % CPU time: 1.090 seconds
%@  J(1) = 20.
%@ Partial restoration of elegance ... [elegant_qcompare].
%@ Warning: overwriting goal_expansion/2
%@  D = 1 ...   % CPU time: 10.436 seconds
%@    % CPU time: 10.440 seconds
%@  J(1) = 20.
%@ false.
%@ Pure code that will go in the paper:
%@  D = 1 ...   % CPU time: 11.039 seconds
%@    % CPU time: 11.044 seconds
%@  J(1) = 20.
%@ Faster code, exploiting several inelegances:
%@ Warning: overwriting goal_expansion/2
%@  D = 1 ...   % CPU time: 1.043 seconds
%@    % CPU time: 1.048 seconds
%@  J(1) = 20.
%@ Partial restoration of elegance ... [elegant_qcompare].
%@ Warning: overwriting goal_expansion/2
%@  D = 1 ...   % CPU time: 9.901 seconds
%@    % CPU time: 9.905 seconds
%@  J(1) = 20.
%@ false.
