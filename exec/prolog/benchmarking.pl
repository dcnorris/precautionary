% Module of benchmarking utilities

:- module(benchmarking, [
	      version/1,
	      time/2
	  ]).

version(Version) :-
    '$scryer_prolog_version'(Version).

%% The following adapted from @triska's time:time/1
:- meta_predicate time(0, -).

time(Goal, Elapsed) :-
        '$cpu_now'(T0),
        (   Goal,
            since_elapsed(T0, Elapsed)
	;   since_elapsed(T0, Elapsed)
	).

since_elapsed(T0, Elapsed) :-
        '$cpu_now'(T),
        Elapsed is T - T0.
