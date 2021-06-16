:- module(inconceivable, [
	      inconceivable/3
	  ]).

:- use_module(library(clpz)).
:- use_module(library(format)).
:- use_module(library(time)).
:- use_module(library(charsio)).

%% TODO: Refine the console output from inconceivable/3, with an understanding
%%       that it will generally be called on a Query that MUST FAIL.
%%       There may even be scope for omitting the separate Var argument
%%       by abstracting it automatically from the Query itself, recognized
%%       perhaps as the only unbound named variable.
:- meta_predicate inconceivable(0, -, +).
inconceivable(Query, Var, Range) :-
    %% TODO: Is there a way to get the name of Var as originally provided,
    %%       instead of the renamed version?
    write_term_to_chars(Var, [], IndexName), % at present, IndexName invariably "A"
    Var in Range,
    indomain(Var),
    format(" % ~s = ~d ...", [IndexName, Var]),
    time(call(Query)).
