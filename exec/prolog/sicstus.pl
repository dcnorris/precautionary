:- use_module(library(clpfd)).
:- use_module(library(lists)).

% Approximate format_//2 as provided by Scryer's library(format):
:- use_module(library(codesio)).

format_(Fs, Args) --> call(format__(Fs, Args)).

format__(Fs, Args, Ls0, Ls) :-
    with_output_to_codes(format(Fs, Args), Codes0, Ls),
    name(Atom0, Codes0),
    atom_chars(Atom0, Ls0).
