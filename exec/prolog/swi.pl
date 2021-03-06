:- use_module(library(clpfd)).

% Approximate format_//2 as provided by Scryer's library(format):
format_(Fs, Args) --> call(format__(Fs, Args)).

format__(Fs, Args, Ls0, Ls) :-
    format(chars(Ls0,Ls), Fs, Args).

