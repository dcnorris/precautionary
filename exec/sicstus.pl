/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Setup file to run esc.pl with SICStus Prolog.

   SICStus Prolog is a state-of-the-art, ISO standard compliant
   Prolog development system, available from: https://sicstus.sics.se/
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

:- use_module(library(clpfd)).
:- use_module(library(lists)).

% Approximate format_//2 as provided by Scryer's library(format):
:- use_module(library(codesio)).

format_(Fs, Args) --> call(format__(Fs, Args)).

format__(Fs, Args, Ls0, Ls) :-
    with_output_to_codes(format(Fs, Args), Codes, []),
    atom_codes(Atom, Codes),
    atom_chars(Atom, Chars),
    append(Chars, Ls, Ls0).
