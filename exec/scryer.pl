/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Setup file to run esc.pl with Scryer Prolog.

   Scryer Prolog is a modern Prolog system implemented in Rust,
   a language that provides desirable safety guarantees and performance.

   Scryer Prolog aims for compliance with the Prolog ISO standard.
   It is available from https://github.com/mthom/scryer-prolog/
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

:- use_module(library(format)).
:- use_module(library(lists)).
:- use_module(library(dcgs)).
:- use_module(library(clpz)).
:- use_module(library(time)).

nth1(N, List, Elt) :- nth0(N, [elt0|List], Elt).
