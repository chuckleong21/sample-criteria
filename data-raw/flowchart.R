## code to prepare `flowchart` grViz object goes here
library(DiagrammeR)
library(DiagrammeRsvg)
library(xml2)
library(magrittr)
library(rsvg)

add_mathjax(
  grViz(
    "
  digraph G {
    # graph statements
    graph [label = '$\\\\Large{R^2_{CS,adj,k,r}\\\\hspace{0.5em}\\\\textrm{derivation}}$',
           fontsize = 10,
           labelloc = 't'];

    # node statements
    node [shape = diamond]
    a [label = '$\\R ^2_{CS,adj,k,r}$\navailable?', height = 1];
    b [label = '$\\R^2_{CS_app,k,r}$\navailable?', height = 1];
    c [label = '(pairwise) C-statistics\navailable?'];
    f [label = 'Adjust?'];

    node [shape = box]
    d [label = '$\\\\phi_{k,r}$\napproach'];
    e [label = 'use directly'];
    g [label = 'simulation\napproach'];
    h [label = '$\\\\texttt{max}(\\R^2_{CS,app,k,r})$'];
    i [label = '$\\R^2_{CS,app,k,r}$'];
    j [label = '$\\R^2_{Nagelkerke}$\n$\\\\times{}$\n$\\\\texttt{max}(R^2_{CS,app,k,r})$'];
    k [label = '$\\R^2_{Nagelkerke}$\n$\\\\times{}$\n$\\\\texttt{max}(R^2_{CS,app,k,r})$'];

    # edge statements
    a -> e [label = 'yes'];
    a -> b [label = 'no'];
    b -> c [label = 'no'];
    c -> d [label = 'no'];
    b -> f [label = 'yes'];
    c -> g [label = 'yes'];
    d -> h; h ->k;
    f -> i [label = 'no'];
    f -> j [label = 'yes'];
  }
  "
  )
) -> flowchart

usethis::use_data(flowchart, internal = TRUE)
