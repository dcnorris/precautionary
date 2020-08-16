// initialize an introjs instance          
var intro = introJs();

/*
 * Hard-coding the help system as below avoids complexities of R -> Javascript
 * communication, as well as awkwardness of typing in a data.frame representation.
 * In theory, however, we do lose the ability to customize help text on-the-fly,
 * which might be useful for commenting on particular outputs in the safety table.
*/
intro.setOptions({steps: [
  {
    element: '#hyperprior',
    intro: "TODO: Explain the hyperprior visualization",
    position: 'bottom'
  },
  {
    element: '#simprogress',
    intro: "TODO: Explain the progress bar, with advice to stop sim at ±0.1 MCSE",
    position: 'top'
  }
]});

// TODO: Understand why this doesn't get invoked properly from server
// handler 1
Shiny.addCustomMessageHandler("setHelpContent",
  
  // callback function. 
  // note: data is passed by shiny and contains the tour data
  function(data){

    // load data 
    intro.setOptions({steps: data});
  }
);

// handler 2
Shiny.addCustomMessageHandler("startHelp", function(message) {

    // start intro.js
    // note: we don't need information from shiny, just start introJS
    intro.start();
  }
);
