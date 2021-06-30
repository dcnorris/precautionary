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
    element: '#dose-levels',
    intro: "This app assumes your dose-escalation study uses a fixed set of 3&ndash;7 doses. You can span the range between lowest and highest doses using an arithmetic or geometric sequence, or manually customize the doses.",
    position: 'bottom'
  },
  {
    element: '#optimal-dose-heterogeneity',
    intro: "A crucial perspective adopted by this app is that the optimal dose of the drug will vary from patient to patient. This app models this variability as a log-normally distributed &lsquo;MTD<sub>i</sub>&rsquo;. You are asked to provide ...",
    position: 'top'
  },
  {
    element: document.querySelector('#median_mtd').parentElement,
    intro: "your best guess about the <i>median</i> MTD<sub>i</sub> in the population, i.e. the dose that will cause a dose-limiting toxicity (DLT) in <i>half</i> of the patient population; ...",
    position: 'bottom'
  },
  {
    element: document.querySelector('#sigma_median').parentElement,
    intro: "an indication of your <i>uncertainty</i> about that median, expressed as a percentage &plusmn; ...",
    position: 'bottom'
  },
  {
    element: document.querySelector('#sigma_CV').parentElement,
    intro: "and your guess as to the <i>coefficient of variation</i> (CV) of MTD<sub>i</sub> within the population.",
//    "<br/><small>For reasons outlined in the <a href=https://cran.r-project.org/web/packages/precautionary/vignettes/Intro.html>introductory vignette</a> for R package 'precautionary', you don't need to further specify your uncertainty about &sigma;<sub>CV</sub>.</small>",
    position: 'bottom'
  },
  {
    element: '#hyperprior',
    intro: "You see here 1000 possible MTD<sub>i</sub> distributions, drawn randomly from the <i>hyperprior</i> you define via <b>median MTD<sub>i</sub></b>, <b>&sigma;<sub>median</sub></b> and <b>&sigma;<sub>CV</sub></b>. To gain a better sense for the meaning of these parameters, experiment with different values while observing how this plot shifts and stretches.",
    position: 'bottom'
  },
  {
    element: '#resample',
    intro: "You can also click this button to resample the hyperprior, and check the stability of the results.",
    position: 'right'
  },
  {
    element: '#dose-escalation-design',
    intro: "You can choose from several popular dose-finding designs, and input values for applicable design parameters.",
    position: 'right'
  },
  {
    element: '#crm-skeleton',
    intro: "If you choose the CRM, you can enter its skeleton here, or accept defaults obtained from the hyperprior sample",
    position: 'right'
  },
  {
    element: document.querySelector('#r0').parentElement,
    intro: "Even when your dose-finding design recognizes only <i>binary</i> DLT's, you can still anticipate <i>graded</i> toxicities provided you are willing to venture a guess as to your drug's <i>therapeutic index</i> (TI). In a dose-finding trial conducted under an &lsquo;MTD heuristic&rsquo;, the most useful TI notion is a <i>dosing ratio</i> <b>r<sub>0</sub></b> that relates adjacent high-level toxicities. Ask yourself, <i><b>Typically, what dose multiplier would push a grade-3 toxicity up to grade 4, or a grade-4 toxicity up to grade 5?</b></i> (Example: If you expect that a 50% increase in dose would typically bump the toxicity grade up by 1 level, then set <span>r<sub>0</sub> = 1.5</span>.)",
    position: 'top'
  },
  {
    element: '#safety',
    intro: "This simple table contains the key safety-related quantities of interest, as estimated from your trial simulations together with <b>r<sub>0</sub></b>. Expected total enrollment is broken down into expected numbers of patients who will experience each grade of toxicity. Obviously, the numbers of patients expected to experience grade-4 and grade-5 toxicities are of primary concern from a safety standpoint.",
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
