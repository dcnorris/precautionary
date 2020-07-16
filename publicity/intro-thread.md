---
title: "intro-thread"
author: "David C. Norris"
date: "7/15/2020"
output:
  html_document:
    keep_md: true
---

<!--
My immediate aim here is to generate an HTML document
suitable for copy-n-paste into Tweetdeck or Twitter.

I regard this as an initial step toward a full-fledged
output format for sciencey threads on Twitter.
-->



New #rstats package #precautionary on CRAN

"The shaft of the arrow had been feathered with one of the eagle's own plumes. We often give our enemies the means of our own destruction."

 --- Aesop

[THREAD] 1/

https://CRAN.R-project.org/package=precautionary

---


By unifying several main dose-escalation designs under a single simulation framework, @CatchTwentyToo's valuable #escalation package has greatly facilitated the development of package #precautionary. 2/

https://twitter.com/CatchTwentyToo/status/1229444639238873089

---


My initial experiments toward #precautionary began by cloning the source for #escalation, and "changing just one line of simulation code" precisely as I intimated in this 2018 letter desk-rejected by @CCR_AACR: https://pubpeer.com/publications/4CEB4434EA94CC0AA3BCAF63617A78#1
 3/

<img src="big-git-diff.png" width="1200" />

---


In R, 'rbinom(size=1, prob=.)' generates a random Bernoulli variate (0 or 1) with the given prob of being a 1. Under the hood, this function almost surely works by (a) drawing a uniformly distributed random number between 0 and 1, then (b) comparing that against the given prob 4/



---


Thus, from a purely *software-engineering* POV, these changes 'unwrap' rbinom() to expose its inner workings, retaining $u_i \in (0,1]$.

But these changes to the underlying formal workings turn out to have *conceptual* correlates. Retaining $u_i$ makes a world of difference! 5/



---


This $u_i$ has been on my radar for quite some time. I believe I've traced it back to its origins nearly 2 decades ago.

From the look of things, the #OneSizeFitsAllogists knew from the outset how dangerous it is. 6/

https://twitter.com/davidcnorrismd/status/1271804500664365056

---


What makes $u_i$ so dangerous is that it acknowledges a latent toxicity threshold characterizing an individual patient's susceptibility or tolerance to the drug's toxicity. No wonder they tried so hard to keep the lid on the Pandora's Box of rbinom()! 7/



---


But when you can't even get the Opening up the Pandora's Box of rbinom() unleashes a lintany of ills I can also include *multiple* links: 8/

https://precisionmethods.guru/

https://www.gnu.org

---


This is slightly obnoxious, that I can't use continuation lines in chunks. Oh, well. I wonder how long it will take to use up my 280 character limit! One way to get there faster, I suppose, is to include a link. Let me try that: 9/

https://precisionmethods.guru

---

