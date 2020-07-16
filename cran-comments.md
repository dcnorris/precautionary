## Test environments
* local OS X install, R 4.0.1
* win-builder (devel)

## R CMD check --as-cran results
There were no ERRORs or WARNINGs

There were no NOTEs.

## Responses to July 16 feedback on initial CRAN submission

1. Added a key reference to Description field of DESCRIPTION file

2. Reduced unnecessary exports of S3 methods not for use by user:
 - num_doses.*
 - print.(hyper|precautionary)
 - prob_administer.precautionary
 - prob_recommend.precautionary

3. Added documentation for summary.precautionary

4. Added examples for plot,mtdi-(distribution|generator) methods

5. Restore old <- options() where altered in examples & vignette

6. Added "cph" role in Authors@R field of DESCRIPTION
