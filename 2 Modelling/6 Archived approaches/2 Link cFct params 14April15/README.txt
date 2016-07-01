The intention for this was to force parameters that applied to cF and ct to be linked between the two. (e.g. same betas, same epsilons etc). Intention:  Linked parameters=  b, p, j, m, e; NOT linked=i, f, w


BUT, given the number of permutations of w for cF vs ct, it is too computationally intensive to grid them all. As a check, I picked good-performing (top 10, BIC-wise) models from the previous approach (that didn't force the cFct link), and ran these models with the forced cFct link, to see if this new approach would deteriorate the BICs. Ls would likely drop (since the forced link is a constraint), but the BIC penalty would be substantially less as well. 

Models were picked that had the same base (i.e. bpj vs bpm), that either had the same or different w parameters for the cF vs ct. In both cases, the BICs deteriorate so badly that it's not worth considering (see picture BIC results). 

As such, this approach is aborted.


D:\Dropbox\SANDISK\4 Explore experiment\3 Analysis\4 Fit computational models\6 Archived approaches\1 Params within condition 12April15

