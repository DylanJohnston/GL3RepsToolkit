# (README last updated: 8th Septemeber 2019)

# A toolkit for the Representation Theory of GL3(Fp)

This repository contains my python code for a GL3 Representation theory toolkit. 
Within this toolkit you are able to (assuming this was accessed via my paper so terminology is established.):
 ~ Decompose a W representation into a linear combination of (not necessarily irreducible) F reps. 
 ~ Perform Steinberg's tensor theorem on a F rep
 ~ Perform the Littlewood-Richardson rule on two W reps
 ~ Write W as a linear combination of irreducible representations (only works for W in p^2-restricted region)
 ~ Write FxF as a linear combination of irreducible representations (only works for both weights in the p-restricted region)
 
**Both methods LR_coeff_finder() and SLP_coeff_finder() struggle massively for "large" weight values due to floating point issues when solving the nxn linear system. You will most likely hit maximum recursion depth if your weights are too big. I don't think an answer will be returned unless it's pretty much certain to be correct. I think getting no answer is better than getting a wrong answer back! I will have to read up on better ways to solve this system. So as for now please only use the LR option for both weights less than (6,6,6) or (7,7,7) maybe, and only use any of the decomposition methods (option 1, option 4, option 5) for primes 7 and below.***

This code makes no effort to be pythonic and may be quite slow and inefficient in places. This is something I will fix and work on as time goes on but it was not in my highest interest to have 'nice' code, just code that worked.

Using this is simple: just download the .py file and open it and compile in Python. You will be greeted by a menu and prompted for inputs within the shell.

-Dylan.
