# (README last updated: 4th Septemeber 2019)

# A toolkit for the Representation Theory of GL3(Fp)  (README last updated: 4th Septemeber 2019)

This repository contains my python code for a GL3 Representation theory toolkit. 
Within this toolkit you are able to (assuming this was accessed via my paper so terminology is established.):
 ~ Decompose a W representation into a linear combination of (not necessarily irreducible) F reps. 
 ~ Perform Steinberg's tensor theorem on a F rep
 ~ Perform the Littlewood-Richardson rule on two W reps (warning: this is slow, as of 4/9/19)
 ~ Write W as a linear combination of irreducible representations (only works for W in p^2-restricted region)
 ~ Write FxF as a linear combination of irreducible representations (only works for both weights in the p-restricted region)
 
This code makes no effort to be pythonic and may be quite slow and inefficient in places. This is something I will fix and work on as time goes on but it was not in my highest interest to have 'nice' code, just code that worked.

Using this is simple: just download the .py file and open it and compile in Python. You will be greeted by a menu and prompted for inputs within the shell.

-Dylan.
