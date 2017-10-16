16-10-17
===
---
Review
---
At large particle numbers, the scaling of $$S_2$$ appeared to tail off for $$m=3,5$$. This is inconsistent with our expectation of that it would obey the relationship $$S_2 = a \sqrt n - \gamma$$

A solution was put forward that the size of the total system ($$L$$) was the cause.

Rerunning the results for $$L=20$$ and comparing with the previous $$L=10$$ results:

https://raw.githubusercontent.com/bjader/qhe/master/images/m3_16-10-17.png
https://raw.githubusercontent.com/bjader/qhe/master/images/m5_16-10-17.png

From this we can conclude that the system size was limiting our results previously.

Future goals
---
Unfortunately, conducting linear regression on our new results still does not align with the expected value of $$\gamma=\frac{1}{2} \ln(m)$$. 

| m        | Theoretical $$\gamma$$           | Actual  |
| ------------- |:-------------:| -----:|
| 1      | 0 | 0.199 $$\pm$$ 0.008 |
| 3      | 0.549 | -0.267 $$\pm$$ 0.014 |
| 5     | 0.805      |   -0.583 $$\pm$$ 0.034 |

A previous suggestion for this discrepency was that our scaling understanding $$S_2 = a \sqrt n - \gamma$$ could be fundamentally flawed. Perhaps this needs to be investigated. 