24-10-17
===
---
Review
---
Fitting both the gradient and intercept of our $$S_2$$ curves fits where the intercept was not consistent with the theoretical value of $$\gamma$$. For recap:

| m        | Theoretical $$\gamma$$           | Actual intercept  |
| ------------- |:-------------:| -----:|
| 1      | 0 | 0.199 $$\pm$$ 0.008 |
| 3      | 0.549 | -0.267 $$\pm$$ 0.014 |
| 5     | 0.805      |   -0.583 $$\pm$$ 0.034 |

A suggestion was put forward to run a new linear regression, but with the intercept fixed to the correct theoretical values above.
Surprisingly, performing linear regression in Python using the `statsmodels` or `scipy.stats` packages proved non trivial to implement. Instead, we can rearrange our equation:

$$S_2 = a \sqrt n - \gamma$$ 
$$\rightarrow$$ $$S_2 +\gamma = a \sqrt n$$

Below are the plots for $$m=1,3,5$$. We expect an intercept with the origin if our results match the theory:

https://raw.githubusercontent.com/bjader/qhe/master/images/m1_24-10-17.png
https://raw.githubusercontent.com/bjader/qhe/master/images/m3_24-10-17.png
https://raw.githubusercontent.com/bjader/qhe/master/images/m5_24-10-17.png

This is fantastic! From first observations, our data seems to fit the theoretical scaling very well. As a reminder, we arrived at the theoretical equation using the topological area law:

$$S_{2} = al_{a} - \gamma$$ (https://arxiv.org/pdf/1106.0015.pdf)

and substituting the relationship obtained from our desnity profiles

$$l_{a} = 2r_{0} = 2\sqrt{n}$$ 


Future goals
---
Supposing we have used Monte Carlo to verify the field theoretic value of Topological Entanglement Entropy $$\gamma$$, this is most similar to the results published in the Zhang paper linked above. However, I don't fully understand their use of "chiral" and "lattice" states. Furthermore they only compute m=1,3,4 with errors of 5-10%.

Proposed steps:
- Generate full statistical properties of the fit. All I have done is eyeballed a "good fit" at the moment
        - Compute the p-value of our fit
        - Allow intercept to vary again to get an uncertainty on $$\gamma$$
- Understand the differences and similarities between Zhang and our work