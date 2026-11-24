# Chapter 10: Inferences Based on Two Samples

*Bilkent University MATH 260 lecture notes, Spring 2026. [Back to contents](../README.md)*

## 10.1 The Two-Sample $z$ CI and Test

**Assumptions.**

1.  $X_1,\dots,X_m$ is a random sample from a population with mean $\mu_1$ and standard deviation $\sigma_1$.

2.  $Y_1,\dots,Y_n$ is a random sample from a population with mean $\mu_2$ and standard deviation $\sigma_2$.

3.  The $X$ and $Y$ samples are independent of each other.

Consider the parameter $\mu_1-\mu_2$. The natural estimator of $\mu_1-\mu_2$ is $\bar{X}-\bar{Y}$; the test statistic is the standardized form of $\bar{X}-\bar{Y}$. We have

$$
E(\bar{X}-\bar{Y}) = E(\bar{X})-E(\bar{Y}) = \mu_1-\mu_2,
$$

$$
V(\bar{X}-\bar{Y}) = V(\bar{X})+V(\bar{Y}) = \frac{V(X_1+\dots+X_m)}{m^2} + \frac{V(Y_1+\dots+Y_n)}{n^2} = \frac{\sigma_1^2}{m}+\frac{\sigma_2^2}{n},
$$

$$
\sigma_{\bar{X}-\bar{Y}} = \sqrt{\frac{\sigma_1^2}{m}+\frac{\sigma_2^2}{n}}.
$$

**Proposition.** The expectation of $\bar{X}-\bar{Y}$ is $\mu_1-\mu_2$, so $\bar{X}-\bar{Y}$ is an unbiased estimator of $\mu_1-\mu_2$. Its standard deviation is $\sqrt{\sigma_1^2/m+\sigma_2^2/n}$.

With $\theta=\mu_1-\mu_2$, $\hat\theta=\bar{X}-\bar{Y}$ and $\sigma_{\hat\theta}=\sigma_{\bar{X}-\bar{Y}}$, the test statistic value is

$$
\frac{\hat\theta - \text{null value}}{\sigma_{\hat\theta}} .
$$

**The two-sample $z$ CI for $\mu_1-\mu_2$ with known $\sigma$’s.** Assumptions: (1) the population distributions are normal; (2) $\sigma_1$ and $\sigma_2$ are known. The pivotal quantity is

$$
Z = \frac{\bar{X}-\bar{Y}-(\mu_1-\mu_2)}{\sqrt{\dfrac{\sigma_1^2}{m}+\dfrac{\sigma_2^2}{n}}} \sim N(0,1).
$$

Then we have the probability statement

$$
P\!\left(-z_{\alpha/2} \lt  \frac{\bar{X}-\bar{Y}-(\mu_1-\mu_2)}{\sigma_{\bar{X}-\bar{Y}}} \lt  z_{\alpha/2}\right) = 1-\alpha
$$

$$
\Rightarrow\; P\!\left(\bar{X}-\bar{Y}-z_{\alpha/2}\,\sigma_{\bar{X}-\bar{Y}} \lt  \mu_1-\mu_2 \lt  \bar{X}-\bar{Y}+z_{\alpha/2}\,\sigma_{\bar{X}-\bar{Y}}\right) = 1-\alpha .
$$

This implies the $100(1-\alpha)\%$ CI for $\mu_1-\mu_2$ with known $\sigma$’s: for two independent samples from normal population distributions, the $100(1-\alpha)\%$ CI has endpoints

$$
\bar{x}-\bar{y} \pm z_{\alpha/2}\,\sigma_{\bar{X}-\bar{Y}},
$$

with lower limit $\bar{x}-\bar{y}-z_{\alpha/2}\sigma_{\bar{X}-\bar{Y}}$ and upper limit $\bar{x}-\bar{y}+z_{\alpha/2}\sigma_{\bar{X}-\bar{Y}}$. The corresponding one-sided confidence bounds are $\bar{x}-\bar{y}-z_\alpha\sigma_{\bar{X}-\bar{Y}}$ (lower bound) and $\bar{x}-\bar{y}+z_\alpha\sigma_{\bar{X}-\bar{Y}}$ (upper bound).

**Sample size determination.** For known $\sigma_1$ and $\sigma_2$ we can assume $m=n$. The CI width is

$$
W = 2z_{\alpha/2}\sqrt{\frac{\sigma_1^2}{m}+\frac{\sigma_2^2}{n}}
\;\overset{m=n}{=}\; 2z_{\alpha/2}\sqrt{\frac{\sigma_1^2+\sigma_2^2}{n}}
\;\Rightarrow\; n = \frac{z_{\alpha/2}^2(\sigma_1^2+\sigma_2^2)}{(W/2)^2}.
$$

**Test procedures for $\mu_1-\mu_2$ with known $\sigma$’s (two-sample $z$ test).** Null hypothesis $H_0:\mu_1-\mu_2=\Delta_0$. Test statistic value:

$$
z = \frac{\bar{x}-\bar{y}-\Delta_0}{\sqrt{\dfrac{\sigma_1^2}{m}+\dfrac{\sigma_2^2}{n}}}.
$$

| Alternative hypothesis       | Rejection region                                                | $P$-value        |
|:-----------------------------|:----------------------------------------------------------------|:-----------------|
| $H_a:\mu_1-\mu_2\gt \Delta_0$   | $z\ge z_\alpha$ (upper-tailed)                                  | $1-\Phi(z)$      |
| $H_a:\mu_1-\mu_2\lt \Delta_0$   | $z\le -z_\alpha$ (lower-tailed)                                 | $\Phi(z)$        |
| $H_a:\mu_1-\mu_2\ne\Delta_0$ | either $z\le -z_{\alpha/2}$ or $z\ge z_{\alpha/2}$ (two-tailed) | $2[1-\Phi(\lvert z\rvert)]$ |

These concepts may be applied to nonnormal distributions when $m,n\gt 40$ (by the CLT).

**$\beta$, power, and sample size determination.** Consider $H_a:\mu_1-\mu_2\gt \Delta_0$ with rejection region $z\ge z_\alpha$:

$$
z = \frac{\bar{x}-\bar{y}-\Delta_0}{\sqrt{\sigma_1^2/m+\sigma_2^2/n}} \ge z_\alpha
\;\Leftrightarrow\; \bar{x}-\bar{y}\ge \Delta_0 + z_\alpha\sqrt{\frac{\sigma_1^2}{m}+\frac{\sigma_2^2}{n}} .
$$

So $H_0$ is not rejected when $\bar{x}-\bar{y}\lt \Delta_0+z_\alpha\sigma_{\bar{X}-\bar{Y}}$. Let $\Delta'$ be a particular value of $\mu_1-\mu_2$ that exceeds $\Delta_0$. Then

$$
\begin{aligned}
\beta(\Delta') &= P(\text{type II error}) = P(H_0 \text{ is not rejected when it is false})\\
&= P\!\left(\bar{X}-\bar{Y}\lt \Delta_0+z_\alpha\sigma_{\bar{X}-\bar{Y}} \text{ when } \mu_1-\mu_2=\Delta'\right)\\
&= P\!\left(\frac{\bar{X}-\bar{Y}-\Delta'}{\sigma_{\bar{X}-\bar{Y}}} \lt  \frac{\Delta_0-\Delta'}{\sigma_{\bar{X}-\bar{Y}}}+z_\alpha \text{ when } \mu_1-\mu_2=\Delta'\right)
= \Phi\!\left(z_\alpha - \frac{\Delta'-\Delta_0}{\sigma_{\bar{X}-\bar{Y}}}\right).
\end{aligned}
$$

Similarly, with $\sigma_{\bar{X}-\bar{Y}}=\sqrt{\sigma_1^2/m+\sigma_2^2/n}$:

| Alternative hypothesis       | $\beta(\Delta')$                                                                                                                                                               |
|:-----------------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| $H_a:\mu_1-\mu_2\gt \Delta_0$   | $\Phi\!\left(z_\alpha - \dfrac{\Delta'-\Delta_0}{\sigma_{\bar{X}-\bar{Y}}}\right)$                                                                                             |
| $H_a:\mu_1-\mu_2\lt \Delta_0$   | $1-\Phi\!\left(-z_\alpha - \dfrac{\Delta'-\Delta_0}{\sigma_{\bar{X}-\bar{Y}}}\right)$                                                                                          |
| $H_a:\mu_1-\mu_2\ne\Delta_0$ | $\Phi\!\left(z_{\alpha/2} - \dfrac{\Delta'-\Delta_0}{\sigma_{\bar{X}-\bar{Y}}}\right) - \Phi\!\left(-z_{\alpha/2} - \dfrac{\Delta'-\Delta_0}{\sigma_{\bar{X}-\bar{Y}}}\right)$ |

**Exercise.** Find the sample size for the one-tailed and two-tailed procedures.

- For the one-tailed sample size, use $\beta(\Delta')=\beta$.

- For the two-tailed sample size formula, replace $\alpha$ with $\alpha/2$ in the one-tailed formula.

*Answer* (with $m=n$): $\displaystyle n=\frac{(\sigma_1^2+\sigma_2^2)(z_\alpha+z_\beta)^2}{(\Delta'-\Delta_0)^2}$ for a one-tailed test, and the same with $z_{\alpha/2}$ in place of $z_\alpha$ for a two-tailed test.

## 10.2 The Two-Sample $t$ CI and Test

**Assumptions.** (1) The population distributions are normal; (2) $\sigma_1$ and $\sigma_2$ are unknown.

Previously we introduced the random variable

$$
Z = \frac{\bar{X}-\bar{Y}-(\mu_1-\mu_2)}{\sqrt{\dfrac{\sigma_1^2}{m}+\dfrac{\sigma_2^2}{n}}}.
$$

To get rid of the unknown $\sigma_1$ and $\sigma_2$ above, we can replace them with their sample standard deviations $S_1$ and $S_2$.

**Welch’s theorem.** For two normal populations, the random variable

$$
T = \frac{\bar{X}-\bar{Y}-(\mu_1-\mu_2)}{\sqrt{\dfrac{S_1^2}{m}+\dfrac{S_2^2}{n}}}
$$

has approximately a $t$ distribution with df $\nu$ estimated by

$$
\nu = \frac{\left(\dfrac{s_1^2}{m}+\dfrac{s_2^2}{n}\right)^2}{\dfrac{(s_1^2/m)^2}{m-1}+\dfrac{(s_2^2/n)^2}{n-1}}
= \frac{\left[(se_1)^2+(se_2)^2\right]^2}{\dfrac{(se_1)^4}{m-1}+\dfrac{(se_2)^4}{n-1}},
$$

where $se_1=s_1/\sqrt{m}$ and $se_2=s_2/\sqrt{n}$. This is called **Welch’s df** (round $\nu$ down to an integer).

Then

$$
P\!\left(-t_{\alpha/2,\nu} \lt  T = \frac{\bar{X}-\bar{Y}-(\mu_1-\mu_2)}{\sqrt{S_1^2/m+S_2^2/n}} \lt  t_{\alpha/2,\nu}\right) \approx 1-\alpha .
$$

**The two-sample CI.** The two-sample CI for $\mu_1-\mu_2$ with confidence level approximately $100(1-\alpha)\%$ is

$$
\bar{x}-\bar{y} \pm t_{\alpha/2,\nu}\sqrt{\frac{s_1^2}{m}+\frac{s_2^2}{n}},
$$

where $\nu$ is Welch’s df.

**Test procedure (two-sample $t$ test).** Null hypothesis $H_0:\mu_1-\mu_2=\Delta_0$. Test statistic value:

$$
t = \frac{\bar{x}-\bar{y}-\Delta_0}{\sqrt{\dfrac{s_1^2}{m}+\dfrac{s_2^2}{n}}}.
$$

| Alternative hypothesis       | Rejection region                                                        |
|:-----------------------------|:------------------------------------------------------------------------|
| $H_a:\mu_1-\mu_2\gt \Delta_0$   | $t\ge t_{\alpha,\nu}$ (upper-tailed)                                    |
| $H_a:\mu_1-\mu_2\lt \Delta_0$   | $t\le -t_{\alpha,\nu}$ (lower-tailed)                                   |
| $H_a:\mu_1-\mu_2\ne\Delta_0$ | either $t\ge t_{\alpha/2,\nu}$ or $t\le -t_{\alpha/2,\nu}$ (two-tailed) |

The $P$-value can be written similarly to the one-sample $t$ test.

**Large-sample procedure.**

- If $m$ and $n$ are large ($m\gt 40$ and $n\gt 40$), by the CLT we get an approximately normal distribution.

- If $m$ and $n$ are large, then $\nu\ge\min(m-1,n-1)$ is also large. Then we can replace the $t$ critical value with the $z$ critical value to get an approximation to the $t$ procedure.

**Pooled $t$ procedures.** Assume $\sigma_1=\sigma_2=\sigma$ (even if it is unknown). Then we have the random variable

$$
Z = \frac{\bar{X}-\bar{Y}-(\mu_1-\mu_2)}{\sqrt{\dfrac{\sigma^2}{m}+\dfrac{\sigma^2}{n}}} = \frac{\bar{X}-\bar{Y}-(\mu_1-\mu_2)}{\sigma\sqrt{\dfrac{1}{m}+\dfrac{1}{n}}}.
$$

To get rid of $\sigma$ we could:

1.  replace $\sigma$ with $S_1$;

2.  replace $\sigma$ with $S_2$;

3.  use the ordinary average $(S_1+S_2)/2$ (very naive);

4.  best choice: use the **pooled** (i.e. combined) estimator of $\sigma^2$,

    $$
    S_p^2 = \frac{m-1}{m+n-2}S_1^2 + \frac{n-1}{m+n-2}S_2^2,
    \qquad\text{for which}\qquad
    \frac{(m+n-2)S_p^2}{\sigma^2}\sim\chi^2_{m+n-2}\ \ (\text{exercise}).
    $$

The random variable is

$$
T = \frac{\bar{X}-\bar{Y}-(\mu_1-\mu_2)}{S_p\sqrt{\dfrac{1}{m}+\dfrac{1}{n}}} \sim t_{m+n-2}.
$$

**Pooled $t$ CI.** The pooled $t$ CI for $\mu_1-\mu_2$ with confidence level $100(1-\alpha)\%$ is

$$
\bar{x}-\bar{y}\pm t_{\alpha/2,m+n-2}\cdot s_p\sqrt{\frac{1}{m}+\frac{1}{n}} .
$$

**Pooled $t$ test.** Null hypothesis $H_0:\mu_1-\mu_2=\Delta_0$. Test statistic value:

$$
t = \frac{\bar{x}-\bar{y}-\Delta_0}{s_p\sqrt{\dfrac{1}{m}+\dfrac{1}{n}}},
$$

with rejection regions as before but with $m+n-2$ df.

- When $\sigma_1=\sigma_2$, the pooled $t$ procedure is stronger than the two-sample (Welch’s) $t$ procedure.

- When $\sigma_1\ne\sigma_2$, using the pooled $t$ gives incorrect results.

- When $\sigma_1\ne\sigma_2$ (or we are not sure), use the two-sample $t$ procedure unless there is strong evidence for using the pooled $t$ procedure.

**Power and type II errors, $\beta(\Delta')$.**

- It is difficult to find the exact power and type II error probability.

- When $m$ and $n$ are large, we can use the formulas obtained for the two-sample $z$ procedure; in all those formulas, replace the $\sigma$’s with their estimates.

## 10.3 Analysis of Paired Data

**Assumptions.** Suppose the data consist of $n$ independently selected pairs $(X_1,Y_1)$, $(X_2,Y_2)$, …, $(X_n,Y_n)$ with $E(X_i)=\mu_1$ and $E(Y_i)=\mu_2$. Let $D_i=X_i-Y_i$. The $D_i$’s are assumed to be normally distributed with mean $\mu_D$ and standard deviation $\sigma_D$.

- The $D_i$’s are independent of each other.

- $X$ and $Y$ within a pair are *not* independent.

Consider $D=X-Y$. Then $E(D)=E(X-Y)=E(X)-E(Y)=\mu_1-\mu_2$, so $\mu_D=\mu_1-\mu_2$, and a CI for $\mu_D$ is equivalent to a CI for $\mu_1-\mu_2$. Let $\bar{D}$ be the sample mean and $S_D$ the sample standard deviation of the differences. Then

$$
T = \frac{\bar{D}-\mu_D}{S_D/\sqrt{n}} \sim t_{n-1}.
$$

**Paired CI.** The paired $t$ CI for $\mu_D=\mu_1-\mu_2$ with confidence level $100(1-\alpha)\%$ has endpoints

$$
\bar{d}\pm t_{\alpha/2,n-1}\cdot\frac{s_D}{\sqrt{n}},
$$

where $\bar{d}$ and $s_D$ are the observed values. The one-sided confidence bounds are $\bar{d}-t_{\alpha,n-1}\,s_D/\sqrt{n}$ (lower) and $\bar{d}+t_{\alpha,n-1}\,s_D/\sqrt{n}$ (upper).

**The paired $t$ test.** Null hypothesis $H_0:\mu_D=\Delta_0$. Test statistic value: $t=\dfrac{\bar{d}-\Delta_0}{s_D/\sqrt{n}}$.

| Alternative hypothesis | Rejection region                                                        |
|:-----------------------|:------------------------------------------------------------------------|
| $H_a:\mu_D\gt \Delta_0$   | $t\ge t_{\alpha,n-1}$ (upper-tailed)                                    |
| $H_a:\mu_D\lt \Delta_0$   | $t\le -t_{\alpha,n-1}$ (lower-tailed)                                   |
| $H_a:\mu_D\ne\Delta_0$ | either $t\le -t_{\alpha/2,n-1}$ or $t\ge t_{\alpha/2,n-1}$ (two-tailed) |

The $P$-value is computed as in the earlier $t$ test procedures.

**Comparison of the two statistics.**

$$
\text{Paired } t:\ \ t=\frac{\bar{d}-\Delta_0}{s_D/\sqrt{n}}
\qquad\qquad
\text{Two-sample } t:\ \ t=\frac{\bar{x}-\bar{y}-\Delta_0}{\sqrt{\dfrac{s_1^2}{m}+\dfrac{s_2^2}{n}}}
$$

## 10.4 Inferences About Two Population Proportions

Let $X\sim\text{Bin}(m,p_1)$ and $Y\sim\text{Bin}(n,p_2)$ be independent. Consider the parameter $p_1-p_2$. The natural estimator of $p_1-p_2$ is $\hat{p}_1-\hat{p}_2$, where $\hat p_1=X/m$, $\hat p_2=Y/n$. For large samples $\hat{p}_1$ and $\hat{p}_2$ have approximately normal distributions, so $\hat{p}_1-\hat{p}_2$ also has an approximately normal distribution.

$$
E(\hat{p}_1-\hat{p}_2) = E(\hat{p}_1)-E(\hat{p}_2) = E\!\left(\frac{X}{m}\right)-E\!\left(\frac{Y}{n}\right) = \frac{E(X)}{m}-\frac{E(Y)}{n} = \frac{mp_1}{m}-\frac{np_2}{n} = p_1-p_2,
$$

$$
V(\hat{p}_1-\hat{p}_2) = V(\hat{p}_1)+V(\hat{p}_2) = \frac{V(X)}{m^2}+\frac{V(Y)}{n^2} = \frac{mp_1q_1}{m^2}+\frac{np_2q_2}{n^2} = \frac{p_1q_1}{m}+\frac{p_2q_2}{n},
$$

where $q_1=1-p_1$, $q_2=1-p_2$. So $\sigma_{\hat{p}_1-\hat{p}_2}=\sqrt{\dfrac{p_1q_1}{m}+\dfrac{p_2q_2}{n}}$, and we have the random variable

$$
Z = \frac{\hat{p}_1-\hat{p}_2-(p_1-p_2)}{\sqrt{\dfrac{p_1q_1}{m}+\dfrac{p_2q_2}{n}}}. \qquad(1)
$$

**A large-sample test procedure.** Let $H_0:p_1-p_2=0$ be the null hypothesis. Let $p$ denote the common proportion, i.e. $p_1=p_2=p$, and $q=1-p$. Then the random variable (1) becomes

$$
Z = \frac{\hat{p}_1-\hat{p}_2-0}{\sqrt{pq\left(\dfrac{1}{m}+\dfrac{1}{n}\right)}}. \qquad(2)
$$

To get rid of $p$ in the denominator, we can use the sample proportion as an estimate. The estimator of $p$ is

$$
\hat{p} = \frac{X+Y}{m+n} = \frac{X}{m+n}+\frac{Y}{m+n} = \frac{m}{m+n}\hat{p}_1+\frac{n}{m+n}\hat{p}_2, \qquad(3)
$$

which is called the *weighted average* of the estimators. Substituting (3) into (2), we obtain

$$
Z^{\ast} = \frac{\hat{p}_1-\hat{p}_2}{\sqrt{\hat{p}\hat{q}\left(\dfrac{1}{m}+\dfrac{1}{n}\right)}},\qquad \hat{q}=1-\hat{p}.
$$

**Two-proportion $z$ test.** Null hypothesis $H_0:p_1-p_2=0$. Test statistic value:

$$
z = \frac{\hat{p}_1-\hat{p}_2}{\sqrt{\hat{p}\hat{q}\left(\dfrac{1}{m}+\dfrac{1}{n}\right)}}.
$$

| Alternative hypothesis | Rejection region                                   |
|:-----------------------|:---------------------------------------------------|
| $H_a:p_1-p_2\gt 0$        | $z\ge z_\alpha$                                    |
| $H_a:p_1-p_2\lt 0$        | $z\le -z_\alpha$                                   |
| $H_a:p_1-p_2\ne0$      | either $z\le -z_{\alpha/2}$ or $z\ge z_{\alpha/2}$ |

The $P$-value is computed as in the previous $z$ procedures. These procedures are valid when $m\hat{p}_1\ge10$, $m\hat{q}_1\ge10$, $n\hat{p}_2\ge10$ and $n\hat{q}_2\ge10$.

**Power, $\beta$, and sample size determination.** Because of the definition of $\beta$, we cannot use the random variables (2) and (3), as they assume $H_0$ is true; in the definition of $\beta$, $H_0$ is false. $\beta$ is a function of both $p_1$ and $p_2$, i.e. $\beta(p_1,p_2)$ (textbook p. 605). For example, for $H_a:p_1-p_2\gt 0$,

$$
\beta(p_1,p_2) = \Phi\!\left(\frac{z_\alpha\sqrt{\bar{p}\bar{q}\left(\frac{1}{m}+\frac{1}{n}\right)}-(p_1-p_2)}{\sigma}\right),
\quad
\bar p=\frac{mp_1+np_2}{m+n},\ \bar q = 1-\bar p,\ \sigma=\sqrt{\frac{p_1q_1}{m}+\frac{p_2q_2}{n}} .
$$

Setting $\beta(p_1,p_2)=\beta$ and solving gives the sample size.

**A large-sample two-proportion CI.** The CI for $p_1-p_2$, using the random variable (1), is

$$
\hat{p}_1-\hat{p}_2 \pm z_{\alpha/2}\sqrt{\frac{\hat{p}_1\hat{q}_1}{m}+\frac{\hat{p}_2\hat{q}_2}{n}} .
$$

This CI sometimes deviates from its nominal level. To overcome this, add one success and one failure to each sample; then the confidence level is approximately correct. This adjusted CI can even be used for small samples.

**Small-sample procedures.**

- Fisher’s exact test.

- Add one success and one failure (adjusted CI).

- Use a software package to simulate the binomial distributions.

## 10.5 Inferences About Two Population Variances

**Theorem.** Let $X_1,\dots,X_m$ and $Y_1,\dots,Y_n$ be independent random samples from normal distributions with standard deviations $\sigma_1$ and $\sigma_2$. Let $S_1$ and $S_2$ be the sample standard deviations. Then the random variable

$$
F = \frac{S_1^2/\sigma_1^2}{S_2^2/\sigma_2^2}
$$

has an $F$ distribution with $\nu_1=m-1$ and $\nu_2=n-1$ degrees of freedom.

**The $F$ test for equality of variances.** Null hypothesis $H_0:\sigma_1=\sigma_2$ ($\sigma_1^2=\sigma_2^2$). Test statistic value: $f=s_1^2/s_2^2$.

| Alternative hypothesis    | Rejection region                                                            |
|:--------------------------|:----------------------------------------------------------------------------|
| $H_a:\sigma_1\gt \sigma_2$   | $f\ge F_{\alpha,\nu_1,\nu_2}$                                               |
| $H_a:\sigma_1\lt \sigma_2$   | $f\le F_{1-\alpha,\nu_1,\nu_2}$                                             |
| $H_a:\sigma_1\ne\sigma_2$ | either $f\le F_{1-\alpha/2,\nu_1,\nu_2}$ or $f\ge F_{\alpha/2,\nu_1,\nu_2}$ |

**$P$-value calculation:** use software.

**Exercise.** Derive the CI for $\sigma_1/\sigma_2$.

*Answer.* From $P\big(F_{1-\alpha/2,\nu_1,\nu_2} \lt  \frac{S_1^2/S_2^2}{\sigma_1^2/\sigma_2^2} \lt  F_{\alpha/2,\nu_1,\nu_2}\big)=1-\alpha$, a $100(1-\alpha)\%$ CI for $\sigma_1^2/\sigma_2^2$ is

$$
\left(\frac{s_1^2}{s_2^2}\cdot\frac{1}{F_{\alpha/2,\nu_1,\nu_2}},\ \frac{s_1^2}{s_2^2}\cdot\frac{1}{F_{1-\alpha/2,\nu_1,\nu_2}}\right),
$$

and the CI for $\sigma_1/\sigma_2$ is obtained by taking square roots of the endpoints.
