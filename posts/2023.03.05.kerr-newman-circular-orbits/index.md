@def title = "Circular orbits of charged particles"
@def date = Date(2023, 03, 05)
@def tags = ["2023", "julia", "relativity"]

# {{fill title}}

In the last couple of days I've been trying to add the [Kerr-Newman](https://en.wikipedia.org/wiki/Kerr%E2%80%93Newman_metric) metric to [Gradus.jl](https://github.com/astro-group-bristol/Gradus.jl). The Kerr-Newman metric describes a black hole spacetime with mass $M$, spin $a$, and charge $Q$. The metric may be expressed in Boyer-Lindquist coordinates as

$$
\d s^2 = 
    -\frac{\Delta - a^2 + \sin^2 \theta}{\Sigma}  \d t^2 
    + \frac{\Sigma}{\Delta} \d r^2 
    + \Sigma \d \theta^2 
    + \frac{\sin^2 \theta}{\Sigma} \left[\left(r^2 + a^2\right)^2 - a^2 \Delta \sin^2 \theta \right] \d \phi^2 
    - 2 \frac{a \sin^2 \theta}{\Sigma} \left(r^2 + a^2 - \Delta \right)  \d t \d \phi,
$$

where 

$$
\Sigma = r^2 + a^2 \cos^2 \theta, \quad \text{and} \quad \Delta = r^2 - 2 M r + a^2 + Q^2.
$$


Implementing the metric itself is straight forward, allowing for a new class of geodesic to be traced, namely those with non-zero charge per unit mass $q$. These are the particles which interact electromagnetically with the central singularity. We consider only time-like geodesics with charge.

These have been discussed in a number of papers, for example [Schroven, Hackmann and Lämmerzahl (2017)](#fndef:kris)[^kris], which studies in detail the innermost stable circular orbit (ISCO) for different charge configuration, and the effect this has on accreting matter of negative (electrons) and positive charge (protons/ions). There is also the earlier paper by [Hackmann and Xu (2013)](#fndef:hackmann)[^hackmann], which explores more generally the orbits of charged particles in the Kerr-Newman spacetime and classifies them to describe their behaviour. These works, and references therein, use the Hamilton-Jacobi formalism with the Carter constant to calculate an effective potential on the four-velocity components, and use this to infer stable (circular) orbits and the ISCO. I would instead like to study this from a step back, directly from the 2nd order geodesic equation, using modern numerical methods. This means adding an Ansatz for the Lorentz force

$$
\label{eq:geod}
\frac{\d^2}{\d \lambda^2} x^\mu = - \tensor{\Gamma}{\mu}{\alpha\beta} \dot{x}^\alpha \dot{x}^\beta + q \tensor{F}{\mu}{\alpha} \dot{x}^\alpha,
$$

where $\tensor{F}{\mu}{\alpha}$ is the Faraday tensor, and the dot denotes differentiation with respect to affine time $\lambda$. The second-order ordinary differential equation may be solved with some initial position $x^\mu$ and 3-velocity $\vec{\dot{x}}$, as $\dot{x}^t$ is constrained via

$$
\label{eq:constraint}
g_{\alpha \beta} \dot{x}^\alpha \dot{x}^\beta = -\mu^2
$$

where $\mu$ is the invariant mass.

In this post, I would like to study the class of circular charged orbits in the equatorial plane, beginning from eq. \eqref{eq:geod}, to determine semi-analytic solutions for the angular momentum and energy of these orbits.

## Faraday tensor

Starting from the four-potential[^kris]

$$
A_\mu = \frac{Q r}{\Sigma} \left(1 \d t - a \sin^2 \theta \d \phi \right),
$$

the Faraday tensor is given by

$$
F_{\mu\nu} = \partial_\mu A_\nu - \partial_\nu A_\mu.
$$

This means there are only 4 unique non-zero components of the anti-symmetric tensor, namely $F_{rt}, F_{\theta t}, F_{r\phi}, F_{r\theta}$, with the other 4 being their anti-symmetric counterparts.

It is worth noting that with this potential, one must only consider the product of the particle and black hole charge $qQ$ in the equations of motion, since there are no terms that depend solely on $q$ or $Q$. There are therefore 2 classes of solution we can expect depending on the sign of $qQ$. Similarly in the metric, only terms with $Q^2$ appear and consequently the parity is also lost. We therefore restrict ourselves to $Q \geq 0$ and allowing $q$ to take both negative and positive values without loss of generality.

## Keplerian angular velocity


Following the method of [my previous blog post](../2022.05.09.plunging-photons), we will rework eq. \eqref{eq:geod} in the equatorial plane, subject to $\dot{x}^r = \dot{x}^\theta = 0$, and calculate the Keplerian angular velocity, defined

$$
\Omega_\phi := \frac{\dot{x}^\phi}{\dot{x}^t}.
$$

The energy and angular momentum already have known solutions in terms of $\Omega_\phi$ and metric terms. After expanding the Christoffel symbols and working through a little algebra (as in the previous post), one obtains

$$
0 = \frac{1}{2}\left( \frac{\partial}{\partial r} g_{tt} (\dot{x}^t)^2
    + 2 \frac{\partial}{\partial r} g_{t\phi} \, \dot{x}^t \dot{x}^\phi
    +\frac{\partial}{\partial r} g_{\phi \phi} (\dot{x}^\phi)^2 \right)
    + q g_{rr} \left(\tensor{F}{r}{t} \dot{x}^t + \tensor{F}{r}{\phi} \dot{x}^\phi \right),
$$

and dividing by $(\dot{x}^t)^2$, one obtains

$$
\label{eq:expanded_geod}
0 = \frac{1}{2}\left( \frac{\partial}{\partial r} g_{tt}
    + 2 \frac{\partial}{\partial r} g_{t\phi} \Omega_\phi
    +\frac{\partial}{\partial r} g_{\phi \phi} \Omega_\phi^2 \right)
    + q \frac{g_{rr}}{\dot{x}^t} \left(\tensor{F}{r}{t} + \tensor{F}{r}{\phi} \Omega_\phi \right),
$$

which bears a factor $1/\dot{x}^t$ on the charge terms. The first term in this equation is a simple quadratic, but with these charge terms added, we require some way of eliminating or determining $\dot{x}^t$ in order to solve for $\Omega_\phi$.

The only thing I can think to do here is to use eq. \eqref{eq:constraint}, and write

$$
\label{eq:ut}
\frac{1}{\dot{x}^t} = \frac{1}{\mu} \sqrt{- \left(g_{tt} + 2g_{t\phi} \Omega_\phi + g_{\phi \phi} \Omega_t^2 \right)}.
$$

One could then balance eq. \eqref{eq:expanded_geod} to have one term either side of the equals sign, substitute the $1/\dot{x}^t$, and square to obtain a quartic equation. Doing so also loses information about the sign of $q$, but that is something we could add in later. Although this is in theory then analytically solvable, in reality the expression is horrendous, and I cannot find a nice simple way of reducing the quartic into e.g. the product of two quadratics, or a depressing it to simplify the expression. There are no doubt other ways to tackle such an equation, but I am not versed in polynomial analysis and would have a hard time approaching the problem.

For now, I'd instead prefer a semi-analytic approach, and use a root finder to solve for the solutions. Since the overall form of the equation is leading order quadratic, and we can _a priori_ select $\dot{x}^t$ that goes forwards in time, we may therefore anticipate classes of solutions for $\Omega_\phi$, corresponding to product set of prograde and retrograde orbits with positively and negatively charged particles.

Implementing a short Julia function:

```julia
using Gradus

function functor_Ω(
    m::KerrNewmanMetric,
    rθ;
    q = 0.0,
    μ = 1.0
)
    g, jacs = Gradus.metric_jacobian(m, rθ)
    # only want the derivatives w.r.t. r
    ∂rg = jacs[:, 1]

    x = SVector(0, rθ[1], rθ[2], 0)
    F = Gradus.maxwell_tensor(m, x)

    function f(ω)
        Δ = (ω^2 * ∂rg[4] + 2 * ω * ∂rg[5] + ∂rg[1])
        arg = -(ω^2 * g[4] + 2 * ω * g[5] + g[1]) / μ^2
        # check to ensure we can preserve forwards-in-time
        ut = if arg >= 0
            sqrt(arg)
        else
            # return something distinct
            return NaN
        end
        # combine with charge terms
        0.5 * Δ + (F[2, 4] * ω + F[2, 1]) * g[2] * q * ut
    end
end
```

This function limits the domain of $\Omega_\phi$ to solutions with $\frac{1}{u^t} \geq 0$. We then set up some initial configuration:

```julia
m = KerrNewmanMetric(M = 1.0, a = 0.0, Q = 0.9)
rθ = SVector(3.0, π/2)
```

I'd like to note that this picture is unphysical. If a black hole _were to have charge_, it is likely that the charge is vanishingly small in the metric[^kris], and that the charge product $qQ$ is instead non-zero. But there is nothing in the mathematics that prevents us from exploring this class of spacetimes, so we will press on.

Plotting the value of the constraint equation for different $\Omega_\phi$ for two different spins:

![](./omega-plot.svg)

Discontinuities are due to limited resolution ($N_{\Omega_\phi}=300$). The two solutions, prograde and retrograde, are clearly defined. When $|q| > 0.0$, the possible range of $\Omega_\phi$ is limited by the positivity constraint in eq. \eqref{eq:ut}. 
- For the Schwarzschild case $a=0.0$, the function is symmetric about the central axis. The roots differ only by parity, $\Omega_{\phi,1} = - \Omega_{\phi, 2}$, for both the charged and uncharged orbits. We identify the effect of charge is to either increase ($q=-1$) or decrease ($q=1$) the angular velocity relative to the uncharged case.
- For the spinning case $a = 0.4$, the function is skewed and loses its symmetry, but the effects of charge remain the same.

We can clearly identify the values of $\Omega_\phi$ corresponding to each of the $q = \{-1, 0, 1\}$ cases. Ambiguity is introduced if we move closer:

![](./omega-plot-2.svg)

Examining the two cases again:
- There are still four solutions to the charged particle orbits in the Schwarzschild case, but only if one considers the $q = \pm 1$ curves to be continuous. Else we must conclude that there are four possible orbits for $q=1$ and none for $q = -1$ at the same radius, which we can rule out, as circular orbits _must_ have a unique $\Omega_\phi$, modulus the direction (prograde or retrograde).
- For the spinning case, we identify only two real roots, corresponding to prograde orbits, both on the $q = 1$ branch. Again, we must consider the $q = \pm 1$ curves to be continuous to rationalize the roots. At this distance, retrograde circular orbits for the spinning case are impossible.

I am unsure what the interpretation here is: one could interpret that the direction of time has changed, as in eq. \eqref{eq:ut}, there is no difference between $q \rightarrow -q$ and $u^t \rightarrow -u^t$. Thus the transformation $(q, u^t) \rightarrow (q, -u^t)$ is equivalent to the $(-q, u^t)$ case, which puts negative charge orbits on the curve corresponding to positive charge orbits when the direction of time is changed.

We can analytically continue the curves to represent this with the modification to our function:
```julia
ut = sign(arg) * sqrt(abs(arg))
```

![](./omega-plot-3.svg)

I have coloured the line corresponding to $q = 1$ for $a = 0.4$ orange to emphasize how the continuation is realized. I have chosen to adopt this in the Gradus.jl implementation. It also helps with the root finder for $\Omega_\phi$ to have the function defined all along the domain (~;

## Energy and ISCO

Here we will use the [method outlined previously](posts/2022.05.09.plunging-photons/#determining_the_isco) for determining the ISCO, using the same minima finding method in the energy expression. The energy and angular momentum are from the canonical four-momenta, which is modified under the influence of the electromagnetic potential, as described in [Tursunov, Zdeněk, and Kološ (2016)](#fndef:tursovon)[^tursunov] for an external magnetic field, and [Carter (1968)](#fndef:carter)[^carter] for the general class of Kerr metrics.

Quickly re-deriving the canonical momenta from the Lagrangian

$$
\mathcal{L} = \frac{1}{2} g_{\mu\nu} \dot{x}^\mu \dot{x}^\nu - q A_\sigma \dot{x}^\sigma,
$$

with the momenta defined

$$
\pi_\nu := \frac{\partial \mathcal{L}}{\partial \dot{x}^\nu},
$$

and consequently

$$
\label{eq:canonical-momentum}
\pi_\nu = g_{\nu\mu}\dot{x}^\mu - q A_\nu = \dot{x}_\nu - q A_\nu.
$$

There is here an ambiguity in the sign of $A_\nu$ as it enters in the above equations. I am using the same convention as Hackmann et al. here, which is in contrast to e.g. Carter. I must admit I do not quite understand the origin of this discrepancy, but the differing results are dramatic.

We identify as usual $E = -\pi_t$ and $L_z = \pi_\phi$ as constants of motion, and 

$$
\therefore \quad E = -(\dot{x}_t - q A_t), \quad \text{and} \quad L_z = (\dot{x}_\phi - q A_\phi).
$$


We find $\dot{x}_t$ and $\dot{x}_\phi$ using $\Omega_\phi$, as in the previous blog post, and can then plot the energy as a function of radius:

![](./energy-charged.svg)

For the case of $q=1$, we have $E \approx 1$ for the majority of orbits. The _marginally bound_ orbit is defined as the circular orbit with $E = 1$, so although these orbits are stable, they are not very tightly bound relative to the other setups. This has implications for a maximal charge product $qQ$ for which circular orbits may be found. Numerical limitations furthermore prevent fully describing the energetic system at small radii, as can be seen by blatant discontinuities in the curves, where the root finder has struggled to find solutions for $\Omega_\phi$.

We can find the ISCO by determining the minima of the energy of the circular orbits, i.e.

$$
\left.\frac{\partial E}{\partial r}\right\rvert_{r=r_\text{ISCO}} = 0.
$$

Plotting the ISCO as a function of charge product $qQ$:

![](./charged-iscos.svg)

At $qQ = 1$, both curves are divergent and $r_\text{ISCO} \rightarrow \infty$. For negative charge products the ISCO increases -- this is equivalent to an additional attractive force, or an effective increase in gravitational potential, pulling on the particle. Orbits that were once stable are no longer, and hence $r_\text{ISCO}$ increases.

For positive charge products, the ISCO initially decreases, reaching a minima, before increasing rapidly. The initial decrease may be thought of as the repulsion holding up the orbit and weakening the gravitational potential, allowing previously unstable orbits to become stable. The subsequent increase is interpreted as the EM repulsion overpowering gravity, suppressing the gravitational term in eq. \eqref{eq:geod}.

## Orbit finding approach

Using the method I previously outlined in [another post](../2022.03.02.circular-orbits-kerr), Gradus.jl can automatically find stable circular orbits in the equatorial plane by minimizing some stability measure. It requires only the modified geodesic equation eq. \eqref{eq:geod}, and solves for some $v^\phi$ which minimizes the variations in $\Delta x^r$ over the course of the integration. It is therefore a good numerical test of the semi-analytic approach to ensure the model is at least self consistent.

Plotting $L_z$ against $E$ is a method for analyzing the stability and classifying different orbits. For example, the ISCO is a minimum of both $E$ and $L_z$, and therefore on the $L_z - E$ plane corresponds to a cusp in the curve representing circular orbits.

Plotting the analytic method (dashed) and the points determined by the circular orbit solver (crosses, connected by solid line):

![](./el-plot-1.svg)

Agreement is seen between the methods. This can also be used to more accurately determine the ISCO for the circular orbit solver, as sometimes the integrator is "too good" and classifies some unstable / marginally bound orbits as stable circular. These plots clearly identify the location of the ISCO energetically.

## Line profiles

The effects of GR on the redshift is calculated using the ratio of energies

$$
g := \frac{\left. k^\nu u_{\nu} \right\rvert_{\text{obs}}}{\left. k^\mu u_\mu \right\rvert_\text{disc}},
$$

where we assume $\left. k^\nu u_{\nu} \right\rvert_{\text{obs}} = -1$ for simplicity, as the observer is assumed to be far enough from the singularity to be in flat space. Our calculation therefore simplifies to

```julia
using Tullio

function redshift_charged(q)
    function _f(m, gp, t)
        x = gp.u2
        g = Gradus.metric(m, x)
        u = CircularOrbits.fourvelocity(m, x[2]; q = q)
        @tullio o := -g[i, j] * gp.v2[i] * u[j]
        inv(o)
    end
end
```

Here it is already implemented in a form that returns a function that we can use to create a `PointFunction` to use with the integrator. We can then calculate moderately high resolution line profiles with:

```julia
# set charge and build point function
q = 1.0
redshift_pf = PointFunction(redshift_charged(q))
# position and disc configuration
u = SVector(0.0, 1000.0, deg2rad(40), 0.0)
d = GeometricThinDisc(0.0, 1000.0, π/2)
# set up binning and calculate over fine grid
bins = collect(range(0.1, 1.3, 200))
g, f = @time lineprofile(
    m,
    u,
    d,
    bins = bins,
    algorithm = BinnedLineProfile(),
    plane = PolarPlane(GeometricGrid(); Nr = 700, Nθ = 1300, r_max = 50.0),
    verbose = true,
    redshift_pf = predshift_pf,
    minrₑ = Gradus.isco(m, q = q),
)
```

Plotting these:

![](./high-charge-line-profiles.svg)

The effect of the energetics of the charged orbits is clear in the line profiles. For the case of $q=1$, the circular orbit energy is pretty much constant for all stable orbits. This is particularly important around the ISCO radius, as that is where we can expect the maximal redshift. Consequently, the redshift of the line profile is much narrower. For the case of $q = -1$ we see the opposite effect, where the energy of the stable orbits varies to a greater extent, and therefore the line profile is significantly broadened. 

The effect is somewhat similar to changing the inclination of the observer:

![](./diff-incl-line-profiles.svg)

The above is for the Schwarzschild case only.

A "more realistic" scenario anticipates $Q$ to be very small, but could still allow $qQ \sim 1$. We can approximate this by setting e.g. $Q \approx 10^{-16}$ and $q \approx 10^{16}$, so that $Q$ does not affect the spacetime, and enters the model only through the Lorentz force in the geodesic equation.

![](./low-charge-line-profiles.svg)

The Schwarzschild and $qQ = 0$ case correspond exactly, but the effects of the charge are still very apparent.


## References

[^kris]: Schroven, Kris, Eva Hackmann, and Claus Lämmerzahl. ‘Relativistic Dust Accretion of Charged Particles in Kerr-Newman Spacetime’. Physical Review D 96, no. 6 (26 September 2017): 063015. https://doi.org/10.1103/PhysRevD.96.063015.
[^hackmann]: Hackmann, Eva, and Hongxiao Xu. ‘Charged Particle Motion in Kerr-Newmann Space-Times’. Physical Review D 87, no. 12 (24 June 2013): 124030. https://doi.org/10.1103/PhysRevD.87.124030.
[^tursunov]: Tursunov, Arman, Zdeněk Stuchlík, and Martin Kološ. ‘Circular Orbits and Related Quasiharmonic Oscillatory Motion of Charged Particles around Weakly Magnetized Rotating Black Holes’. Physical Review D 93, no. 8 (7 April 2016): 084012. https://doi.org/10.1103/PhysRevD.93.084012.
[^carter]: B Carter (1968). _Global structure of the Kerr family of gravitational fields_. PHYS REV **174**:5.
