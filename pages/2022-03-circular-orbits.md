@def title = "Circular orbits in static, axis-symmetric spacetimes"
@def description = "I'd like to be able to calculate programmatically key results related to circular orbits in static, axis-symmetric spacetimes. To this end, I want to see how feasible a computer aided approach is, using only a handful of constraints, in re-deriving known analytic solutions, and compare this to automated circular orbit discovery using our geodesic integrator."
@def image = "/assets/stablecirc.svg"
@def date = "02/03/2022"
@def tags = ["2022", "geodesics", "julia", "relativity"]

# Circular orbits in static, axis-symmetric spacetimes

<!--
Johannsen and Psaltis 2011
Bardeen et al 1972
Cunningham 1975
Gates et al 2021
Reynolds and Begelman 1997
Patra et al 2022
Ingram et al 2019
Igata et al 2020
-->

I'm going to try and re-derive 50 year old results in a day, programmatically, with modern computers, and see how close I can get. 

We'll look at results for orbits in the Kerr metric, using a mathematical approach that can be scripted for other static, axis-symmetric spacetimes, which is the underlying motivation; can we get a good enough of an approximation of tailored, analytic solutions to merit writing generic code for alternative metrics?

Of particular importance when studying metrics are the permitted circular orbits, especially for creating observational tests, as the orbits are used to model processes in e.g. accreting matter. Such trajectories also allow for the calculation of redshift values for discs when ray tracing, which in turn allows us to visualize and simulate other processes, such as flux. The little exploration here will hopefully lead to automated stable orbit discovery with [our integrator](https://github.com/astro-group-bristol), and (hopefully) allow us to eventually simulate arbitrary disc profiles for different metrics.

## Introduction

- Notation: *italic* is basis-vector, serif is Euler's exponential $\e$.

The main reference for this analysis is Bardeen, Press and Teukolsky (1972)[^1], which explores orbits in rotating black holes, and various emergent processes. The paper covers an awful lot, but we are only interested (for now) in Section II.

The standard, static, axis-symmetric, and asymptotically flat spacetime line element is

\begin{equation}
\label{eq:general-metric}
\d s^2  = 
    - \e^{2\nu} \d t^2 
    + \e^{2\psi} \left( \d\phi - \omega \d t \right)^2
    + \e^{2\mu_1} \d r^2
    + \e^{2\mu_2} \d \theta^2,
\end{equation}

where $\omega = -\frac{\tensor{g}{}{\phi t}}{\tensor{g}{}{\phi \phi}}$. 

We differentiate *stationary* and static *static*, where a stationary means the components of the metric *do not change with time* ($\frac{\partial}{\partial t} g_{\mu\nu} = 0$), and static adds symmetry under time parity, or simply that $t \rightarrow -t$ does not change the metric.

Then axis-symmetric implies invariance under translations in $\phi$ (things look the same when you rotate), and asymptotically flat means the metric reduces to the Minkowski metric, $\eta$, for $r \rightarrow \infty$. 

These metrics govern the structure of spacetime in the vicinity of some central singularity, and are all we need to study a black hole and its orbits.

This metric becomes Kerr when

\begin{equation}
\label{eq:metric-comp}
    \e^{2\nu} = \frac{\Sigma \Delta}{A}, 
        \quad
    \e^{2\psi} = \frac{A}{\Sigma} \sin^2 \theta,
        \quad
    \e^{2\mu_1} = \frac{\Sigma}{\Delta},
        \quad
    \e^{2\mu_2} = \Sigma
\end{equation}

and $\omega = \frac{2Mar}{A}$. The symbols have their usual meaning, e.g. from Cunningham (1975)[^2], namely

\begin{equation}
\Sigma = r^2 + a^2 \cos^2 \theta, \quad \Delta = r^2 - 2 M r + a^2, \quad A = (r^2+a^2)^2 - a^2 \Delta \sin^2 \theta.
\end{equation}

There's a future blog post in my take on the physical significance of all of these symbols, but that's not for today.

## Circular orbits in the Kerr metric

It is a known result that circular orbits exist in the equatorial plane of static axis-symmetric spacetimes, i.e. where $\theta = \frac{\pi}{2}$ and $\dot{\theta} = 0$, as shown by Carter (1968)[^3]. We can further constrain our mathematical picture by noting that circular orbits have no motion in the radial components,

\begin{equation}
    \label{eq:ur_dotur}
    \dot{u}^r = 0, \quad \text{and} \quad \frac{\partial}{\partial r}\dot{u}^r = 0.
\end{equation}

We can use the invariance of four velocities to express these constraints in terms of the metric,

\begin{equation}
    \label{eq:velconstraint}
    \tensor{g}{}{\sigma \nu} \tensor{\dot{u}}{\sigma}{} \tensor{\dot{u}}{\nu}{} = - \mu^2,
\end{equation}

where, for time-like geodesics, the effective mass parameter $\mu = 1$. In order to rearrange for the radial component, this equation may be expanded,

\begin{equation}

\tensor{g}{}{rr} \tensor{\dot{u}}{r}{} \tensor{\dot{u}}{r}{} 
+ \tensor{g}{}{tt} \tensor{\dot{u}}{t}{} \tensor{\dot{u}}{t}{} 
+ 2 \tensor{g}{}{t\phi} \tensor{\dot{u}}{t}{} \tensor{\dot{u}}{\phi}{} 
+ \tensor{g}{}{\phi\phi} \tensor{\dot{u}}{t}{} \tensor{\dot{u}}{t}{}
= - \mu^2,

\end{equation}

and, using the freedom to raise and lower indices under contraction, rewritten

\begin{equation}
    \left(\tensor{\dot{u}}{r}{} \right)^2 = 
        \frac{-1}{ \tensor{g}{}{rr} } \left(
            \tensor{g}{tt}{} \left(\tensor{\dot{u}}{}{t} \right)^2
            + 2 \tensor{g}{t\phi}{} \tensor{\dot{u}}{}{t} \tensor{\dot{u}}{}{\phi}
            + \tensor{g}{\phi\phi}{} \left(\tensor{\dot{u}}{}{\phi} \right)^2
            + \mu^2
        \right).
\end{equation}

In this form, the usual identification is permitted for energy $\tensor{\dot{u}}{}{t} = - E$ and angular momentum $\tensor{\dot{u}}{}{\phi} = L_z$ from the covariant momenta,

\begin{equation}
    \therefore \quad \left(\tensor{\dot{u}}{r}{} \right)^2 = 
        \frac{-1}{ \tensor{g}{}{rr} } \left(
            \tensor{g}{tt}{} E^2
            - 2 \tensor{g}{t\phi}{} E L_z
            + \tensor{g}{\phi\phi}{} L_z^2
            + \mu^2
        \right).
\end{equation}

The first constraint from eq. \eqref{eq:ur_dotur} may now be calculated from this expression, namely

\begin{equation}
    \label{eq:initial_con1}
    \tensor{\dot{u}}{r}{} = \frac{\pm 1}{ \e^{\mu_1} } \sqrt{
        \e^{-2\nu} \left(
            E + \omega L_z
        \right)^2
        - \e^{-2\Psi} L_z^2
        - \mu^2
    }
    = 0.
\end{equation}

If we rephrase the above constraint in terms of _effective potentials_, we find the expression cleans up to

\begin{equation}
    \tensor{\dot{u}}{r}{} = \pm \sqrt{ \mathcal{R} \, \e^{-2\mu_1} } = 0.
\end{equation}

Our constraints may be modified now to $\mathcal{R} = 0$ and $\frac{\partial}{ \partial r}\mathcal{R} = 0$, the latter of which is directly

\begin{equation}
\label{eq:sage-dont-like}
    \left( E + \omega L_z \right)^2 \frac{\partial}{\partial r} \e^{-2\nu}
    + 2 \left( E + \omega L_z \right) \e^{-2\nu} L_z \frac{\partial}{\partial r} \omega
    - L_z^2 \frac{\partial}{\partial r} \e^{-2 \Psi}
    = 0,
\end{equation}

calculated from eq. \eqref{eq:initial_con1}.

The term $\mathcal{R} \frac{\partial}{\partial r} \e^{-2\mu_1}$ is omitted, as it is manifestly zero, since $\mathcal{R} = 0$. Implicitly it is assumed that $E$ and $L_z$ are constant with respect to $r$ over an orbit. Indeed, this assumption is valid as $E$ and $L_z$ are constants of motion for any given geodesic, whereas the $E(r)$ and $L_z(r)$ which we later solve for are functions in the phase space. This distinction is subtle but important to note.  

Here we will specialize to the Kerr metric, with subsequent manipulations performed using SageMath, and the occasional small bit of hand manipulation. Substituting in the metric components from eq. \eqref{eq:metric-comp}, we obtain

\begin{equation}
    \label{eq:constraint-1}
    \mathcal{R} = 2 M r {\left(\mu^{2} r^{2} + {\left(E a - L_z\right)}^{2}\right)} + {\left({\left(E^{2} - \mu^{2}\right)} {\left(a^{2} + r^{2}\right)} - L_z^{2}\right)} \Sigma = 0,
\end{equation}

and

\begin{equation}
    \label{eq:constraint-2}
    \frac{\partial}{\partial r} \mathcal{R} = 6 \, M \mu^{2} r^{2} + 2 \, {\left(E a - L_z\right)}^{2} M + 2 \, {\left(E^{2} - \mu^{2}\right)} r \Sigma + {\left({\left(E^{2} - \mu^{2}\right)} {\left(a^{2} + r^{2}\right)} - L_z^{2}\right)} \frac{\partial}{\partial r}\Sigma = 0.
\end{equation}

It is worth noting that $\mu$ has reappeared in the derivative term. This is due to the complexity of eq. \eqref{eq:sage-dont-like} when taking the derivative in SageMath -- instead, SageMath _prefers_ to work with eq. \eqref{eq:constraint-1} to produce more visually friendly results. Both will eventually solve the system identically, so the choice is largely stylistic.

In the equatorial plane, $\Sigma = r^2$; therefore

\begin{equation}
    \label{eq:constraint_1}
    \mathcal{R} = 2 M r {\left(\mu^{2} r^{2} + {\left(E a - L_z\right)}^{2}\right)} + r^{2} {\left({\left(E^{2} - \mu^{2}\right)} {\left(a^{2} + r^{2}\right)} - L_z^{2}\right)} = 0
\end{equation}

and

\begin{equation}
    \label{eq:constraint_2_partial}
    \frac{\partial}{\partial r} \mathcal{R} =
    6 M \mu^{2} r^{2} 
    + 2 M {\left(E a - L_z\right)}^{2}
    + 2 r^3 \left(
        E^2 - \mu^2
    \right)
    + 2r \left(
        {\left(E^{2} - \mu^{2}\right)} {\left(a^{2} + r^{2}\right)} - L_z^{2}
    \right)
    = 0.
\end{equation}

Noticing a similar term in the two equations, we can rearrange eq. \eqref{eq:constraint_1} for

\begin{equation}
    2 r {\left({\left(E^{2} - \mu^{2}\right)} {\left(a^{2} + r^{2}\right)} - L_z^{2}\right)} = -4 M {\left(\mu^{2} r^{2} + {\left(E a - L_z\right)}^{2}\right)},
\end{equation}

and substitute into eq. \eqref{eq:constraint_2_partial} to obtain

\begin{equation}
    \frac{\partial}{\partial r} \mathcal{R} =   
    6 M \mu^{2} r^{2} 
    + 2 M {\left(E a - L_z\right)}^{2} 
    + 2 r^3 \left(
        E^2 - \mu^2
    \right)
    - 4 M {\left(\mu^{2} r^{2} + {\left(E a - L_z\right)}^{2}\right)}
    = 0,
\end{equation}

which may be simplified further:

\begin{equation}
\label{eq:constraint_2}
\therefore \quad 
\frac{\partial}{\partial r} \mathcal{R} =
    2 M \mu^{2} r^{2} 
    + 2 r^{3} {\left(E^{2} - \mu^{2}\right)} 
    - 2 M {\left(E a - L_z\right)}^{2}
    = 0.
\end{equation}

Our constrains in eq. \eqref{eq:ur_dotur} are now eq. \eqref{eq:constraint_1} and \eqref{eq:constraint_2} respectively, which may be solved simultaneously to find $E$ and $L_z$ satisfying circular orbits. In SageMath, this is fiddly, and sometimes SageMath will fail after a minute, or run for seemingly hours without terminating. But one soon develops a sense for what SageMath can and cannot do, and after this, it is generally quite straight forward to solve, even for other metrics. We find there are four solutions in total, representing the intersection of the positive and negative energy regimes with prograde and retrograde orbits,

\begin{equation}
\label{eq:Emu}
\frac{E}{\mu} =
    \pm \sqrt{
        \frac{
            r^5 - 7Mr^4 + 16M^2r^3-
            3(4M^3 + Ma^2)r^2 + 5M^2a^2r 
            \pm 2 M a \Delta \sqrt{M r}
        }
        {   
            r^5 - M r^2 \left(2 r^2 - Mr + 4 \Delta \right)
        }
    },
\end{equation}

where the first $\pm$ denotes the positive and negative energies, and the second $\pm$ the prograde and retrograde respectively. The angular momentum can be rearranged as

\begin{equation}
\label{eq:Lzmu}
    \frac{Lz}{\mu} =
    \frac{E}{\mu} a
    \pm \frac{r}{\sqrt{M}} \sqrt{
        1 - r + r \left(\frac{E}{\mu} \right)^2
    }.
\end{equation}

At the time of writing, this is the most simplified form I have been able to manipulate these equations into, which seems pretty okay since we haven't had to even touch any advanced techniques.

As an additional note, with SageMath, many of the preceding steps in preparing the constraints are not actually necessary, but would provide results that are unwieldy for such a blog post. I was also partially curious to see if simplifying the constraint equations dramatically would make any difference to the solutions (it didn't really).

### A comparison to Bardeen et al. (1972)

The Bardeen et al. paper quotes a more beautiful analytic solution

\begin{equation}
\label{eq:bardeen_e}
\frac{E}{\mu} = \frac{
        r^\frac{3}{2} - 2M r^\frac{1}{2} \pm a M^\frac{1}{2}
    }{
        r^{3/4}\left( r^\frac{3}{2} - 3Mr^\frac{1}{2} \pm 2 a M^\frac{1}{2} \right)^\frac{1}{2}
    },
\end{equation}

and

\begin{equation}
\frac{L_z}{\mu} = \frac{
        \pm M^\frac{1}{2} \left( 
            r^2 \mp 2 a M^\frac{1}{2} r^\frac{1}{2} + a^2    
        \right)
    }{
        r^{3/4}\left( r^\frac{3}{2} - 3Mr^\frac{1}{2} \pm 2 a M^\frac{1}{2} \right)^\frac{1}{2}
    },
\end{equation}

both of which are now ubiquitous in the literature.

They are also substantially different from eq. \eqref{eq:Emu} and eq. \eqref{eq:Lzmu}. We can try to find the degree of similarity by plotting the equations with $r > 0.1$ to avoid the central singularity:

![energy-diff](/assets/angmom-energy.svg)

The plots for energy (left panel) are in good agreement up until the energy becomes negative in the Bardeen case (though such energies are allowed in e.g. the [Penrose process](https://en.wikipedia.org/wiki/Penrose_process)). The SageMath derivation fails here, presumably due to an ambiguity when squaring roots, which loses sign information. Indeed, we do have a negative energy solution from SageMath, but we would require _a priori_ knowledge of which solution to apply in which domain. This wouldn't be ideal, and without a more refined solution, the SageMath approach would not even suggest negative energies are possible.

The consequence of the energy sign is also pronounced in the angular momentum (right panel). For the $a=0.998$ case, the plots are still in good agreement, but the divergence for $a=0.0$ is quite apparent, and worse still for $a=-0.988$ (not depicted, as the diagram became difficult to annotate). Even when using the Bardeen et al. energy solution in the SageMath expression for $L_z$, some of the features of the curves do not disappear, such as the sharp cusp at $r\sim2$, which, again, is problematic.

Fortunately, however, the differences in the two solutions is for our purposes pedagogical, since the divergence lies within the _photon orbit_, at which radii stable circular massive orbits are long since prohibited.


## Special orbit boundaries

There are many special boundaries in the Kerr metric. For purely circular orbits, we consider three of these boundaries and the respective possible orbits within them, as they are arguably most relevant to our study -- these are the _innermost stable circular orbits_ (ISCO), the _marginally bound orbits_, and the _photon orbit_. 

There is of course also the _event horizon_, but this is not an orbital boundary. We will draw reference to this boundary, but only with the need of understanding that we cease to consider geodesics after they would pass this horizon.

### ISCO

The _innermost stable circular orbit_, as the name suggests, are the set of orbits which are both circular and energetically parabolic (stable). Stability requires the energy be a minimum; this may be rephrased as

\begin{equation}
    \frac{\partial^2}{\partial r^2} \mathcal{R} \geq 0,
\end{equation}

which solves for an inequality $r > r_\text{ISCO}$, which is in the Bardeen paper. Physically, these are the orbits where perturbations to the potential will not cause the orbit to diverge; it is the energetically _preferred state_, if you will. 

### Marginally bound orbits

In direct contrast to the ISCO, we consider unstable circular orbits at $r < r_\text{ISCO}$. All such _unbound_ circular orbits have hyperbolic energies (will diverge from circular orbits when perturbed), with a specific radius hosting the _marginally bound_ orbits, defined as $r$ where $\frac{E}{\mu} = 1$. Any small perturbation of such a marginally bound orbit outwards will escape the black hole, but with asymptotically zero velocity.

The _marginally bound_ radius is may be solved for the Kerr metric as

\begin{equation}
    r_{\text{mb}} = 2M \mp a + 2 \sqrt{M^2 \mp aM },
\end{equation}

corresponding to first (infalling) radius at which the condition $\frac{E}{\mu} = 1$ is satisfied.

### Photon orbit

Photons are null geodesics with $\mu=0$, and, due to this, cannot form circular orbits anywhere except at a one particular radius, known as the _photon radius_, denoted $r_\text{ph}$. This radius is where $\frac{E}{\mu} \rightarrow \infty$, or, equivalently, where

\begin{equation}
r^2 - 3Mr \pm 2 a \sqrt{M r} = 0, 
\end{equation}

i.e. the denominator from eq. \eqref{eq:bardeen_e} vanishes. This equation is solved by Bardeen, with the solution

\begin{equation}
r_\text{ph} = 2M \left[ 
    1 + \cos\left(
        \frac{2}{3} \arccos \left( \frac{ \mp a}{M} \right)
    \right)
\right].
\end{equation}

The photon orbit represents, if you will, the innermost possible (unstable) circular orbit, which requires velocity equal to the speed of light. Orbits between with $r_\text{ph} < r < r_\text{ISCO}$ can only form unstable circular orbits, which may quickly and catastrophically collapse.

Analogous to the marginally bound orbits, the unstable orbits at the photon orbit, if perturbed outwards, will arrive at infinity still traveling at the speed of light.

## Tracing circular orbits

With solutions for the energy and angular momentum, and a brief overview of relevant orbital boundaries, we can begin to simulate some of these orbits.

Tracing stable circular orbits when $E$ and $L_z$ is known is relatively simple. We only need to calculate $\dot{u}^\phi$, since our integrator automatically constrains $\dot{u}^t$ via eq. \eqref{eq:velconstraint}. We find

\begin{equation}
    \dot{u}^\phi = g^{\phi\phi} L_z - g^{tt}E,
\end{equation}

using our eq. \eqref{eq:Emu} and eq. \eqref{eq:Lzmu} with $\mu=1.0$ for time-like geodesics. Setting up the integrator is then quite straight forward:

```julia
using GeodesicTracer
using ComputedGeodesicEquations
using GeodesicBase

# define `energy` and `angmom` functions
Δ(M, a, r) = r^2 + a^2 - 2 * M * r
function energy(M, a, r)
    var = 2 * M * a * Δ(M, a, r) * sqrt(M * r)
    denom = r^5 - M * r^2 * (4 * Δ(M, a, r) + 2r^2 - M * r)
    en = (
        5 * M^2 * a^2 * r + 16 * M^2 * r^3 - 7 * M * r^4 + r^5 -
        3 * (4 * M^3 + M * a^2) * r^2 - var
    )
    sqrt(abs(en / denom))
end

angmom(E, M, a, r) = E * a + (r / √M) * sqrt(abs(1 - r + r * E^2))
angmom(M, a, r) = angmom(energy(M, a, r), M, a, r)


# function for calculating the initial radial velocity
function phi_vel(m, u)
    metric = GeodesicBase.metric(m, u)
    inv_metric = inv(metric)

    E = my_energy(m.M, m.a, u[2])
    Lz = my_angmom(E, m.M, m.a, u[2])

    inv_metric[4, 4] * Lz - inv_metric[4, 1] * E
end
```

This code snippet defines utility functions for calculating the initial radial velocity component for a given metric and position. The actual tracing code is then a single function call:

```julia
using StaticArrays

m = BoyerLindquist(M=1.0, a=0.0)
u = @SVector [0.0, 6.0, π/2, 0.0]
v = @SVector [0.0, 0.0, 0.0, phi_vel(m, u)]

sol = tracegeodesics(
    m,
    u,
    v,
    (0.0, 1200.0),
    μ = 1.0,
    abstol = 1e-10,
    reltol = 1e-10,
)
```

Tracing a couple of configurations and adding a boundary for the event horizon:

![stable-circ](/assets/stablecirc.svg)

**Note**: for performance reasons, the integrator terminates at $1.1 \, r_\text{s}$ to avoid wasting time on _doomed_ orbits.

Above are two orbits corresponding to circular $E$ and $L_z$ for $M=1.0$ and $a=0.0$. The dashed line has an orbit at $r = 5.6$, whereas the solid line $r = 6.0$, i.e. an orbit just _within_ the ISCO, and one _exactly on_ the ISCO. 

The integration is run for a total affine time of $1200\,\text{s}$, as up until around $1000\,\text{s}$ both orbits are circular. In theory, the inner orbit is still above $r_\text{mb}$ and should be able to remain a circular orbit, but due to the long integration times, we see the machine error contribute _enough_ that the orbit becomes unstable, collapses and falls into the black hole. Other values of $r < r_\text{ISCO}$ produce radically different unstable orbits, which all diverge in one way or another if the integrator is run for long enough.

### Orbit tracing regimes

We can investigate the stability of the integrator by setting the integration time to $10,000\,\text{s}$, and tracing all possible circular orbits for a specific $M$ and $a$. A _measure of deviation_ we employ is designed to tell us if the orbit collapsed inwards or outwards, and is calculated

\begin{equation}
    \mathcal{Q}_\text{d} = \frac{1}{r_\text{init}^2} \sum_{i=1} r_i^2,
\end{equation}

such that we can distinguish three cases:
- $\mathcal{Q}_\text{d} = 1$ implies the orbit is stable and remained at the same radius, 
- $\mathcal{Q}_\text{d} < 1$; the orbit collapsed and spiraled inwards (as in the above figure), and
- $\mathcal{Q}_\text{d} > 1$; the orbit spiralled outwards.

Note that the final two cases _are not_ sufficient conditions for the orbit to fall into or completely escape the black hole respectively. This measure instead tells us the majority of the orbit was either within or beyond the initial radius, but some orbits may have crossed back and forth.

![error-in-circular](/assets/circ-deviations.svg)

This figure is for $M=1.0$ and $a=-0.6$, since the effect of $a$ is only to translate the distribution left and right, and this configuration illustrates the regimes clearly. The radii considered are in steps of $\Delta r = 0.01$, and the pattern of $\mathcal{Q}_\text{d}$ is highly dependent on this, indicating that whether the orbit collapses or expands is due to machine precision. A few interesting things to note; the integrator is able to hold stable orbits all the way until the ISCO for long integration times. Between the marginally bound radius and the ISCO, deviations from the orbit are relatively small, which would indicate that the orbit remains circular for most of the integration time, before diverging -- which is what we would expect from hyperbolic energies when dealing with numeric simulations. Between the photon orbit and the marginally bound orbit is where we see orbits that either collapse or expand extremely far away.

Within the photon horizon, the integration terminates after one or two steps towards the event horizon, which is why we gradually see $\mathcal{Q}_\text{d}$ approach unity, as only a single point remains in the integration result.

Particularly pronounced are the different boundaries, causing clearly different general behavior for the integrator. To investigate this further, we'll vary $a$ and introduce a _measure of stability_, coined as simply the normalized squared sum of the residuals in $r$, expressed 

\begin{equation}
    \label{eq:Qs}
    \mathcal{Q}_\text{s} = \sqrt{ \frac{1}{N} \sum_{i=1} \left( \frac{r_i}{r_\text{init}} - 1 \right)^2},
\end{equation}

chosen to not penalize greater radii orbits for small fluctuations -- hence divide by initial radius $r_\text{init}$ -- and to not penalize orbits with many steps by dividing by total integrated points $N$. The motivation for using $\mathcal{Q}_\text{s}$ over $\mathcal{Q}_\text{d}$ is to be able to clearly differentiate regions where circular orbits are possible without penalizing eventual divergence in the same way $\mathcal{Q}_\text{d}$ does.

Applying this measure for a range of $a$ and $r$:

![error-in-circular](/assets/error-in-circular.svg)

Here, $\mathcal{Q}_\text{s} < 10^{-9}$ have been coloured transparent.

Overlayed are the different orbit boundaries. These boundaries align well with the different integration regimes, which is always nice to see. There are however a few anomalies, namely those orbits close to the event horizon for highly retrograde black holes, $a \rightarrow -0.998$, where we see very poor stability, as the geodesics are launched out of the system from small energetic perturbations. 

The small region of white within the event horizon in the bottom left corner is due to the integration being able to take a single step before it realizes any violation and terminates; it is, of course, a pure pathology.


Up until this point we have been able to use SageMath to derive results for $E$ and $L_z$ for circular orbits which hold all the way until $r_\text{ph}$. We have also confirmed that these values lead to circular orbits, and that the stability (and presence) of these orbits aligns with the different special boundaries. But what if we don't have good solutions for $E$ and $L_z$?

## Automated circular orbit discovery

Generalizing the last section; can we reliably and quickly find circular orbits without any _a priori_ knowledge of the velocity vectors, other than that we are searching in the equatorial plane?

The temptation here is to brute-force search over the $\dot{u}^\phi$ parameter space until some stability condition is met, but we can be more elegant with this by using [Optim.jl](https://julianlsolvers.github.io/Optim.jl/stable/), by optimizing univariate functions. Optim.jl will find a minima over a bounded range of any single argument function using the [golden-section search](https://en.wikipedia.org/wiki/Golden-section_search) algorithm.

For our purposes, the implementation needs a function which describes the degree of _circularity_ of our orbits, which is effectively just $\mathcal{Q}_\text{s}$ from eq. \eqref{eq:Qs}. 

The boilerplate is a few utility functions:

```julia
function trace_single_geodesic(m, u, vϕ)
    v = @SVector [0.0, 0.0, 0.0, vϕ]
    tracegeodesics(
        m, u, v, 
        # use a fairly long, but not excessively long integration time
        (0.0, 1000.0), 
        μ=1.0, 
        abstol=1e-10, 
        reltol=1e-10
    )
end

Qs(rs) = sqrt(sum((rs ./ rs[1] .- 1.0).^2) / length(rs))

# utility function to trace and evaluate some trial `vϕ`
function estimate_stability(m, u, vϕ)
    sol = trace_single_geodesic(m, u, vϕ)
    rs = selectdim(sol, 1, 6)
    Qs(rs)
end
```

Now we come to the crux, which is wrapping `estimate_stability` into a univariate function for Optim.jl. This is achieved through
```julia
using Optim

function find_vϕ_for_orbit(m, r_init; lower_bound=0.0, upper_bound=0.1)
    u = @SVector [0.0, r_init, deg2rad(90.0), 0.0]
    res = optimize(
        vϕ -> estimate_stability(m, u, vϕ),
        lower_bound,
        upper_bound,
        GoldenSection()
    )
    Optim.minimizer(res)
end
```

The limitation of the univariate optimizer is the bounding window of possible values; estimating the bounding box too low risks excluding the _true_ optimum, and estimating the bounding box too large risks local minima and bad high iteration count, which in turns leads to poor convergence. To try and mitigate this, we'll use a rolling window that continuously rescales these boundaries after each discovered orbit. 

We can bake an assumption into this rolling window, which is that stable orbits closer to the singularity are more likely to have a higher initial $\dot{u}^\phi$ than those further out, and solve the system backwards, starting with a conservative flat-space estimate for the velocity, and progressively work inwards towards the black hole. Half of this assumption is pretty robust, since the metrics we are dealing with are asymptotically flat, but it may not always be the case that orbits closer to the singularity rotate faster. 

In code, we express this rolling window as

```julia
function find_vϕ_for_orbit_range(m, rs; lower=0.0, upper=0.1)    
    map(reverse(rs)) do r
        vϕ = find_vϕ_for_orbit(m, r; lower_bound=lower, upper_bound=upper)
        # continuously adjust based on heuristic
        lower = vϕ * 0.9
        upper = vϕ * 2.0
        vϕ
    # call reverse twice so the order is maintained 
    end |> reverse
end
```

We can then solve for circular orbits over a range of $r$ with
```julia
m = BoyerLindquist(M=1.0, a=-0.4)
vϕ = find_vϕ_for_orbit(m, 10.0)

r_range = 1.1:0.1:10.0
vϕs = @time find_vϕ_for_orbit_range(m, r_range)
```

The solver takes on average 35 iterations (36 function calls) to converge to a minima of $\mathcal{Q}_\text{s}$ per orbit. Shorter integration times will also reliably converge in approximately the same number of function calls, with the routine being able to estimate stability after as little as a quarter of a full orbit reliably, up until approximately the ISCO. Within the ISCO, higher integration times are needed to accurately find solutions which illustrate the different orbital boundaries.

Calculating the energy for each orbit and plotting these again eq. \eqref{eq:Emu} shows very good agreement to within $10^{-9}$ for numerical accuracy:

![simulated-circs](/assets/simulated-circs.svg)

Of course, the region within the photon orbit is impossible to probe with this technique, since the integrator cannot handle this domain. The negative energy regimes are therefore also impossible to discover with this method. Nevertheless, it proves reliable for finding stable orbits, up to the initial domain of the bounding window. 

I'd also note here that this prototype code is able to those 100 stable orbits in just under a second, with approximately O($n$) scaling. Some aspects of this method could be optimized, but that is beyond the intention of this blog post.

## Discussion

We can get good approximations of analytic solutions using both SageMath and integrator-driven automated circular orbit discovery. Both methods are faced with limitations, however, in how they (fail to) probe within the photon radius. 

In the case of the SageMath, the negative energy description is entirely lost, and the angular momenta within the photon radius are inaccurate at best. When using the integrator, we cannot even probe within this radius reliably. It also has a practical limitation; SageMath is an unreliable solver at the best of times. Often the expressions it produces when invoking its `solve` function are a mess, if it is able to solve them at all. I do not know how _well_ this work will apply to metrics more complex than Kerr (which is the intended purpose after all). Though fortunately it seems that some fairly generic operations can go some distance to remedying this.

SageMath merits further investigation, and I will make the relevant SageMath notebooks for recreating the analysis in this blog available soon. 

The integrator suffers from an initial guess problem for the rolling window bounds on the possible range of $\dot{u}^\phi$, but other than that is able to work without issue directly from the input metric. It is also able to describe _enough_ of the space of circular orbits that the different orbital boundaries are visible, and will no doubt prove to be useful when studying other metrics.

Overall, I am pleased with how well the integrator performs for such orbits, and am confident these methods should be applicable to other spacetimes too.

### Future work

The next step is to apply this work to modified Kerr metrics, and to assemble a full pipeline for either deriving and wrapping symbolic expressions into Julia functions (related to cleaning the SageMath notebooks), or a pipeline for automatic orbit discovery in the integrator ecosystem (related to committing the above Julia code). This will then be used to calculate the redshift values for geometrically thin discs in the equatorial plane, to move the integrator in the direction of _observational results purely from metrics_.

Another avenue to later explore are the class of spherical orbits, which are not constrained to the equatorial plane, but nevertheless are able to maintain a constant radius.

A future post may also investigate ways of creating or even machine learning angular velocity tables for circular orbits, and how these could be shared and use in integration models.

## References

[^1]: J M Bardeen, W H Press, S A Teukolsky (1972). _Rotating black holes: locally nonrotating frames, energy extraction, and scalar synchrotron radiation_. APJ, **178**:347-169.

[^2]: C T Cunningham (1975). _The effects of redshifts and focusing on the spectrum of an accretion disk around a Kerr black hole_. APJ, **202**:788-802.

[^3]: B Carter (1968). _Global structure of the Kerr family of gravitational fields_. PHYS REV **174**:5.