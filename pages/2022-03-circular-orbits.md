@def title = "Circular orbits in static, axis-symmetric spacetimes"
@def description = "TODO"
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

The coverage of this blog post quite broad, but tries to investigate and re-derive results for the Kerr spacetime using an approach that can be programmatically computed for other static, axis-symmetric spacetimes. 

Circular orbits are of particular importance, as they are used to model processes in thin accretion discs. Circular orbit trajectories allows for the calculation of e.g. redshift values for thin discs when ray tracing using simple models. The work here will lead to automated stable orbit discovery with our integrator, and (hopefully) allow us to eventually simulate arbitrary disc profiles for different metrics.

## Introduction

- Notation: *italic* is basis-vector, serif is Euler's exponential $\e$.

From Bardeen et al. (1972)[^1], the standard, static, axis-symmetric, and asymptotically flat spacetime line element is

\begin{equation}
\label{eq:general-metric}
\d s^2  = 
    - \e^{2\nu} \d t^2 
    + \e^{2\psi} \left( \d\phi - \omega \d t \right)^2
    + \e^{2\mu_1} \d r^2
    + \e^{2\mu_2} \d \theta^2,
\end{equation}

where generally $\omega = -\frac{\tensor{g}{}{\phi t}}{\tensor{g}{}{\phi \phi}}$. In this context, the word *static* simply means *doesn't change with time* ($\frac{\partial}{\partial t} g_{\mu\nu} = 0$), axis-symmetric implies invariance under translations in $\phi$, and asymptotically flat means the metric reduces to the Minkowski metric, $\eta$, for $r \rightarrow \infty$. This metric becomes Kerr if

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

## Circular orbits in the Kerr metric

It is a known result that circular orbits exist in the equatorial plane of static axis-symmetric spacetimes, i.e. $\theta = \frac{\pi}{2}$ and $\dot{\theta} = 0$, shown by Carter (1968)[^3]. These circular orbits are further constrained by motion in the radial components,

\begin{equation}
    \label{eq:ur_dotur}
    \dot{u}^r = 0, \quad \text{and} \quad \frac{\partial}{\partial r}\dot{u}^r = 0.
\end{equation}

We can use the invariance of four velocities to express these constraints further,

\begin{equation}
    \label{eq:velconstraint}
    \tensor{g}{}{\sigma \nu} \tensor{\dot{u}}{\sigma}{} \tensor{\dot{u}}{\nu}{} = - \mu^2,
\end{equation}

where, for time-like geodesics, $\mu = 1$. In order to rearrange for the radial component, this equation may be expanded,

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

Writing this equation this way allows the usual identification for energy $\tensor{\dot{u}}{}{t} = - E$ and angular momentum $\tensor{\dot{u}}{}{\phi} = L_z$ from the covariant momenta,

\begin{equation}
    \therefore \quad \left(\tensor{\dot{u}}{r}{} \right)^2 = 
        \frac{-1}{ \tensor{g}{}{rr} } \left(
            \tensor{g}{tt}{} E^2
            - 2 \tensor{g}{t\phi}{} E L_z
            + \tensor{g}{\phi\phi}{} L_z^2
            + \mu^2
        \right).
\end{equation}

The constraints from eq. \eqref{eq:ur_dotur} may now be calculated from this expression. 

\begin{equation}
    \tensor{\dot{u}}{r}{} = \frac{\pm 1}{ \e^{\mu_1} } \sqrt{
        \e^{-2\nu} \left(
            E + \omega L_z
        \right)^2
        - \e^{-2\Psi} L_z^2
        - \mu^2
    }
    = 0.
\end{equation}

If we rephrase the above constraint in terms of effective potentials, we find the expression may be written much more simply as

\begin{equation}
    \tensor{\dot{u}}{r}{} = \pm \sqrt{ \mathcal{R} \, \e^{-2\mu_1} } = 0.
\end{equation}

Our constraints are now $\mathcal{R} = 0$ and $\frac{\partial}{ \partial r}\mathcal{R} = 0$, the latter of which is directly

\begin{equation}
\label{eq:sage-dont-like}
    \left( E + \omega L_z \right)^2 \frac{\partial}{\partial r} \e^{-2\nu}
    + 2 \left( E + \omega L_z \right) \e^{-2\nu} L_z \frac{\partial}{\partial r} \omega
    - L_z^2 \frac{\partial}{\partial r} \e^{-2 \Psi}
    = 0.
\end{equation}

The term $\mathcal{R} \frac{\partial}{\partial r} \e^{-2\mu_1}$ is omitted, as it is manifestly zero, since $\mathcal{R} = 0$. Implicitly it is assumed that $E$ and $L_z$ are constant with respect to $r$ over an orbit. Indeed, this assumption is valid as $E$ and $L_z$ are constants of motion for any given geodesic, whereas the $E(r)$ and $L_z(r)$ which we later solve for are functions on the phase space. This distinction is important to note.  

Here we will specialise to the Kerr metric, with subsequent manipulations performed using SageMath, and the occasional small bit of hand manipulation. Substituting in the metric components from eq. \eqref{eq:metric-comp}, we obtain

\begin{equation}
    \label{eq:constraint-1}
    \mathcal{R} = 2 M r {\left(\mu^{2} r^{2} + {\left(E a - L_z\right)}^{2}\right)} + {\left({\left(E^{2} - \mu^{2}\right)} {\left(a^{2} + r^{2}\right)} - L_z^{2}\right)} \Sigma = 0
\end{equation}
\begin{equation}
    \label{eq:constraint-2}
    \frac{\partial}{\partial r} \mathcal{R} = 6 \, M \mu^{2} r^{2} + 2 \, {\left(E a - L_z\right)}^{2} M + 2 \, {\left(E^{2} - \mu^{2}\right)} r \Sigma + {\left({\left(E^{2} - \mu^{2}\right)} {\left(a^{2} + r^{2}\right)} - L_z^{2}\right)} \frac{\partial}{\partial r}\Sigma = 0.
\end{equation}

It is worth noting that $\mu$ has reappeared in the derivative term. This is due to the complexity of eq. \eqref{eq:sage-dont-like} when used in SageMath -- instead, SageMath prefers to just take the derivative with respect to $r$ of eq. \eqref{eq:constraint-1}. Both will eventually solve the system identically, but this version is more friendly for a blog post.

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

Our constrains in eq. \eqref{eq:ur_dotur} are now eq. \eqref{eq:constraint_1} and \eqref{eq:constraint_2} respectively, which may be solved simultaneously to find $E$ and $L_z$ satisfying circular orbits. In SageMath, this is fiddly, but generally quite straight forward to solve. We find there are four solutions in total, representing the intersection of the positive and negative energy regimes with prograde and retrograde orbits,

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

where $\pm$ denotes prograde and retrograde respectively, and

\begin{equation}
\label{eq:Lzmu}
    \frac{Lz}{\mu} =
    \frac{E}{\mu} a
    \pm \frac{r}{\sqrt{M}} \sqrt{
        1 - r + r \left(\frac{E}{\mu} \right)^2
    }.
\end{equation}

At the time of writing, this is the most simplified form I have been able to manipulate these equations into. 

As an additional note, with SageMath, many of the preceding steps in preparing the constraints are not actually necessary, but would provide results that are unwieldy for such a blog post. At the end of this post I have included links to my relevant Jupyter notebooks for recreating these solutions.

### A comparison to Bardeen (1972)

Bardeen (1972) quotes

\begin{equation}
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

which are substantially different from eq. \eqref{eq:Emu} and eq. \eqref{eq:Lzmu}. We can try to find the degree of similarity by plotting the equations with $r > 0.1$ to avoid the central singularity:

![energy-diff](/assets/angmom-energy.svg)

The plots for energy (left panel) are in good agreement up until the energy becomes negative in the Bardeen case (such energies are allowed, e.g. in the [Penrose process](https://en.wikipedia.org/wiki/Penrose_process)). The SageMath derivation fails in this regime, as it is the square root of a squared quantity, which loses sign information. Indeed, we do have a negative energy solution from SageMath, but we would require _a priori_ knowledge of which solution to apply in which domain.

The consequence of the energy sign is even more pronounced in the angular momentum (right panel), as the angular momentum is calculated from the energy. For the $a=0.998$ case, the plots are still in good agreement, but the divergence for $a=0.0$ is pronounced, and worse still for $a=-0.988$ (not depicted, as the diagram became difficult to annotate). Some of the features of the SageMath $L_z$ curves do not disappear when using a sign-concious approach, such as the sharp cusp at $r\sim2$.

Fortunately, the differences in the two solutions is for our purposes pedagogical, since the divergence lies within the _marginally bound_ orbits, at which radii circular orbits are prohibited.


## Special orbit boundaries

There are many special boundaries in the Kerr metric, including the event-horizon and the ergosphere, setting the boundary of the static limit, at which the frame dragging effect appears which originally motivated the LNRF. For purely circular orbits, we consider two of these boundaries and the respective possible orbits within them, as they are arguably most relevant to our study -- these are the _innermost stable circular orbit_ (ISCO) and the _marginally bound orbits_.

### ISCO

Stability requires the energy be a minimum; this may be rephrased as

\begin{equation}
    \frac{\partial^2}{\partial r^2} \mathcal{R} \geq 0
\end{equation}

### Marginally bound orbits

All unbound circular orbits have hyperbolic energies $\frac{E}{\mu} > 1$. An unbound orbit is then a parabolic orbit at $r > r_{\text{mb}}$, with the _marginally bound_ radius defined, with the Kerr metric,

\begin{equation}
    r_{\text{mb}} = 2M \mp a + 2 \sqrt{M^2 \mp aM },
\end{equation}

corresponding to first (infalling) radius at which $\frac{E}{\mu} = 1$.

### Photon orbit

## Tracing circular orbits

Tracing stable circular orbits when $E$ and $L_z$ is known is relatively simple. We only need to calculate $\dot{u}^\phi$, since the integrator automatically constrains $\dot{u}^t$ via eq. \eqref{eq:velconstraint}. We find

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

This defined a utility function for calculating the initial radial velocity component for a given metric and position. The actual tracing code is then:

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

Above are two orbits corresponding to circular $E$ and $L_z$ for $M=1.0$ and $a=0.0$. The dashed line has an orbit at $r = 5.6$, whereas the solid line $r = 6.0$, i.e. an orbit just _within_ the ISCO, and one _exactly on_ the ISCO. 

**Note**: for performance reasons, the integrator terminates at $1.1 \, r_\text{s}$ to avoid wasting time on _doomed_ orbits.

The integration is run for a total affine time of $1200\,\text{s}$, as up until around $1000\,\text{s}$ both orbits are circular. In theory, the inner orbit is still above $r_\text{mb}$ and should remain a circular orbit, but due to the long integration times, we see the machine error contribute _enough_ that the orbit collapses and falls into the black hole. Other values of $r < r_\text{ISCO}$ produce radically different unstable orbits, which all collapse in one way or another if the integrator is run for long enough.

### Orbit tracing regimes

We can investigate the stability of the integrator by setting the integration time to $10,000\,\text{s}$, and tracing all possible circular orbits for a specific $M$ and $a$. A _measure of deviation_ we employ is designed to tell us if the orbit collapsed inwards or outwards, and is calculated

\begin{equation}
    \mathcal{Q}_\text{d} = \frac{1}{r_\text{init}^2} \sum_{i=1} r_i^2,
\end{equation}

such that we can distinguish three cases:
- $\mathcal{Q}_\text{d} = 1$ implies the orbit is stable and remained at the same radius, 
- $\mathcal{Q}_\text{d} < 1$; the orbit collapsed and spiraled inwards (as in the above figure), and
- $\mathcal{Q}_\text{d} > 1$; the orbit spiralled outwards.

Note that the final two cases are not sufficient conditions to state that the orbit fell into or escaped the black hole respectively. This measure instead tells us the majority of the orbit was either within or beyond the initial radius, but it also may well have crossed back and forth.

![error-in-circular](/assets/circ-deviations.svg)

This figure is for $M=1.0$ and $a=-0.6$, since the effect of $a$ is only to translate the distribution left and right, and this configuration illustrates the regimes clearly. The radii considered are in steps of $\Delta r = 0.01$, and the pattern of $\mathcal{Q}_\text{d}$ is highly dependent on this, indicating that whether the orbit collapses or expands is due to machine precision. A few interesting things to note; the integrator is able to hold stable orbits well all the way until the ISCO. Between the marginally bound radius and the ISCO, deviations from the orbit are relatively small, which would indicate that the orbit remains circular for most of the integration time, before diverging -- which is what we would expect from hyperbolic energies when dealing with numeric simulations. Between the photon orbit and the marginally bound orbit is where we see orbits that either collapse or expand far away.

Below the photon horizon the integration terminates after one or two steps towards the event horizon, which is why we gradually see $\mathcal{Q}_\text{d}$ approach unity, as only a single point remains in the integration result.

Particularly pronounced are the different boundaries, causing clearly different general behavior for the integrator. To investigate this further, we'll vary $a$ and introduce a _measure of stability_, coined as simply the normalized squared sum of the residuals in $r$, expressed 

\begin{equation}
    \label{eq:Qs}
    \mathcal{Q}_\text{s} = \sqrt{ \frac{1}{N} \sum_{i=1} \left( \frac{r_i}{r_\text{init}} - 1 \right)^2},
\end{equation}

chosen to not penalize greater radii orbits for small fluctuations -- hence divide by initial radius $r_\text{init}$ -- and to not penalize orbits with many steps by dividing by total integrated points $N$. The motivation for using $\mathcal{Q}_\text{s}$ over $\mathcal{Q}_\text{d}$ is to be able to be able to clearly differentiate regions where circular orbits are possible without penalizing eventual divergence in the same way $\mathcal{Q}_\text{d}$ does.

Applying this measure for a range of $a$ and $r$:

![error-in-circular](/assets/error-in-circular.svg)

Overlayed are the different orbit boundaries. These boundaries align well with the different integration regimes, which is always nice to see. There are however a few anomalies, namely those orbits close to the event horizon for highly retrograde black holes, $a \rightarrow -0.998$, where we see very poor stability, as the geodesics are launched out of the system from small energetic perturbations. The small region within the event horizon of white is due to the integration being able to take a step before it terminates, but is not physical.

## Automated circular orbit discovery

Generalizing the last section; can we reliably and quickly find circular orbits without any _a priori_ knowledge of the velocity vectors, other than that we are searching in the equatorial plane?

The temptation here is to brute for search the $\dot{u}^\phi$ parameter space until some stability condition is met, but we can be more elegant with this by using [Optim.jl](https://julianlsolvers.github.io/Optim.jl/stable/), for optimizing univariate functions. Optim.jl is guaranteed to find a minima over a bounded range of any single argument function using the [golden-section search](https://en.wikipedia.org/wiki/Golden-section_search) algorithm.

For our purposes, we just need a function which described the degree of stability for our circular orbits, which is just $\mathcal{Q}_\text{s}$ from eq. \eqref{eq:Qs}. The boilerplate we need are a few utility functions:

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

Qs(rs) = sqrt(sum((rs .- rs[1]).^2) / length(rs))

# utility function to trace and evaluate some trial `vϕ`
function estimate_stability(m, u, vϕ)
    sol = trace_single_geodesic(m, u, vϕ)
    rs = selectdim(sol, 1, 6)
    Qs(rs)
end
```

Now we come to the crux, which is wrapping the this evaluating function for Optim.jl:
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

The limitation of the univariate optimizer is the bounding window of possible values; we'll use a rolling window that continuously rescales these boundaries after each successive integration. Estimating the bounding box too low risks excluding the true optimum, and estimating the bounding box too large risks local minima and bad high iteration count, which in turns leads to poor convergence.

We can bake an assumption into this rolling window, which is that stable orbits closer to the singularity are more likely to have a higher initial $\dot{u}^\phi$ than those further out, and solve the system backwards, starting with a conservative flat-space estimate for the velocity, and progressively work inwards towards the black hole. In code, this is

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

The solver takes on average 35 iterations (36 function calls) to converge to a minima of $\mathcal{Q}_\text{s}$. Shorter integration times will converge similarly, and the routine seems to be able to estimate stability after as little as a quarter of a full orbit reliably, up until approximately the ISCO. Within the ISCO, higher integration times are needed to accurately find solutions.

Calculating the energy for each orbit and plotting these again eq. \eqref{eq:Emu} shows very good agreement to within $10^{-9}$ for numerical accuracy:

![simulated-circs](/assets/simulated-circs.svg)

Of course, the region within the photon orbit is impossible to probe with this technique, since the integrator cannot handle this domain. Never the less, this method proves reliable for finding stable orbits, up to the accuracy of the initial bounding window. 

I'd also note here that this prototype code is able to those 100 stable orbits in just under a second, with approximately O($n$) scaling. Some aspects of this method could be optimized, but that is beyond the intention of this blog post.

## Discussion

### SageMath limitations

SageMath is an unreliable solver at the best of times. Often the expressions it produces when invoking its `solve` function are a mess, if it is able to solve them at all. I do not know how _well_ this work will apply to metrics more complex than Kerr (which is the intended purpose after all). Fortunately, it seems the integrator is able to find stable (circular) orbits without any _a priori_ knowledge, which is very ideal.

A future post will investigate ways of creating or even machine learning velocity banks for sharing and use in integrating.

### Future work

The next steps are to apply this work to modified Kerr metrics, and to assemble a full pipeline for either deriving and wrapping symbolic expressions into Julia functions, or for automatically discovering these stable orbits. This will then be used to calculate the redshift values for geometrically thin discs in the equatorial plane, using only a metric input.

Another avenue to later explore are the class of spherical orbits, which are not constrained to the equatorial plane, but never the less are able to maintain a constant radius.

## References

[^1]: J M Bardeen, W H Press, S A Teukolsky (1972). _Rotating black holes: locally nonrotating frames, energy extraction, and scalar synchrotron radiation_. APJ, **178**:347-169.

[^2]: C T Cunningham (1975). _The effects of redshifts and focusing on the spectrum of an accretion disk around a Kerr black hole_. APJ, **202**:788-802.

[^3]: B Carter (1968). _Global structure of the Kerr family of gravitational fields_. PHYS REV **174**:5.