@def title = "Plunging photons out of stable orbits"
@def description = "I'd like to be able to calculate programmatically key results related to circular orbits in generic spacetimes. To this end, I want to see how feasible a computer aided approach is, using only a handful of constraints, in re-deriving known analytic solutions, and compare this to automated circular orbit discovery using our geodesic integrator, focusing on the Kerr metric."
@def image = "/assets/plunging/dilaton-cover.png"
@def date = "09/05/2022"
@def tags = ["2022", "geodesics", "julia", "relativity"]

# Plunging photons out of stable orbits

Following work from [a previous blog post](./2022-03-circular-orbits) concerning circular orbits in the Kerr metric, this post will investigate the region within the innermost stable circular orbit (ISCO). I'd also like to re-derive results for general static, axis-symmetric metrics, such that the energy $E$ and angular momentum $L_z$ from the previous blog post may be more compactly expressed. This will differ from the SageMath method by introducing the _Keplerian frequency_, and manipulating the equations by hand. 

I aim to calculate velocity vectors of massive geodesics in the so-called _plunging region_, i.e. the region within the ISCO, and use this to calculate redshifts in our class of spacetimes. In the literature, these results rely on having explicit first-order equations of motion. Since $L_z$, Carter's constant $Q$, and $E$ are known for orbits at the ISCO, the plunging region trajectories are calculated by substituting these into the first-order equations, thereby lifting the condition on the radial component of the velocity. 

In practice, this is equivalent to pushing the photon an _infinitesimal_ amount towards the black hole. Consequently, if these first-order equations are unknown, we can instead create numerical approximations for the analytic solutions by tracing this _slightly offset_ trajectory with [Gradus.jl](https://github.com/astro-group-bristol/Gradus.jl) (I finally settled on a name). 

## Stable orbits in static, axis-symmetric metrics

Following the previous blog post, and elaborating on the method of Johannsen (2013)[^1], we start by restricting ourselves to $\theta = \frac{\pi}{2}$, with $\dot{\theta} = 0$, and require, for stable orbits, that

\begin{equation}
    \label{eq:ur_dotur}
    \dot{u}^r = 0, \quad \text{and} \quad \frac{\partial}{\partial r}\dot{u}^r = 0.
\end{equation}

Here, as in the previous post, the dot denotes differentiation with respect to the affine parameter $\lambda$.

We write the geodesic equation

\begin{equation}
    \label{eq:geodesic_eq}
    \frac{\d}{\d \lambda} \dot{u}^\sigma + \tensor{\Gamma}{\sigma}{\nu \varepsilon} \dot{u}^\nu \dot{u}^\varepsilon = 0,
\end{equation}

and use the definition of the Christoffel symbols,

\begin{equation}
\tensor{\Gamma}{\sigma}{\nu \varepsilon} = \frac{1}{2} g^{\sigma \alpha} \left(
    \frac{\partial}{\partial u^\nu} g_{\alpha \varepsilon} 
    + \frac{\partial}{\partial u^\varepsilon} g_{\alpha \nu}
    - \frac{\partial}{\partial u^\alpha} g_{\nu \varepsilon}
\right).
\end{equation}

to expand the geodesic equation as

\begin{equation}
    g_{\sigma \alpha} \frac{\d}{\d \lambda} \dot{u}^\sigma
    + \frac{1}{2} \left(
        \frac{\partial}{\partial u^\nu} g_{\alpha \varepsilon} 
        + \frac{\partial}{\partial u^\varepsilon} g_{\alpha \nu}
        - \frac{\partial}{\partial u^\alpha} g_{\nu \varepsilon}
    \right) \dot{u}^\nu \dot{u}^\varepsilon
    = 0.
\end{equation}

Notice the indices of the terms in bracket may be interchanged in the sum over the velocity vectors given the symmetry of the tensors involved, permitting us to combine terms and rearrange

\begin{equation}
    g_{\sigma \alpha} \frac{\d}{\d \lambda} \dot{u}^\sigma
    + 
    \dot{u}^\nu \dot{u}^\varepsilon \frac{\partial}{\partial u^\nu} g_{\alpha \varepsilon}
    = \frac{1}{2} \dot{u}^\mu \dot{u}^\beta \frac{\partial}{\partial u^\alpha} g_{\mu \beta}.
\end{equation}

Next, we (ab)use the chain rule to express

\begin{equation}
    \dot{u}^\varepsilon \frac{\d}{\d \lambda} g_{\alpha \varepsilon}
    = \dot{u}^\varepsilon \frac{\d}{\d \lambda} u^\nu \frac{\partial}{\partial u^\nu} g_{\alpha \varepsilon}
    = \dot{u}^\varepsilon \dot{u}^\nu \frac{\partial}{\partial u^\nu} g_{\alpha \varepsilon},
\end{equation}

and therefore

\begin{equation}
    g_{\sigma \alpha} \frac{\d}{\d \lambda} \dot{u}^\sigma
    + 
    \dot{u}^\varepsilon \frac{\d}{\d \lambda} g_{\alpha \varepsilon}
    = \frac{1}{2} \dot{u}^\mu \dot{u}^\beta \frac{\partial}{\partial u^\alpha} g_{\mu \beta},
\end{equation}

or, using the product rule,

\begin{equation}
    \frac{\d}{\d \lambda} \left( g_{\sigma\alpha} \dot{u}^\sigma \right)
    = \frac{1}{2} \dot{u}^\mu \dot{u}^\beta \frac{\partial}{\partial u^\alpha} g_{\mu \beta}.
\end{equation}

Examining just the radial derivative components ($\alpha = r$), along with the condition that $\dot{u}^r = 0$, forces the LHS to disappear for our class of metrics, with the only non-zero terms of the RHS for stable orbits being

\begin{equation}
    \label{eq:geodesic_rewrite}
    0 = \frac{\partial}{\partial r} g_{tt} (\dot{u}^t)^2
    + 2 \frac{\partial}{\partial r} g_{t\phi} \, \dot{u}^t \dot{u}^\phi
    +\frac{\partial}{\partial r} g_{\phi \phi} (\dot{u}^\phi)^2.
\end{equation}

This is just an alternative method to achieve a result from the previous blog post, which I wanted to derive motivated by the Johannsen paper.

### Keplerian frequency

We introduce the _Keplerian frequency_, derived from Kepler's third law, which is the rate at which the angular coordinate changes with coordinate time. As an aside, the Keplerian frequency also allows us to calculate the rotation limits of various astrophysical bodies, as a star rotating faster than this frequency would promote particles into higher orbits, effectively ripping itself apart (something for another post).

The Keplerian frequency for our purposes may be defined 

\begin{equation}
    \label{eq:keplerian_freq}
    \Omega_\phi := \frac{\dot{u}^\phi}{\dot{u}^t},
\end{equation}

and therefore Eq. \eqref{eq:geodesic_rewrite} may be expressed

\begin{equation}
    \Omega_\phi^2 \, \frac{\partial}{\partial r} g_{\phi \phi}
    + 2 \Omega_\phi \, \frac{\partial}{\partial r} g_{t \phi}
    + \frac{\partial}{\partial r} g_{t t}
    = 0,
\end{equation}

inviting us to solve for the Keplerian frequency

\begin{equation}
    \label{eq:keplerian_freq_sol}
    \Omega_\phi 
    =
    \left(
        - \frac{\partial}{\partial r} g_{t \phi} \pm
        \sqrt{
            \left( \frac{\partial}{\partial r} g_{t \phi} \right)^2 
            - \frac{\partial}{\partial r} g_{t t} \frac{\partial}{\partial r} g_{\phi \phi}
        }
    \,\right)
    \left( \frac{\partial}{\partial r} g_{\phi \phi} \right)^{-1}.
\end{equation}

As an exercise, substituting in the Kerr metric terms yields

\begin{equation}
    \Omega_\phi = \frac{1}{ M a \pm r^{3/2} }.
\end{equation}

### Energy and angular momentum

To solve for the energy and angular momentum of a given stable circular orbit, we use that

\begin{equation}
    g^{\sigma \nu} \dot{u}_\sigma \dot{u}_\nu = - \mu^2,
\end{equation}

expanding components to find

\begin{equation}
    \label{eq:ur_squared}
    \quad \left(\tensor{\dot{u}}{r}{} \right)^2 = 
        \frac{-1}{ \tensor{g}{}{rr} } \left(
            \tensor{g}{tt}{} E^2
            - 2 \tensor{g}{t\phi}{} E L_z
            + \tensor{g}{\phi\phi}{} L_z^2
            + \mu^2
        \right),
\end{equation}

where we have identified covariant momenta $\dot{u}_t = -E$ and $\dot{u}_\phi = L_z$, and used to the freedom to raise and lower indices under contraction. We will write an effective potential $(\dot{u}^r)^2 = V_\text{eff}(r)$, and solve 

\begin{equation}
    \label{eq:eff_pot_conditions}
    V_\text{eff}(r) = 0, \quad \text{and} \quad \frac{\d}{\d r} V_\text{eff}(r) = 0,
\end{equation}

simultaneously for $E$ and $L_z$. This approach is very similar to last time, except now we will use the Keplerian frequency in Eq. \eqref{eq:keplerian_freq}, with the observation that $\dot{u}^\phi = g^{\phi \phi} L_z - g^{t\phi} E$ and $\dot{u}^t = g^{t \phi} L_z - g^{t t} E$, inferred from the operation of the metric. Consequently, the Keplerian frequency is

\begin{equation}
    \label{eq:kepler_e_l}
    \Omega_\phi = \frac{g^{\phi \phi} L_z - g^{t\phi} E}{g^{t \phi} L_z - g^{t t} E}.
\end{equation}

Johannsen uses the definition of $\Omega_\phi$ in Eq. \eqref{eq:keplerian_freq_sol} along with Eq. \eqref{eq:eff_pot_conditions} and Eq. \eqref{eq:kepler_e_l} to derive the energy and angular momentum, but I was unable to recreate this. Instead, I came up with method to derive these using only the latter two and the metric condition

\begin{equation}
    g^{\mu \nu} g_{\mu \rho} = \tensor{\delta}{\nu}{\rho},
\end{equation}

manifest for any metric tensor. In particular, for the class of static, axis-symmetric metrics, with $\nu = \rho$ one obtains

\begin{equation}
    \label{eq:g_lower_1}
    g^{t \phi} \, g_{t \phi} + g^{\phi \phi} \, g_{\phi \phi}
    = g^{t \phi} \, g_{t \phi} + g^{t t} \, g_{t t} = 1,
\end{equation}

and with $\nu \ \cancel{=} \ \rho$,

\begin{equation}
    \label{eq:g_lower_2}
    g^{t \phi}\, g_{tt} + g^{\phi\phi} \, g_{t \phi}
    = g^{t \phi}\, g_{\phi \phi} + g^{tt}\, g_{t \phi} = 0.
\end{equation}

We begin by re-arranging Eq. \eqref{eq:kepler_e_l} into

\begin{equation}
    L_z = \frac{g^{t \phi} - g^{t t} \Omega_\phi}{g^{\phi \phi} - g^{t\phi} \Omega_\phi} E.
\end{equation}

Multiplying in the fraction by $g_{t\phi}$

\begin{equation}
    L_z = \frac{g_{t\phi} \, g^{t \phi} - g_{t\phi} \, g^{t t} \Omega_\phi}{g_{t\phi} \, g^{\phi \phi} - g_{t\phi} \, g^{t\phi} \Omega_\phi} E,
\end{equation}

and using Eq. \eqref{eq:g_lower_2} to substitute relevant symbols, we find

\begin{equation}
    \label{eq:lz_e_a_b}
    L_z = - \left( \frac{g_{t \phi} + g_{\phi \phi} \Omega_\phi}{g_{tt} + g_{t \phi} \Omega_\phi} \right) E = - \, \frac{\mathcal{A}}{\mathcal{B}}\, E,
\end{equation}

where we use $\mathcal{A}$ and $\mathcal{B}$ for convenience to denote the numerator and denominator respectively. Substituting this expression for $L_z$ into Eq. \eqref{eq:eff_pot_conditions} yields

\begin{equation}
    E^2 \left[
        g^{tt} + 2 g^{t\phi} \left( \frac{\mathcal{A}}{\mathcal{B}} \right) + g^{\phi \phi}\left( \frac{\mathcal{A}}{\mathcal{B}} \right)^2
    \right] + \mu^2 = 0,
\end{equation}

or equivalently

\begin{equation}
    \frac{E^2}{\mu^2} = \frac{
        \mathcal{B}^2
    }{
        - g^{tt} \mathcal{B}^2  - 2 g^{t\phi} \mathcal{A} \mathcal{B} - g^{\phi \phi} \mathcal{A}^2
    },
\end{equation}

and, by symmetry of Eq. \eqref{eq:lz_e_a_b}, one can immediately write

\begin{equation}
    \frac{L_z^2}{\mu^2} = \frac{
        \mathcal{A}^2
    }{
        - \left( g^{tt} \mathcal{B}^2  + 2 g^{t\phi} \mathcal{A} \mathcal{B} + g^{\phi \phi} \mathcal{A}^2 \right)
    }.
\end{equation}

Examining this common denominator by expanding all $\mathcal{A}$ and $\mathcal{B}$ and collecting powers of $\Omega_\phi$ gives us

\begin{align}
    \Omega_\phi^2 &\left( 
        g^{tt} \, g_{t\phi} \, g_{t\phi} + 2 g^{t\phi} \, g_{t\phi} \, g_{\phi\phi} + g^{\phi\phi} \, g_{\phi\phi} \, g_{\phi\phi}
    \right)\\
    +\ 2 \,\Omega_\phi &\left( 
        g^{\phi\phi} \, g_{t\phi} \, g_{t\phi} + g^{t \phi} \, g_{\phi\phi} \, g_{tt} + g^{tt} \, g_{tt} \, g_{t\phi} + g^{\phi \phi} \, g_{\phi\phi} \, g_{t\phi}
    \right)\\
    +\ \phantom{2 \,\Omega_\phi} &\left( 
        g^{tt} \, g_{tt} \, g_{tt} + 2g^{t\phi} \, g_{t\phi} \, g_{tt} + g^{\phi\phi} \, g_{\phi\phi} \, g_{\phi\phi}
    \right),
\end{align}

dropping the negative sign temporarily.

The next steps are pretty methodological, so I will demonstrate the mechanics for $\Omega_\phi^2$ terms, and then state the results for the other powers in the interest of brevity. First we use Eq. \eqref{eq:g_lower_1} to simplify the first term

\begin{equation}
    \Omega_\phi^2 \left(
        -g^{t\phi} \, g_{t\phi} \, g_{\phi\phi} + 2 g^{t\phi} \, g_{t\phi} \, g_{\phi\phi} + g^{\phi\phi} \, g_{\phi\phi} \, g_{\phi\phi}
    \right),
\end{equation}

which may then be combined to write

\begin{equation}
    \Omega_\phi^2 \, g_{\phi\phi} \left(
        g^{t\phi} \, g_{t\phi} + g^{\phi\phi} \, g_{\phi\phi}
    \right),
\end{equation}

and, via Eq. \eqref{eq:g_lower_2},

\begin{equation}
    \implies \Omega_\phi^2 \, g_{\phi\phi}.
\end{equation}

Similarly for the other powers. After working our way through this piece of algebra, we obtain

@@result
\begin{equation}
    \label{eq:stable_energy}
    \frac{E}{\mu} = \pm \frac{g_{tt} + g_{t\phi} \, \Omega_\phi}{
        \sqrt{
            -g_{\phi \phi} \, \Omega_\phi^2 - 2 g_{t\phi} \, \Omega_\phi - g_{tt}
        }
    },
\end{equation}

and

\begin{equation}
    \label{eq:stable_angmom}
    \frac{L_z}{\mu} = \pm \frac{g_{t\phi} + g_{\phi\phi} \, \Omega_\phi}{
        \sqrt{
            -g_{\phi \phi} \, \Omega_\phi^2 - 2 g_{t\phi} \, \Omega_\phi - g_{tt}
        }
    }.
\end{equation}
@@

### Tracing orbits with the analytic constants

To verify these results, we can attempt to trace orbits around a Kerr black hole. We first define the relevant formulae

```julia
using Gradus
using StaticArrays

# implementation for energy
# fix ourselves to the negative energy solution arbitrarily
function energy(g, Ωϕ)
    -(g[1] + g[5] * Ωϕ) / √(-g[1] - 2g[5]*Ωϕ - g[4]*Ωϕ^2)
end

# implementation for angular momentum
# include a flag for whether we want the prograde or retrograde orbit
function angmom(g, Ωϕ, prograde::Bool)
    res = (g[5] + g[4] * Ωϕ) / √(-g[1] - 2g[5]*Ωϕ - g[4]*Ωϕ^2)
    if prograde
        res
    else
        -res
    end
end

# keplerian frequency
# requires choice of positive or negative domain
function Ω(m, rθ, positive::Bool)
    # require derivatives with respect to radial coordinate
    jacs = Gradus.GeodesicTracer.metric_jacobian(m, rθ)
    ∂rg = jacs[:, 1]

    Δ = √(∂rg[5]^2 - ∂rg[1] * ∂rg[4])
    if positive 
        -(∂rg[5] + Δ) / ∂rg[4]
    else
        -(∂rg[5] - Δ) / ∂rg[4]
    end
end

# calculate all relevant metric quantities for stable circular orbit system
function stable_system(m, r; contra_rotating=false, prograde=true)
    let rθ = @SVector([r, π/2])
        g = Gradus.GeodesicTracer.metric_components(m, rθ)
        ginv = Gradus.GeodesicTracer.inverse_metric_components(g)

        Ωϕ = Ω(m, rθ, contra_rotating)
        E = energy(g, Ωϕ)
        L = angmom(g, Ωϕ, prograde)

        (g,ginv,E,L)
    end
end
```

This code is all relatively straight forward, using Gradus.jl to calculate the metric components for us, and just using them where needed. We can then write a utility function for creating stable 4-velocities:

```julia
function fourvelocity(m, r; kwargs...)
    g, ginv, E, L = stable_system(m, r; kwargs...) 
    vt = -E * ginv[1] + L * ginv[5]
    vϕ = -E * ginv[5] + L * ginv[4]

    @SVector[vt, 0.0, 0.0, vϕ]
end
```

and trace as usual

```julia
# choice of metric
m = BoyerLindquistAD(M=1.0, a=0.0)

# range of radii
rs = range(Gradus.isco(m), 10.0, 7)

# generate position and init velocities
us = [@SVector([0.0, r, deg2rad(90), 0.0]) for r in rs]
vs = map(r -> fourvelocity(m, r), rs)

# trace the geodesics for a fairly long duration to ensure stability
simsols = @time tracegeodesics(
    m,
    us, vs,
    (0.0, 2000.0);
    abstol=1e-10, reltol=1e-10,
    μ=1.0
)
```

![analytic-stable-orbits](/assets/plunging/analytic-stable-orbits.svg)

Excellent, the orbits are indeed perfectly stable!

## Plunging orbits

The plunging region is the region of unstable circular orbits, where the energies are hyperbolic -- any small deviation from the trajectory will cause the particle to spiral inwards towards the black hole. As stated in the introduction, these trajectories are traditionally calculated using analytic first-order equations of motion, by substituting $E$ and $L_z$ at the ISCO, along with the condition of equatorial orbits $Q=0$, and calculating $\dot{u}^\nu$.

### Four-velocity calculations

Instead, we will attempt to calculate $\dot{u}^\nu$ from metric terms at the ISCO.

We need minimally $\dot{u}^r$ and $\dot{u}^\phi$, since $\dot{u}^\theta = 0$ and $\dot{u}^t$ is constrained by the integrator. We already have $\dot{u}^\phi$ from the previous section, and $\dot{u}^r$ may be found from Eq. \eqref{eq:ur_squared}. 

A possible implementation for our plunging four-velocity may then be

```julia
function plunging_fourvelocity(m, r; kwargs...) where {T}
    g, ginv, E, L = stable_system(m, r; kwargs...) 
    vt = -E * ginv[1] + L * ginv[5]
    vϕ = -E * ginv[5] + L * ginv[4]
    
    nom = ginv[1] * E^2 - 2ginv[5] * E * L + ginv[4] * L^2
    denom = -g[2]

    @SVector[vt, -sqrt(abs(nom/denom)), 0.0, vϕ]
end
```
Note that we have hardcoded the negative choice of radial component, since we are infalling. Tracing this with

```julia
u_plunging = @SVector[0.0, Gradus.isco(m), deg2rad(90), 0.0]
v_plunging = CircularOrbits.plunging_fourvelocity(m, Gradus.isco(m))

sol = @time tracegeodesics(
    m,
    u_plunging, v_plunging,
    # need a very large integration time, else we never leave the ISCO
    (0.0, 18_000.0);
    abstol=1e-11, reltol=1e-11,
    μ=1.0
)
```
and plotting:

![plunging-four-velocity](/assets/plunging/plunging-four-velocity.svg)

This is good, but the long integration times are perhaps not ideal, since most of the time is just spent tracing the circular orbit.

### Dropping photons from sub-ISCO radii

An alternative method to the above is to _drop_ the photon into the black hole, by moving the orbiting radius a _tiny fraction_ closer. This only requires the $\dot{u}^\phi$ component of the velocity at the ISCO, and introduces $\delta r$ as a free parameter which we can tweak. In code

```julia
u_dropped = @SVector[0.0, Gradus.isco(m) - δr, deg2rad(90), 0.0]
v_dropped = fourvelocity(m, Gradus.isco(m))

sol = @time tracegeodesics(
    m,
    u_dropped, v_dropped,
    (0.0, 500.0);
    abstol=1e-11, reltol=1e-11,
    μ=1.0,
)
```

Tracing for different $\delta r$:

![plunging-dropped-photon](/assets/plunging/plunging-dropped-photon.svg)

Here, the solid line is $\delta r = 10^{-11}$, whereas the dashed line has $\delta r = 0.001$. The latter has significantly lower integration time, as one might expect, but falls nearly identically, modulo some initial rotation $u^\phi$. Plotting the components of velocity as a function of $r$:

![plunging-velocities](/assets/plunging/plunging-velocities.svg)

The above plot requires the keyword argument `closet_approach=1.0001` to be passed to `tracegeodesics`, as default termination is 1% beyond the inner horizon radius. Also, I used $a=0.2$ in the above plot to include non-zero $g_{t\phi}$ terms.

The Cunningham values are from Cunningham (1975)[^2], calculated analytically from first-order equations of motion, and evenly sampled for the sake of illustration.

Note that all three traced curves agree: the one directly from `plunging_fourvelocity`, and the two dropped from a slightly closer radius. Setting $\delta r$ too high causes the paths to diverge from expectation.

## Interpolating plunging velocities

We can attempt to check the deviation from Cunningham's calculations of the above methods -- to do this, we need to interpolate the velocity component calculated by the tracer as a function of radius. This may be done with [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl):

```julia
using Interpolations

# need to sort our values
# drop the first index as it is duplicated
I = sortperm(pl_sol[2,:])[2:end]
r_sol = sol[2,:][I]
ut_sol = sol[5,:][I]

interp_ut = LinearInterpolation(r_sol, ut_sol)
```

Comparing the values with

```julia
using Statistics

isco = Gradus.isco(m)
# range of points to compute error over
# need to shrink the domain a little to avoid 
# stepping outside of the interpolated domain
r_range = range(
    Gradus.inner_radius(m) + 0.1, 
    isco - 0.1, 
    2000
)

# Cunningham analytic calculation
ut_true = map(
    r -> Gradus.AccretionFormulae.RedshiftFunctions.uᵗ(m.M, isco, r, m.a), 
    r_range
)

# find interpolated values
test_t = map(interp_ut, r_range)

# absolute error
err_t = abs.(test_t .- ut_true)
```

The same is done for $u^r$ and $u^\phi$. These components are then used to create a measure of error distance between the true and interpolated velocities, using the $l_2$ norm,

\begin{equation}
    \text{err} = \sqrt{(\Delta u^t)^2 + (\Delta u^r)^2 + (\Delta u^\phi)^2},
\end{equation}

where $\Delta u^\nu$ denotes the difference between the velocity component of the two methods.

As an example, for various offsets, this error by radius is:

![error-scaling-interpolation-single](/assets/plunging/error-scaling-interpolation-single.svg)

The above figure shows the rolling mean, with window size $10^3$ over a total $10^5$ points. Interesting to note is that the error close to the event horizon is left relatively unchanged for different $\delta r$. Instead, perhaps unsurprisingly, the offset primarily impacts the initial error at $r = r_\text{isco}$. To me, this implies there may be some error in the integration over this region that I will need to check, but investigating this is beyond the scope of this post.


To examine how this error changes with offset and integrator tolerance, we will plot the mean and standard deviation over the lifetime of the plunging trajectory. For $\delta r$ in some range between $10^{-11}$ and $10^{-1}$, we find these scale:

![error-scaling-interpolation](/assets/plunging/error-scaling-interpolation.svg)

Both measures scale relatively linearly up until some cutoff $\delta r$, at which point we find a minimum deviation in error, before a plateau. This plateau and the location of the cutoff scales with the `reltol` of the integrator, implying the mean error is indeed a numerical error, and not directly a methodological error. 

Note that the overall best mean error is very approximately $\sqrt{\texttt{reltol}}$, which is useful for anticipating expected errors, and that this best-mean plateau also occurs at approximately $\delta r = \frac{\sqrt{\texttt{reltol}}}{10}$, giving us some way to tune the offset to the integrator tolerances.

We can, however, overall expect the relative error to be no worse than 1 part in 100 using this interpolation method, which should be sufficient for our purposes (though I will note this seems quite high to me). We can test this by recreating a known rendering result, namely the relative redshift of an equatorial accretion disc.

## Plunging redshift

A standard result in relativity is that the energy of a particle with four-momentum $p^\nu$, measured by an observer moving with $\dot{u}^\mu$ is

\begin{equation}
    E = - g_{\mu \nu} p^\mu \dot{u}^\nu > 0.
\end{equation}

If we consider the redshift of the accretion disc as viewed by a distant observer as the ratio of energies measured at the disc and at the observer,

\begin{equation}
    \label{eq:redshift}
    z = \frac{E_\text{obs}}{E_\text{disc}} = \frac{
        g_{\sigma \rho} \, p^\sigma \, \dot{u}_\text{obs}^\rho
    }{
        g_{\mu \nu} \, p^\mu \, \dot{u}_\text{disc}^\nu
    },
\end{equation}

then we can colour the redshift for each pixel in our traced image. Using [Tullio.jl](https://github.com/mcabbott/Tullio.jl) for the Einsum notation, we can easily express a `PointFunction` which will do this for us:

```julia
using Tullio

# need the ISCO radius so we can differentiate plunging and regular orbits
isco = Gradus.isco(m)
# fixed observer velocity
v_obs = @SVector [1.0, 0.0, 0.0, 0.0]
# fixed observer position
# have to put ourselves further away than usual to ensure the metric-flat space
# approximation holds
u_obs = @SVector [0.0, 10_000.0, deg2rad(85), 0.0]

# an interpolation which returns the plunging velocities for us
# (this is part of the standard Gradus.jl exports)
pl_int = interpolate_plunging_velocities(m)

# metric matrix at observer
g_obs = metric(m, u_obs)

redshift = PointFunction((m, gp, mt) -> begin
    let r = gp.u2[2]
        v_disc = if r < isco
            # plunging region
            vtemp = pl_int(r)
            # we have to reverse radial velocity due to backwards tracing convention
            # see https://github.com/astro-group-bristol/Gradus.jl/issues/3
            @SVector [vtemp[1], -vtemp[2], vtemp[3], vtemp[4]]
        else
            # regular circular orbit
            CircularOrbits.fourvelocity(m, r)
        end

        # get metric matrix at position on disc
        g = metric(m, gp.u2)
        @tullio E_disc := - g[i,j] * gp.v2[i] * v_disc[j]
        @tullio E_obs := - g_obs[i,j] * gp.v1[i] * v_obs[j]
        E_obs / E_disc
    end
end)
```
Note that we position the observer at $10^4$ instead of $10^3 \, r_\text{g}$. We do this to ensure the metric is closer to being flat space, and therefore that $E_\text{obs} \approx 1$ for all geodesics. At $10^3 \, r_\text{g}$, the deviation in $E$ is in the order of $1$ part in $100$ for geodesics directed towards the edges of the accretion disc, but at this alternate distance, the deviation is approximately 1 part in 1000. In the literature, the common practice is to assume $E_\text{obs} = 1$, and hard-code this result into the redshift calculation, as is done with common Cunningham implementations, but we will instead use the full definition of redshift, precisely as in Eq. \eqref{eq:redshift}, and suffer the consequences.

We can then render an image in the usual way:
```julia
# set inner radius within the ISCO
d = GeometricThinDisc(
    1.0, 50.0, deg2rad(90)
)

# small utility function
render_redshift(m, rspf) = rendergeodesics(
    m, u_obs, d, 20_000.0;
    pf = rspf ∘ ConstPointFunctions.filter_intersected,
    image_width=750,
    image_height=400,
    fov_factor=16.0,
    effective_infinity=11_000,
    abstol=1e-11, reltol=1e-11,
    verbose=true
)

# comparison metric
m_cmp = BoyerLindquistAD(M=m.M, a=m.a)

# render the redshift images
img_ad = render_redshift(m, redshift)
img_cmp = render_redshift(m_cmp, ConstPointFunctions.redshift)
```
Subtracting these images to compare the difference between our interpolated redshift and the Cunningham analytic values for $a=0$:

![redshift-method-comp-ad-a=-0.0](/assets/plunging/redshift-method-comp-ad-a=0.0.png)

We see the error is in the order of $10^{-4}$, though fortunately this error is relatively robust to integrator tolerances, retaining the order of magnitude even at tolerances of $10^{-8}$.

For first-order method comparison, where `m_cmp = BoyerLindquistFO(M=m.M, a=m.a)`, the situation is a little worse, but not dire:

![redshift-method-comp-fo-a=0.0](/assets/plunging/redshift-method-comp-fo-a=0.0.png)

This is all pretty good news. Sort of. A difference of $10^{-4}$ is still about 5 orders of magnitude larger than I would desire. The main difference is clearly within the ISCO, and judging by the features in the noise, is caused by some geometrical anomaly. I will investigate this more at a later date.

I have a strong suspicion the differences in the first-order method outside of the ISCO, which look considerably more like noise, are correlated with numerical differences in how the $\beta$ impact parameter is mapped, but I would need to test to confirm -- which, again, is beyond the scope of this post.

But, with a proof of concept, we just have one missing ingredient before we can calculate this redshift for other metrics!

## Determining the ISCO

If we want this interpolation method to work on general static, axis-symmetric methods, the only thing left to define is a function to find the ISCO.

The ISCO is the radius corresponding to

\begin{equation}
    \frac{\partial}{\partial r} E = 0,
\end{equation}

which may be found relatively easily via Eq. \eqref{eq:stable_energy}. The strategy here is to use automatic differentiation to compute the gradient, and use a root finding tool to optimize finding the minima. With the help of [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) and [Roots.jl](https://github.com/JuliaMath/Roots.jl), this is pretty easy and fast:

```julia
using ForwardDiff
using Roots

function isco(m, lower_bound upper_bound)
    dE(r) = ForwardDiff.derivative(
        x -> CircularOrbits.energy(m, x), 
        r
    )
    d2E(r) = ForwardDiff.derivative(dE, r)

    find_zero((dE, d2E), (lower_bound, upper_bound))
end
```

We have to use the bounded solvers, since our domain is undefined below a certain $r$, and throws a `DomainError`. The upper bound may be quite confidently set to something large, e.g. $100 \, r_\text{g}$, however there is generally no safe lower bound we can use that works for all $a$. Instead, an approach is determine the lower bound coarsely as a radius for which $E / \mu > 1$. This may be programmatically found with a simple brute-force:

```julia
function find_lower_isco_bound(m; upper_bound = 100.0, step=-0.005)
    # iterate in reverse with a negative step
    for r in upper_bound:step:1.0
        if CircularOrbits.energy(m, r) > 1.0
            return r
        end
    end
    # for type stability
    return 0.0
end
```

This function, despite its crudeness, performs fine, and seems to be able to find a suitable solution for all currently implemented metrics:
```julia
julia> @btime find_lower_isco_bound($m)
  92.847 μs (0 allocations: 0 bytes)
```

Adding a dispatch to make our `isco` API a little friendlier:
```julia
isco(m; upper_bound = 100.0) = isco(m, find_lower_isco_bound(m; upper_bound), upper_bound)
```

Testing this for a few different metric setups:

![iscos-various-metrics](/assets/plunging/iscos-various-metrics.svg)

Seems to work well! We should now be able to visualise redshift patterns for any of our static, axis-symmetric metrics.

## Redshifts in a non-GR spacetime

Just for a laugh, let's see what the redshift of an equitorial accretion disc for an Einstein-Maxwell-Dilaton-Axion singularity would look like. The metric components for this I have taken from García et al. (1995)[^3], allowing reuse of all of our above code, with just one alteration:

```julia
m = DilatonAxionAD(M=1.0, a=0.998, β=1.0, b=1.0)
```

This metric is a little fiddly to render at the moment, because I probably have some mathematical typos in it, and, since the Schwarzschild equivalent radius is currently not implemented, the cutoff on the inner horizon is underestimated. To avoid uneccessarily long integrations, we can just fix the minimum timestep to something an order of magnitude smaller than our tolerances:

```julia
img = @time rendergeodesics(
    m, u_obs, d, 2000.0;
    pf = redshift ∘ ConstPointFunctions.filter_intersected,
    image_width=1400,
    image_height=800,
    fov_factor=18.0,
    abstol=1e-10, reltol=1e-10,
    verbose=false,
    closest_approach=1.01,
    dtmin=1e-9
)
```

We set `verbose=false` to avoid excessive warnings concerning iteration counts, and produce the following render: 

![dilaton-redshift](/assets/plunging/dilaton-redshift.png)

Gorgeous! But no-doubt highly unphysical, and probably litered with mathematical typos. There are also clearly numerical anomalies and instabilities visible in the render, which appears to be related to close-to-ISCO redshifts -- plently of bug hunting and optimizing to be done.

At some point I will endevour to properly follow the EMDA derivation and motivations, and maybe there's a blog post in that, but this brief tourist's curiousity outside of GR is sufficient for the spectacle, and shows the method of calculating redshift does indeed "work" for static, axis-symmetric metrics.

## References

[^1]: T Johannsen (2013). _Regular black hole metric with three constants of motion_. PHYS REV D **88**:4.

[^2]: C T Cunningham (1975). _The effects of redshifts and focusing on the spectrum of an accretion disk around a Kerr black hole_. APJ, **202**:788-802.

[^3]: A García, D Galtsov, and O Kechkin (1995). _Class of stationary axisymmetric solutions of the Einstein-Maxwell-Dilaton-Axion field equations_. PHYS REV L **74**:1276.
