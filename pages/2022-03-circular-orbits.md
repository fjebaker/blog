@def title = "Circular orbits in static, axis-symmetric spacetimes"
@def description = "TODO"
@def image = ""
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

## Circular orbit constraints in the Kerr metric

It is a known result that circular orbits exist in the equatorial plane of static axis-symmetric spacetimes, i.e. $\theta = \frac{\pi}{2}$ and $\dot{\theta} = 0$, shown by Carter (1968)[^3]. These circular orbits are further constrained by motion in the radial components,

\begin{equation}
    \label{eq:ur_dotur}
    \dot{u}^r = 0, \quad \text{and} \quad \frac{\partial}{\partial r}\dot{u}^r = 0.
\end{equation}

We can use the invariance of four velocities to express these constraints further,

\begin{equation}
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

If we rephrase the above constraint in terms of potentials, we find the expression may be written much more simply as

\begin{equation}
    \tensor{\dot{u}}{r}{} = \frac{\pm 1}{ \e^{\mu_1} } \sqrt{ \mathcal{R} } = 0.
\end{equation}

Our constraints are now $\mathcal{R} = 0$ and $\frac{\partial}{ \partial r}\mathcal{R} = 0$, the latter of which is directly

\begin{equation}
\label{eq:sage-dont-like}
    \left( E + \omega L_z \right)^2 \frac{\partial}{\partial r} \e^{-2\nu}
    + 2 \left( E + \omega L_z \right) \e^{-2\nu} L_z \frac{\partial}{\partial r} \omega
    - L_z^2 \frac{\partial}{\partial r} \e^{-2 \Psi}
    = 0.
\end{equation}

Implicitly it is assumed that $E$ and $L_z$ are constant with respect to $r$ over an orbit. 

Here we will specialise to the Kerr metric, with subsequent manipulations performed using SageMath, and the occasional small bit of hand manipulation. Substituting in the metric components from eq. \eqref{eq:metric-comp}, we obtain

\begin{equation}
    \label{eq:constraint-1}
    \mathcal{R} = 2 M r {\left(\mu^{2} r^{2} + {\left(E a - L_z\right)}^{2}\right)} + {\left({\left(E^{2} - \mu^{2}\right)} {\left(a^{2} + r^{2}\right)} - L_z^{2}\right)} \Sigma = 0
\end{equation}
\begin{equation}
    \label{eq:constraint-2}
    \frac{\partial}{\partial r} \mathcal{R} = 6 \, M \mu^{2} r^{2} + 2 \, {\left(E a - L_z\right)}^{2} M + 2 \, {\left(E^{2} - \mu^{2}\right)} r \Sigma + {\left({\left(E^{2} - \mu^{2}\right)} {\left(a^{2} + r^{2}\right)} - L_z^{2}\right)} \frac{\partial}{\partial r}\Sigma = 0.
\end{equation}

It is worth noting that $\mu$ has reappeared in the derivative term. This is due to the complexity of eq. \eqref{eq:sage-dont-like} when used in SageMath. Though they both solve identically, this version if more friendly for a blog post.

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

Our constrains in eq. \eqref{eq:ur_dotur} are now eq. \eqref{eq:constraint_1} and \eqref{eq:constraint_2} respectively, which may be solved simultaneously to find $E$ and $L_z$ satisfying circular orbits. In SageMath, this is fiddly, but generally quite straight forward to solve. We find there are four solutions in total, representing the intersection of the positive and negative energy regimes with prograde and retrograde orbits.

\begin{equation}
\label{eq:E2mu2}
\frac{E^{2}}{\mu^{2}} =
    \frac{
        {\left(
           M^2 \left( 5 \, a^{2} + 16 \, r^{2} \right) - 7 \, M r^{3} + r^{4} 
            - 3 \, {\left(4M^{2} +a^{2}\right)} M r 
            \pm 2 a (Mr)^\frac{3}{2} {\left( - 1 + \frac{2 M}{r} - \frac{a^{2}}{r^2} \right)} 
        \right)}
    }{
        {\left( r^3 - 6 \, M r^{2} + 9 \, M^{2} r - 4 \, M a^{2} \right)} r
    },
\end{equation}

where $\pm$ denotes prograde and retrograde respectively, and

\begin{equation}
\label{eq:Lz2mu2}
    \frac{L_z^{2}}{\mu^{2}} = 
        \left(
            \frac{E^{2}}{\mu^2} 
            - 1
        \right)
        \left(
            a^{2} + 3r^2
        \right)
        + 4 \, M r.
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

which are substantially different, in form, from eq. \eqref{eq:E2mu2} and eq. \eqref{eq:Lz2mu2}. We can try to find the degree of similarity by plotting the equations with $r > 0.1$ to avoid the central singularity:

![energy-diff](/assets/angmom-energy.svg)


The plots for energy (left panel) are in good agreement up until the energy becomes negative in the Bardeen case (such energies are allowed, e.g. in the [Penrose process](https://en.wikipedia.org/wiki/Penrose_process)). The SageMath derivation fails in this regime, as it is the square root of a squared quantity, which loses sign information. Indeed, we do have a negative energy solution from SageMath, but we would require _a priori_ knowledge of which solution to apply in which domain.

The consequence of the energy sign is even more pronounced in the angular momentum (right panel), as the angular momentum is calculated from the energy. For the $a=0.998$ case, the plots are still in good agreement, but the divergence for $a=0.0$ is pronounced, and worse still for $a=-0.988$ (not depicted, as the diagram became difficult to annotate). Some of the features of the SageMath $L_z$ curves do not disappear when using a sign-concious approach, such as the sharp cusp at $r\sim2$.

Fortunately, the differences in the two solutions is for our purposes pedagogical, since the divergence lies within the _marginally bound_ orbits, at which radii circular orbits are prohibited.


## Special orbit boundaries

There are many special boundaries in the Kerr metric, including the event-horizon and the ergosphere, setting the boundary of the static limit, at which the frame dragging effect appears which originally motivated the LNRF. For purely circular orbits, we consider two of these boundaries and the respective possible orbits within them, as they are arguably most relevant to our study -- these are the _innermost stable circular orbit_ (ISCO) and the _marginally bound orbits_.

### ISCO

### Marginally bound orbits


## Plunging region

No longer have $\dot{r} = 0$.


## References

[^1]: J M Bardeen, W H Press, S A Teukolsky (1972). _Rotating black holes: locally nonrotating frames, energy extraction, and scalar synchrotron radiation_. APJ, **178**:347-169.

[^2]: C T Cunningham (1975). _The effects of redshifts and focusing on the spectrum of an accretion disk around a Kerr black hole_. APJ, **202**:788-802.

[^3]: B Carter (1968). _Global structure of the Kerr family of gravitational fields_. PHYS REV **174**:5.