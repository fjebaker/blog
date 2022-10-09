## Frames

## The Gramm-Schmidt process

The Gramm-Schmidt process gives us a method for calculating orthonormal vectors in any Euclidean space with an inner-product. This may naturally be extended:

@@theorem
**Gramm-Schmidt:**

It is always possible to construct an orthonormal frame on any Riemannian manifold, i.e. a space equipped with a smoothly varying metric, and a standard inner-product operation on the tangent vector space.
@@

The proof lies in the process: given a manifold $\mathcal{M} \subset \mathbb{R}^n$, and a set of linearly independent vectors tangent to $x \in \mathcal{M}$, denoted $S = \left\{ v^{\mu}_{( 1 )}, \ldots, v^{\mu}_{( n )} \right\}$, it is possible to construct an orthonormal frame for the subspace spanned by $S$, denoted $\left\{ e^{\mu}_{( 1 )}, \ldots, e^{\mu}_{( n )} \right\}$, by iteratively subtracting the projections of the vectors.

Defining the projection operator as 

\begin{equation}
    \text{P}_{w^{\mu}} \left(
        v^{\mu}
    \right) = \frac{
        \iprod{v^{\mu}}{w^{\mu}}
    }{
        \iprod{w^{\mu}}{w^{\mu}}
    } w^{\mu},
\end{equation}

the Gramm-Schmidt procedure is then

\begin{equation}
    \label{eq:gs}
    w^{\mu}_{( i )} = v^{\mu}_{( i )} - \sum_{j=1}^{i-1} \text{P}_{w^{\mu}_{(j)}} \left( v^{\mu}_{(i)} \right)
\end{equation}

normalized through $e^{\mu}_{( i )} = \frac{w^{\mu}_{( i )}}{ \lVert w^{\mu}_{( i )} \rVert}$. 

The Riemannian inner-product requires the metric, and so $\iprod{w^{\mu}}{v^{\mu}}$ is perhaps more conveniently written directly as $g_{\mu \nu} w^{\mu} v^{\nu}$. A Julia implementation of the Gramm-Schmidt process may then be

```julia
using Tullio
using StaticArrays

# inner product
iprod(g, v, w) = @tullio p := g[i,j] * v[i] * w[j]
# projection operator
P(g, v, w) = (iprod(g, v, w) / iprod(g, w, w)) .* w

function grammschmidt(g::SMatrix{N,N,T}, S::NTuple{N,V}) where {N,T,V}
    M = zeros(MMatrix{N,N,T})
    for (i, v) in enumerate(S)
        w = if i == 1
            v
        else
            v - sum(j -> P(g, v, M[:,j]), 1:i-1; init=zero(SVector{N,T}))
        end
        @inbounds M[:,i] .= w ./ √abs(iprod(g, w, w))
    end
    SMatrix{N,N,T}(M)
end
```

Note: Julia already has a function which calculates the orthonormal space for a given matrix, similar to this Gramm-Schmidt procedure, however this assumes Euclidean inner-product -- see [`LinearAlgebra.nullspace`](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#LinearAlgebra.nullspace) for more.

The above implementation will return a matrix where the columns are the orthonormal frame spanning the subspace given by $S$.

## Tetrads of the Kerr metric

In order to identify the direction of proper time in the tetrad frame, we may follow the prescription of Krolik et al. (2005), and set $e^\mu_{(t)} = \dot{u}^\mu$, where $\dot{u}^\mu$ is the four-velocity of a geodesic at $u^\mu$. 



Infact, the result of analytically applying the Gramm-Schmidt method is given in Appendix A of Beckwith, Hawley and Krolik (2008). We can do the same, using Symbolics.jl

