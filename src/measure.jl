# The temperature ``T`` is measured by computing the average kinetic energy per degree of freedom.
# ```math
# k_B T = \frac{\langle 2 \mathcal{K} \rangle}{f}.
# ```

for (step, it) in enumerate(
            molecule_dynamics(; lattice_pos, box, temperature, rc, Δt)
        )
    
end