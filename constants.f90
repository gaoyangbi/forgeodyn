!""" Physical constants needed for algo of the diffusion among others. """
module constants
    !# Distance between the core-mantle boundary (CMB) and the layer just below the radial grid of the Coupled-Earth model
    real(kind=8):: delta = 2.7033
    !# Effective magnetic diffusion coefficient of the Coupled-Earth model
    real(kind=8):: eta_mag = 36.577129


    !# Radius of core
    real(kind=8):: r_core = 3485.0
    !# Radius of Earth
    real(kind=8):: r_earth = 6371.2
    !# Radius usde in covobs decomposition
    real(kind=8):: r_covobs = 6371.2

    !# Julian Date of Jan 1, 2000 12h
    real(kind=8):: jd2000 = 2451545.0
end module   
