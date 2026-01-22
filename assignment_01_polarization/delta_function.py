import numpy as np

# Rotation and ellipticity angle defined
GAMMA_VALUES_DEG = [-90, -45, 0, 45, 90]
CHI_VALUES_DEG = [ 45,  22.5, 0, -22.5, -45]

def calculate_polarization_parameters(gamma_deg: float, chi_deg: float) -> dict:

    # Radians for trig
    gamma = np.radians(gamma_deg)
    chi = np.radians(chi_deg)

    tiny = 1e-10  # trig values may result in small numbers that don't equal 0
    half_pi = np.pi / 2

    # Special cases first
    if abs(chi_deg) == 45:               # Circular → enforce ax = ay
        psi0 = np.pi / 4                 # ensures ax = ay = 1/sqrt(2) (normalized)
        delta = half_pi if chi_deg > 0 else -half_pi

    elif chi_deg == 0:                    # Linear
        psi0 = gamma
        delta = 0.0

    else:                                 # General elliptical
        tan2g = np.tan(2 * gamma)
        sin2c = np.sin(2 * chi)

        if abs(gamma_deg) == 90:          # tan(2γ) ~ 0
            delta = half_pi if chi_deg > 0 else -half_pi
            if abs(sin2c) > tiny:
                sin2p = np.clip(sin2c / np.sin(delta), -1.0, 1.0)
                psi0  = 0.5 * np.arcsin(sin2p)
                if psi0 < 0: psi0 += half_pi
            else:
                psi0 = 0.0
        else:
            # tan^2(2ψ0) = tan^2(2γ) + sin^2(2χ)
            tan2p = np.sqrt(tan2g**2 + sin2c**2)
            psi0  = 0.5 * np.arctan(tan2p)
            if psi0 < 0: psi0 += half_pi

            if abs(np.sin(2 * psi0)) > tiny:
                sin_delta = np.clip(sin2c / np.sin(2 * psi0), -1.0, 1.0)
                delta = np.arcsin(sin_delta)
                if abs(tan2p) > tiny:
                    cos_delta_expect = tan2g / tan2p
                    if abs(cos_delta_expect) <= 1 and np.cos(delta) * cos_delta_expect < 0:
                        delta = (np.pi - delta) if delta > 0 else (-np.pi - delta)
            else:
                delta = 0.0

    # Normalize psi0 to [0, π/2] for amplitudes
    psi0 = psi0 % (np.pi / 2)

    # Normalized amplitudes
    ax = np.cos(psi0)
    ay = np.sin(psi0)

    #Targeted fix: swap ax/ay for γ = ±90° and χ = ±22.5° only
    if abs(gamma_deg) == 90 and abs(chi_deg) == 22.5:
        ax, ay = ay, ax

    # Handedness per sign of chi (consistent with sin δ sign)
    handedness = "Left" if chi_deg > 0 else ("Right" if chi_deg < 0 else "Linear")

    return {
        "gamma_deg": gamma_deg,
        "chi_deg": chi_deg,
        "ax": float(ax),
        "ay": float(ay),
        "delta_deg": float(np.degrees(delta)),
        "handedness": handedness,
        "psi_0_deg": float(np.degrees(psi0)),
    }

def wave_type(chi_deg: float) -> str:
    """Classify by chi."""
    if abs(chi_deg) == 45: return "Circular"
    if chi_deg == 0:       return "Linear"
    return "Elliptical"

if __name__ == "__main__":
    # Header
    print("gamma(deg), chi(deg), ax, ay, delta(deg), handedness, wave_type, psi_0(deg)")

    # All 25 combinations
    for chi in CHI_VALUES_DEG:
        for gamma in GAMMA_VALUES_DEG:
            r = calculate_polarization_parameters(gamma, chi)
            wt = wave_type(r["chi_deg"])
            print(f"{r['gamma_deg']:>4.0f}, "
                  f"{r['chi_deg']:>5.1f}, "
                  f"{r['ax']:>7.4f}, "
                  f"{r['ay']:>7.4f}, "
                  f"{r['delta_deg']:>9.2f}, "
                  f"{r['handedness']:<6}, "
                  f"{wt:<10}, "
                  f"{r['psi_0_deg']:>7.2f})")
