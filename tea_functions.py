"""
tea_functions.py — Pure functions extracted from bioprocess_tea_calculator.ipynb
All functions are side-effect-free; no print/plot calls except inside verbose guards.
"""
import re
import math
import numpy as np

# ── Atomic weights & molecular weights ──────────────────────────────────────
ATOMIC_WEIGHTS = {'C': 12.011, 'H': 1.008, 'O': 15.999, 'N': 14.007, 'S': 32.06}
MW_GLUCOSE = 180.156
MW_NH3     = 17.031
MW_O2      = 31.998
MW_H2SO4   = 98.079
MW_MGSO4   = 120.368
MW_BIOMASS = 95.37   # empirical formula C3.85H6.69O1.78N (Battley 1987)

# ── Biomass stoichiometry constants ─────────────────────────────────────────
BIOMASS_YIELD_COEFF = 0.8 * MW_BIOMASS / (0.84 * MW_GLUCOSE)  # gCDW / g glucose
BIOMASS_O2_COEFF    = MW_BIOMASS / (1.212 * MW_O2)             # g biomass / g O2

# ── Fermentation model constants ─────────────────────────────────────────────
INOCULUM_FRACTION   = 0.01    # starting biomass = 1% of final
O2_COOLING_COEFF    = 0.460   # kJ per mmol O2 consumed (Doran 1995)
O2_SATURATION       = 0.2     # mmol O2/L at air saturation
DO_SETPOINT         = 0.25    # dissolved O2 maintained at 25% of saturation

# ── Plant constants ───────────────────────────────────────────────────────────
WORKING_VOL_RATIO   = 0.85    # working volume / total tank volume
ASPECT_RATIO        = 3.0     # fermenter height:diameter ratio
O2_MOLES_PER_M3_AIR = 9.375  # mol O2 / m³ air at STP
DO_UTILISATION      = 0.75    # fraction of O2 in sparged air that is consumed

# ── DSP presets ───────────────────────────────────────────────────────────────
DSP_PRESETS = {
    'small_molecule':      {'yield': 0.90, 'opex_frac': 0.20, 'capex_frac': 0.20},
    'industrial_enzyme':   {'yield': 0.70, 'opex_frac': 0.30, 'capex_frac': 0.35},
    'therapeutic_protein': {'yield': 0.50, 'opex_frac': 0.60, 'capex_frac': 0.50},
}

# ── Organism presets ─────────────────────────────────────────────────────────
ORGANISM_PRESETS = {
    'Generic (model default)': {
        'biomass_yield_coeff': BIOMASS_YIELD_COEFF,   # ~0.504 gCDW/g glucose
        'media_cost':          0.40,
        'note': 'Battley 1987 empirical average (C₃.₈₅H₆.₆₉O₁.₇₈N). '
                'Reasonable for E. coli, B. subtilis on glucose.',
    },
    'E. coli (aerobic)': {
        'biomass_yield_coeff': 0.48,
        'media_cost':          0.25,
        'note': 'Simple mineral salts medium; aerobic growth on glucose.',
    },
    'S. cerevisiae (yeast)': {
        'biomass_yield_coeff': 0.45,
        'media_cost':          0.40,
        'note': 'Aerobic; mineral medium. Lower yield than bacteria on glucose.',
    },
    'B. subtilis': {
        'biomass_yield_coeff': 0.50,
        'media_cost':          0.30,
        'note': 'Aerobic; simple mineral medium. Similar to E. coli.',
    },
    'Mammalian (CHO-like)': {
        'biomass_yield_coeff': 0.20,
        'media_cost':          5.00,
        'note': 'Very rough approximation only. Complex medium required; '
                'high media cost. Model assumes aerobic single-substrate '
                'fermentation — CHO bioreactors are substantially more complex.',
    },
}

# ── Financial defaults (used as defaults in app.py) ──────────────────────────
RAMP_FRACTIONS    = [0.50, 0.75, 1.00]
CAPEX_YR1_FRAC    = 0.70
CAPEX_YR2_FRAC    = 0.30
ONGOING_CAPEX_FRAC = 0.10
DEPRECIATION_YR   = 10


# ════════════════════════════════════════════════════════════════════════════════
# SECTION 1 — CHEMISTRY
# ════════════════════════════════════════════════════════════════════════════════

def parse_formula(formula_str):
    """Parse a molecular formula string into atom counts {C, H, O, N, S}."""
    atoms = {'C': 0, 'H': 0, 'O': 0, 'N': 0, 'S': 0}
    pattern = r'([A-Z][a-z]?)(\d*)'
    matches = re.findall(pattern, formula_str)
    for atom, count in matches:
        if not atom:
            continue
        if atom in atoms:
            atoms[atom] = int(count) if count else 1
        else:
            raise ValueError(
                f"Atom '{atom}' is not supported. "
                f"This calculator handles C, H, O, N, S only."
            )
    if atoms['C'] == 0:
        raise ValueError("Product must contain carbon (C).")
    return atoms


def calculate_molecular_weight(a, b, c, d, e=0):
    """Molecular weight of CaHbOcNdSe in g/mol."""
    return (ATOMIC_WEIGHTS['C'] * a + ATOMIC_WEIGHTS['H'] * b +
            ATOMIC_WEIGHTS['O'] * c + ATOMIC_WEIGHTS['N'] * d +
            ATOMIC_WEIGHTS['S'] * e)


def calculate_h2_co2_ratio(a, b, c, d, e=0):
    """
    H2:CO2 ratio for product CaHbOcNdSe.
    >2 = more reduced than glucose; =2 = neutral; <2 = more oxidised.
    +3*(e/a) accounts for the redox cost of sulfate reduction to organic S
    (8 electrons per S, but the coefficient is 3 after accounting for O/H
    contributed by H2SO4 to the atom balance).
    From Equation 3, Lynch 2021 (extended for S).
    """
    return 0.5*(b/a) - 1.0*(c/a) - 1.5*(d/a) + 3.0*(e/a) + 2.0


def balance_equation(a, b, c, d, e=0):
    """
    Balance the stoichiometric equation for glucose -> product (per 1 mol glucose).
    H2SO4 is the stoichiometric S source (same atom balance as MgSO4 or ammonium sulfate).
    Returns dict with keys: case, description, ratio, glucose, NH3, H2SO4, O2,
    product, H2O, CO2.
    From Supplemental Materials Section 1, Lynch 2021 (extended for S).
    """
    ratio = calculate_h2_co2_ratio(a, b, c, d, e)
    tol   = 1e-9

    if abs(ratio - 2.0) < tol:
        X = 6 / a
        W = X * d
        V = X * e
        Z = (12 + 3*W + 2*V - X*b) / 2
        return {'case': 1, 'description': 'Neutral relative to glucose',
                'ratio': ratio, 'glucose': 1.0,
                'NH3': W, 'H2SO4': V, 'O2': 0.0, 'product': X, 'H2O': Z, 'CO2': 0.0}

    elif ratio > 2.0:
        X = 6 / (a + 0.25*b - 0.5*c - 0.75*d + 1.5*e)
        W = X * d
        V = X * e
        Q = 6 - X * a
        Z = (12 + 3*W + 2*V - X*b) / 2
        return {'case': 2, 'description': 'More reduced than glucose',
                'ratio': ratio, 'glucose': 1.0,
                'NH3': W, 'H2SO4': V, 'O2': 0.0, 'product': X, 'H2O': Z, 'CO2': Q}

    else:
        Y = 6 / a
        W = Y * d
        V = Y * e
        Z = (12 + 3*W + 2*V - Y*b) / 2
        O2 = (Y*c + Z - 6 - 4*V) / 2
        return {'case': 3, 'description': 'More oxidised than glucose',
                'ratio': ratio, 'glucose': 1.0,
                'NH3': W, 'H2SO4': V, 'O2': O2, 'product': Y, 'H2O': Z, 'CO2': 0.0}


def calculate_theoretical_yield(a, b, c, d, e=0):
    """Theoretical max yield in g product per g glucose."""
    MW_product = calculate_molecular_weight(a, b, c, d, e)
    eq = balance_equation(a, b, c, d, e)
    return eq['product'] * MW_product / MW_GLUCOSE


def calculate_yield_coefficients(a, b, c, d, e=0):
    """
    Yield coefficients for all consumed raw materials (g product / g feedstock).
    Returns dict with keys: glucose, NH3 (None if no N), H2SO4 (None if no S),
    O2 (None if Case 1 or 2).
    """
    MW_product = calculate_molecular_weight(a, b, c, d, e)
    eq = balance_equation(a, b, c, d, e)
    g_product = eq['product'] * MW_product
    return {
        'glucose': g_product / MW_GLUCOSE,
        'NH3':     g_product / (eq['NH3']   * MW_NH3)   if eq['NH3']   > 1e-9 else None,
        'H2SO4':   g_product / (eq['H2SO4'] * MW_H2SO4) if eq['H2SO4'] > 1e-9 else None,
        'O2':      g_product / (eq['O2']    * MW_O2)    if eq['O2']    > 1e-9 else None,
    }


def calc_txl_glucose_overhead(MW_product, avg_residue_mw=110.0, atp_per_residue=5.0):
    """
    Extra glucose (g) consumed per gram of protein for the ATP cost of tx/tl.
    Assumes 1 glucose -> 32 ATP. Default 5 ATP-equiv/residue
    (4 translation + 1 mRNA amortised).
    """
    return atp_per_residue * MW_GLUCOSE / (avg_residue_mw * 32.0)


def run_chemistry(formula_str, is_protein=False,
                  avg_residue_mw=110.0, atp_per_residue=5.0, verbose=False):
    """
    Run full chemistry section for a given molecular formula.

    Returns dict with keys: formula, atoms, MW, equation, theoretical_yield,
    txl_glucose_g_per_g, yield_coeffs, is_protein.
    Raises ValueError on unsupported atoms or invalid formula.
    """
    atoms = parse_formula(formula_str)
    a, b, c, d = atoms['C'], atoms['H'], atoms['O'], atoms['N']
    e = atoms.get('S', 0)

    MW  = calculate_molecular_weight(a, b, c, d, e)
    eq  = balance_equation(a, b, c, d, e)
    thy = calculate_theoretical_yield(a, b, c, d, e)
    yc  = calculate_yield_coefficients(a, b, c, d, e)
    txl_overhead = calc_txl_glucose_overhead(MW, avg_residue_mw, atp_per_residue) if is_protein else 0.0

    if verbose:
        print(f"{'='*60}")
        print(f"  CHEMISTRY: {formula_str}")
        print(f"{'='*60}\n")
        print(f"  Molecular weight : {MW:.2f} g/mol")
        print(f"  H2:CO2 ratio     : {eq['ratio']:.4f}")
        print(f"  Redox class      : {eq['description']} (Case {eq['case']})\n")
        print(f"  Theoretical yield : {thy:.4f} g {formula_str} / g glucose")
        if is_protein:
            eff_yield = 1.0 / (1.0/thy + txl_overhead)
            print(f"  Tx/tl overhead    : {txl_overhead:.4f} g glucose / g protein")
            print(f"  Effective yield   : {eff_yield:.4f} g {formula_str} / g glucose")

    return {
        'formula': formula_str, 'atoms': atoms, 'MW': MW,
        'equation': eq, 'theoretical_yield': thy,
        'txl_glucose_g_per_g': txl_overhead,
        'yield_coeffs': yc, 'is_protein': is_protein,
    }


# ════════════════════════════════════════════════════════════════════════════════
# SECTION 2 — FERMENTATION MODEL
# ════════════════════════════════════════════════════════════════════════════════

def run_fermentation_model(titer, rate, yield_fraction, chemistry,
                            biomass_yield_coeff=None, biomass_o2_coeff=None):
    """
    Run the complete fermentation model.

    Parameters
    ----------
    titer               : float  Final product concentration (g/L).
    rate                : float  Average volumetric production rate (g/L/hr).
    yield_fraction      : float  Fraction of theoretical yield achieved (0–0.99).
    chemistry           : dict   Output from run_chemistry().
    biomass_yield_coeff : float  gCDW per g glucose. Defaults to BIOMASS_YIELD_COEFF
                                 (~0.504). Use ORGANISM_PRESETS to get organism values.
    biomass_o2_coeff    : float  g biomass per g O2. Defaults to BIOMASS_O2_COEFF.

    Returns
    -------
    dict with glucose partitioning, biomass, kinetics, oxygen, cooling, and
    time-course arrays (t_points, biomass_curve, product_curve).
    """
    _biomass_yield = biomass_yield_coeff if biomass_yield_coeff is not None else BIOMASS_YIELD_COEFF
    _biomass_o2    = biomass_o2_coeff    if biomass_o2_coeff    is not None else BIOMASS_O2_COEFF

    theoretical_yield = chemistry['theoretical_yield']
    eq                = chemistry['equation']
    txl_overhead      = chemistry.get('txl_glucose_g_per_g', 0.0)

    sugar_to_product = titer / theoretical_yield
    stoich_total     = sugar_to_product / yield_fraction
    sugar_for_txl    = titer * txl_overhead
    total_sugar      = stoich_total + sugar_for_txl
    sugar_to_biomass = stoich_total - sugar_to_product

    final_biomass    = sugar_to_biomass * _biomass_yield
    starting_biomass = INOCULUM_FRACTION * final_biomass

    ferm_time            = titer / rate
    A                    = (final_biomass - starting_biomass) / starting_biomass
    logistic_growth_rate = -math.log(0.01 / A) / ferm_time
    product_to_cell_ratio = titer / (final_biomass - starting_biomass)
    logistic_prod_rate   = product_to_cell_ratio * logistic_growth_rate
    specific_rate        = rate / final_biomass

    O2_for_biomass = final_biomass / _biomass_o2 * (1000.0 / MW_O2)

    if eq['case'] == 3 and eq['O2'] > 1e-9:
        product_O2_yield_coeff = chemistry['yield_coeffs']['O2']
        O2_for_product = (titer / product_O2_yield_coeff) * (1000.0 / MW_O2)
    else:
        O2_for_product = 0.0

    theoretical_biomass_yield_100pct = MW_BIOMASS / (0.84 * MW_GLUCOSE)
    glucose_at_100pct = final_biomass / theoretical_biomass_yield_100pct
    waste_glucose     = sugar_to_biomass - glucose_at_100pct
    O2_for_waste      = waste_glucose * (6.0 * MW_O2 / MW_GLUCOSE) * (1000.0 / MW_O2)
    O2_for_txl        = sugar_for_txl * (6.0 * MW_O2 / MW_GLUCOSE) * (1000.0 / MW_O2)

    cumulative_O2 = O2_for_biomass + O2_for_product + O2_for_waste + O2_for_txl

    max_OTR_biomass = ((1.0 / _biomass_o2) * (1000.0 / MW_O2)
                       * (logistic_growth_rate / 4.0) * final_biomass)

    if eq['case'] == 3 and O2_for_product > 0:
        max_product_rate = (product_to_cell_ratio * logistic_growth_rate
                            * final_biomass / 4.0)
        max_OTR_product  = ((1.0 / chemistry['yield_coeffs']['O2'])
                            * (1000.0 / MW_O2) * max_product_rate)
    else:
        max_OTR_product = 0.0

    max_OTR         = max_OTR_biomass + max_OTR_product
    max_O2_gradient = O2_SATURATION * (1.0 - DO_SETPOINT)
    max_kla         = (max_OTR / max_O2_gradient) / 3600.0

    max_cooling_rate = O2_COOLING_COEFF * max_OTR
    ave_cooling_rate = O2_COOLING_COEFF * (cumulative_O2 / ferm_time)

    t_points      = np.linspace(0, ferm_time, 300)
    K             = final_biomass
    biomass_curve = K / (1.0 + A * np.exp(-logistic_growth_rate * t_points))
    product_curve = product_to_cell_ratio * (biomass_curve - starting_biomass)

    return {
        'titer': titer, 'rate': rate, 'yield_fraction': yield_fraction,
        'theoretical_yield': theoretical_yield,
        'sugar_to_product': sugar_to_product,
        'total_sugar': total_sugar,
        'sugar_to_biomass': sugar_to_biomass,
        'final_biomass': final_biomass,
        'starting_biomass': starting_biomass,
        'ferm_time': ferm_time,
        'A': A,
        'logistic_growth_rate': logistic_growth_rate,
        'product_to_cell_ratio': product_to_cell_ratio,
        'logistic_prod_rate': logistic_prod_rate,
        'specific_rate': specific_rate,
        'O2_for_biomass': O2_for_biomass,
        'O2_for_product': O2_for_product,
        'O2_for_waste': O2_for_waste,
        'cumulative_O2': cumulative_O2,
        'max_OTR': max_OTR,
        'max_kla': max_kla,
        'max_cooling_rate': max_cooling_rate,
        'ave_cooling_rate': ave_cooling_rate,
        'overall_yield': titer / total_sugar,
        'sugar_for_txl': sugar_for_txl,
        'O2_for_txl': O2_for_txl,
        't_points': t_points,
        'biomass_curve': biomass_curve,
        'product_curve': product_curve,
    }


# ════════════════════════════════════════════════════════════════════════════════
# SECTION 3 — OPEX
# ════════════════════════════════════════════════════════════════════════════════

def calculate_plant_logistics(capacity_kta, annual_uptime, batches_on_spec,
                               tank_volume_L, turnaround_time, fermentation):
    """
    Calculate annual fermentation volume and number of tanks required.

    Parameters
    ----------
    capacity_kta     : float  target annual production in kilotonnes
    annual_uptime    : float  fraction of year plant operates (e.g. 0.90)
    batches_on_spec  : float  fraction of batches meeting quality spec
    tank_volume_L    : float  main fermenter total volume in litres
    turnaround_time  : float  hours between batches (clean + refill)
    fermentation     : dict   output from run_fermentation_model()

    Returns
    -------
    dict of plant logistics quantities
    """
    annual_uptime_hr        = annual_uptime * 8760
    annual_production_g     = capacity_kta * 1e9
    annual_production_gross = annual_production_g / batches_on_spec
    tank_working_vol        = tank_volume_L * WORKING_VOL_RATIO

    batch_cycle_time  = fermentation['ferm_time'] + turnaround_time
    batches_per_tank  = annual_uptime_hr / batch_cycle_time
    n_tanks = math.ceil(
        annual_production_gross
        / (fermentation['titer'] * tank_working_vol * batches_per_tank)
    )

    annual_ferm_vol   = n_tanks * batches_per_tank * tank_working_vol
    total_batches     = n_tanks * batches_per_tank
    annual_biomass_kg = annual_ferm_vol * fermentation['final_biomass'] / 1e6
    annual_O2_mmol    = fermentation['cumulative_O2'] * annual_ferm_vol

    return {
        'annual_uptime_hr':        annual_uptime_hr,
        'annual_production_g':     annual_production_g,
        'annual_production_gross': annual_production_gross,
        'tank_working_vol':        tank_working_vol,
        'batch_cycle_time':        batch_cycle_time,
        'batches_per_tank':        batches_per_tank,
        'n_tanks':                 n_tanks,
        'annual_ferm_vol':         annual_ferm_vol,
        'total_batches':           total_batches,
        'annual_biomass_kgCDW':    annual_biomass_kg,
        'annual_O2_mmol':          annual_O2_mmol,
    }


def calculate_raw_material_costs(logistics, fermentation, chemistry,
                                  DSP_yield, price_glucose_per_g,
                                  price_ammonia_per_g, media_cost_per_kgCDW,
                                  price_NaOH_per_kg, price_peracetic_per_L,
                                  price_mgso4_per_kg=0.0):
    """Calculate annual raw material costs ($/yr)."""
    annual_production_g  = logistics['annual_production_g']
    annual_ferm_vol      = logistics['annual_ferm_vol']
    annual_biomass_kgCDW = logistics['annual_biomass_kgCDW']
    total_batches        = logistics['total_batches']
    tank_working_vol     = logistics['tank_working_vol']

    glucose_per_g_product = 1.0 / (fermentation['overall_yield'] * DSP_yield)
    cost_glucose = annual_production_g * glucose_per_g_product * price_glucose_per_g

    d = chemistry['atoms']['N']
    if d > 0 and chemistry['yield_coeffs']['NH3'] is not None:
        NH3_for_product_g = annual_production_g / chemistry['yield_coeffs']['NH3']
    else:
        NH3_for_product_g = 0.0
    NH3_for_biomass_g = annual_biomass_kgCDW * 1000 * (MW_NH3 / MW_BIOMASS)
    cost_ammonia = (NH3_for_product_g + NH3_for_biomass_g) * price_ammonia_per_g

    if chemistry['yield_coeffs'].get('H2SO4') is not None:
        H2SO4_for_product_g  = annual_production_g / chemistry['yield_coeffs']['H2SO4']
        MgSO4_for_product_kg = H2SO4_for_product_g * (MW_MGSO4 / MW_H2SO4) / 1000.0
        cost_sulfate = MgSO4_for_product_kg * price_mgso4_per_kg
    else:
        cost_sulfate = 0.0

    cost_media = annual_biomass_kgCDW * media_cost_per_kgCDW

    NaOH_per_batch_kg     = tank_working_vol * 0.02
    peracetic_per_batch_L = tank_working_vol * (0.02 / 20.0)
    cost_CIP = total_batches * (NaOH_per_batch_kg * price_NaOH_per_kg
                                + peracetic_per_batch_L * price_peracetic_per_L)

    return {
        'glucose': cost_glucose,
        'ammonia': cost_ammonia,
        'sulfate': cost_sulfate,
        'media':   cost_media,
        'CIP':     cost_CIP,
    }


def calculate_utility_costs(logistics, fermentation, ferm_temp_C,
                             tank_volume_L, price_electricity,
                             price_natural_gas, CEPCI, cost_of_fuel):
    """Calculate annual utility costs ($/yr)."""
    annual_uptime_hr  = logistics['annual_uptime_hr']
    annual_ferm_vol   = logistics['annual_ferm_vol']
    annual_O2_mmol    = logistics['annual_O2_mmol']
    annual_biomass_kg = logistics['annual_biomass_kgCDW']

    annual_media_m3   = annual_ferm_vol / 1000
    demand_rate_m3s   = (annual_media_m3 / annual_uptime_hr) / 3600
    water_cost_per_m3 = ((0.0007 + 0.00003 * demand_rate_m3s**-0.6) * CEPCI
                         + 0.02 * cost_of_fuel)
    cost_water = annual_media_m3 * water_cost_per_m3

    annual_air_m3   = annual_O2_mmol / (O2_MOLES_PER_M3_AIR * 1000 * DO_UTILISATION)
    ave_airflow_m3s = annual_air_m3 / (annual_uptime_hr * 3600)

    tank_total_vol_m3      = tank_volume_L / 1000
    tank_radius            = (tank_total_vol_m3 / (math.pi * ASPECT_RATIO)) ** (1/3)
    tank_height            = ASPECT_RATIO * 2 * tank_radius
    max_ferm_pressure_psig = ((1 + tank_height * 1000 * 9.81 / 101325) - 1) * 14.696

    air_cost_per_m3 = ((0.00005 * ave_airflow_m3s**-0.3
                        * math.log10(max_ferm_pressure_psig) * CEPCI)
                       + 0.0009 * math.log10(max_ferm_pressure_psig) * cost_of_fuel)
    cost_compressed_air = air_cost_per_m3 * annual_air_m3

    annual_air_kg       = annual_air_m3 * 1.225
    cost_mass_transfer  = 1.8 * 0.233 * annual_air_kg * price_electricity

    cum_cooling_kJ = 0.460 * annual_O2_mmol
    if ferm_temp_C > 33:
        dT = (ferm_temp_C - 4.0) - 29.0
    else:
        dT = 30.0 - 29.0
    annual_cooling_m3    = cum_cooling_kJ / (4.184 * dT * 1000)
    ave_cooling_flow_m3s = annual_cooling_m3 / (annual_uptime_hr * 3600)
    cost_cooling = annual_cooling_m3 * ((0.0001 + 0.00003 / ave_cooling_flow_m3s)
                                        * CEPCI + 0.003 * cost_of_fuel)

    steri_energy_kJ     = 4.184 * annual_ferm_vol * 1.05 * 0.2 * (120 - 25)
    cost_sterilisation  = (steri_energy_kJ / 1.055056e6) * price_natural_gas

    heat_kill_vol_L     = annual_ferm_vol * fermentation['final_biomass'] * 1.6 / 1e6
    heat_kill_energy_kJ = 4.184 * heat_kill_vol_L * 1.1 * 0.2 * (60 - 25)
    cost_heat_kill      = (heat_kill_energy_kJ / 1.055056e6) * price_natural_gas

    centrifuge_flow_rate = 10000
    hours_centrifuge     = annual_ferm_vol / centrifuge_flow_rate
    n_centrifuges        = math.ceil(hours_centrifuge / annual_uptime_hr)
    actual_uptime_cent   = annual_ferm_vol / (n_centrifuges * centrifuge_flow_rate)
    power_per_cent_kW    = 0.3 * 10
    cost_centrifuge      = (power_per_cent_kW * n_centrifuges
                            * actual_uptime_cent * price_electricity)

    return {
        'water':              cost_water,
        'compressed_air':     cost_compressed_air,
        'mass_transfer':      cost_mass_transfer,
        'cooling_water':      cost_cooling,
        'sterilisation':      cost_sterilisation,
        'heat_kill':          cost_heat_kill,
        'centrifugation':     cost_centrifuge,
        'n_centrifuges':      n_centrifuges,
        'annual_air_m3':      annual_air_m3,
        'annual_cooling_m3':  annual_cooling_m3,
        'steri_energy_kJ':    steri_energy_kJ,
        'heat_kill_energy_kJ': heat_kill_energy_kJ,
        'max_ferm_pressure_psig': max_ferm_pressure_psig,
        'ave_airflow_m3s':    ave_airflow_m3s,
    }


def calculate_labour_cost(n_tanks):
    """Estimate annual labour cost ($/yr) based on fermenter count."""
    operators_per_tank  = 2.5
    loaded_cost_per_FTE = 60_000
    overhead_factor     = 1.5
    return n_tanks * operators_per_tank * loaded_cost_per_FTE * overhead_factor


def calculate_opex(logistics, fermentation, chemistry,
                   DSP_yield, DSP_OPEX_frac,
                   price_glucose_per_g, price_ammonia_per_g,
                   media_cost_per_kgCDW, price_NaOH_per_kg,
                   price_peracetic_per_L, price_mgso4_per_kg,
                   price_electricity,
                   price_natural_gas, CEPCI, cost_of_fuel,
                   ferm_temp_C, tank_volume_L,
                   other_fixed_costs=0.0):
    """
    Calculate total annual OPEX broken down by category ($/yr).

    Pass other_fixed_costs=0 on first call; update after CAPEX is known
    (other_fixed = 3.7% of TIC).
    """
    rm   = calculate_raw_material_costs(
               logistics, fermentation, chemistry, DSP_yield,
               price_glucose_per_g, price_ammonia_per_g,
               media_cost_per_kgCDW, price_NaOH_per_kg, price_peracetic_per_L,
               price_mgso4_per_kg)
    util = calculate_utility_costs(
               logistics, fermentation, ferm_temp_C, tank_volume_L,
               price_electricity, price_natural_gas, CEPCI, cost_of_fuel)
    labour = calculate_labour_cost(logistics['n_tanks'])

    ferm_opex = (rm['glucose'] + rm['ammonia'] + rm['sulfate'] + rm['media'] + rm['CIP']
                 + util['water'] + util['compressed_air'] + util['mass_transfer']
                 + util['cooling_water'] + util['sterilisation'] + util['heat_kill']
                 + util['centrifugation'] + labour + other_fixed_costs)

    total_opex = ferm_opex / (1 - DSP_OPEX_frac)
    DSP_opex   = DSP_OPEX_frac * total_opex

    annual_production_g = logistics['annual_production_g']
    capacity_kg         = annual_production_g / 1000

    return {
        'glucose':         rm['glucose'],
        'ammonia':         rm['ammonia'],
        'sulfate':         rm['sulfate'],
        'media':           rm['media'],
        'CIP':             rm['CIP'],
        'water':           util['water'],
        'compressed_air':  util['compressed_air'],
        'mass_transfer':   util['mass_transfer'],
        'cooling_water':   util['cooling_water'],
        'sterilisation':   util['sterilisation'],
        'heat_kill':       util['heat_kill'],
        'centrifugation':  util['centrifugation'],
        'labour':          labour,
        'other_fixed':     other_fixed_costs,
        'ferm_opex':       ferm_opex,
        'DSP_opex':        DSP_opex,
        'total_opex':      total_opex,
        'opex_per_kg':     total_opex / capacity_kg,
        'n_centrifuges':   util['n_centrifuges'],
        'annual_air_m3':   util['annual_air_m3'],
        'annual_cooling_m3': util['annual_cooling_m3'],
        'steri_energy_kJ':    util['steri_energy_kJ'],
        'heat_kill_energy_kJ': util['heat_kill_energy_kJ'],
        'max_ferm_pressure_psig': util['max_ferm_pressure_psig'],
        'ave_airflow_m3s': util['ave_airflow_m3s'],
    }


# ════════════════════════════════════════════════════════════════════════════════
# SECTION 4 — CAPEX
# ════════════════════════════════════════════════════════════════════════════════

# Equipment cost database — Table S4.1, Lynch 2021
# Format: (quoted_cost $, quoted_size, scaling_exp, inflation_factor, install_factor)
# quoted_size = None -> fixed cost per unit (no scaling)
EQUIP_DB = {
    'fermenter':      (176000,   None,   None, 1.13, 2.0),
    'agitator':       (36000,    36,     0.50, 1.00, 1.5),
    'main_pump':      (3900,     0.01,   0.80, 1.17, 2.3),
    'feed_tank':      (70000,    100,    0.70, 1.17, 2.6),
    'feed_pump':      (3900,     0.01,   0.80, 1.17, 2.3),
    'base_tank':      (98000,    100,    0.70, 1.13, 1.5),
    'base_pump':      (3900,     0.01,   0.80, 1.17, 2.3),
    'acid_tank':      (196000,   100,    0.70, 1.13, 2.0),
    'acid_pump':      (3900,     0.01,   0.80, 1.17, 2.3),
    'dry_chem':       (100000,   None,   None, 1.00, 2.0),
    'media_prep':     (91200,    100,    0.70, 1.17, 2.6),
    'media_pump':     (3900,     0.01,   0.80, 1.17, 2.3),
    'CIP_tank':       (98000,    10,     0.70, 1.13, 2.0),
    'CIP_pump':       (3900,     0.01,   0.80, 1.17, 2.3),
    'centrifuge':     (325000,   None,   None, 1.59, 1.8),
    'broth_tank':     (1317000,  1000,   0.70, 1.00, 1.8),
    'broth_pump':     (3900,     0.01,   0.80, 1.17, 2.3),
    'cooling_tower':  (1375000,  0.1,    0.60, 1.12, 1.5),
    'cooling_pump':   (283671,   0.1,    0.80, 1.12, 3.1),
    'boiler':         (100000,   1000,   0.60, 1.59, 2.0),
    'air_compressor': (1318600,  1318.6, 1.00, 1.03, 1.6),
    'air_receiver':   (17000,    10,     0.70, 1.12, 3.1),
    'air_dryer':      (15000,    0.1,    0.60, 1.17, 1.8),
    'water_tank':     (250000,   500,    0.60, 1.17, 1.7),
    'water_softener': (78000,    500,    0.60, 1.17, 1.8),
    'water_pump':     (15292,    0.01,   0.60, 1.17, 3.1),
    'ww_tank':        (1317000,  1000,   0.60, 1.00, 1.8),
    'ww_pump':        (3900,     0.01,   0.60, 1.17, 2.3),
    'potable_water':  (75000,    None,   None, 1.00, 1.7),
    'heat_exchanger': (15000,    10,     0.50, 1.59, 3.1),
}


def equipment_TIC(name, actual_size=None, n_units=1):
    """
    Total Installed Cost for a piece of equipment.
    Implements Equations S4.1–S4.2, Lynch 2021.
    """
    qc, qs, sf, inf_f, inst_f = EQUIP_DB[name]
    if qs is None or sf is None:
        purchase = qc
    else:
        purchase = qc * (actual_size / qs) ** sf
    return inf_f * purchase * inst_f * n_units


def size_equipment(logistics, fermentation, opex_results, tank_volume_L, ferm_temp_C):
    """
    Calculate sizes for all major equipment items.
    Returns dict of sizes and intermediates for CAPEX calculations.
    """
    n_tanks          = logistics['n_tanks']
    tank_working_vol = logistics['tank_working_vol']
    annual_uptime_hr = logistics['annual_uptime_hr']
    annual_ferm_vol  = logistics['annual_ferm_vol']

    tank_total_vol_m3 = tank_volume_L / 1000
    tank_radius = (tank_total_vol_m3 / (math.pi * ASPECT_RATIO)) ** (1/3)
    tank_height = ASPECT_RATIO * 2 * tank_radius

    max_ferm_P_Pa   = 101325 + tank_height * 1000 * 9.81
    max_ferm_P_psig = (max_ferm_P_Pa / 101325 - 1) * 14.696

    pump_flow_m3s = tank_working_vol / 1000 / 3600
    agitator_kW   = 2.0 * (tank_working_vol / 1000)

    feed_vol_m3  = tank_working_vol * 0.02 * 12 / 1000
    base_vol_m3  = tank_working_vol * 0.005 * 12 / 1000
    CIP_vol_m3   = tank_total_vol_m3 / 100
    broth_vol_m3 = n_tanks * tank_total_vol_m3 * 0.5

    ave_cool_m3s = opex_results['annual_cooling_m3'] / (annual_uptime_hr * 3600)

    steri_time_hr  = 0.05 * annual_uptime_hr
    steri_kW       = opex_results['steri_energy_kJ'] / (steri_time_hr * 3600)
    steam_kg_hr    = steri_kW / 1910 * 3600
    steam_lb_hr    = steam_kg_hr * 2.2046

    y       = 1.4
    T_in    = 298.15
    P_in    = 101325
    P_out   = max_ferm_P_Pa
    MW_air  = 28.96
    mol_flow = (opex_results['ave_airflow_m3s'] * 1.225 * 1000) / MW_air
    W_s_kW   = ((y/(y-1)) * 8.314 * T_in * ((P_out/P_in)**((y-1)/y) - 1)
                * mol_flow / 1000)
    comp_kW  = W_s_kW / 0.70

    denom       = max(50 - 1.5 * max_ferm_P_psig, 1.0)
    receiver_m3 = (30 * opex_results['ave_airflow_m3s'] * 14.7) / denom

    U_kW_m2K = 1.0
    def hx_area(energy_kJ, T_hot, time_hr):
        if energy_kJ <= 0 or T_hot <= 26:
            return 1.0
        Q_kW = energy_kJ / (time_hr * 3600)
        dT1 = 201 - T_hot
        dT2 = T_hot - 25
        lmtd = (dT1 - dT2) / math.log(dT1/dT2) if abs(dT1-dT2) > 0.01 else dT1
        return max(Q_kW / (U_kW_m2K * lmtd), 1.0)

    hx_steri_m2    = hx_area(opex_results['steri_energy_kJ'], 120, steri_time_hr)
    hx_heatkill_m2 = hx_area(opex_results['heat_kill_energy_kJ'], 60, steri_time_hr)

    water_m3s    = (annual_ferm_vol / 1000) / (annual_uptime_hr * 3600)
    water_4hr_m3 = water_m3s * 3600 * 4
    water_2hr_m3 = water_m3s * 3600 * 2
    ww_vol_m3    = max(fermentation['final_biomass'] * annual_ferm_vol * 1.6 / 1e9, 10)

    return {
        'n_tanks': n_tanks, 'n_centrifuges': opex_results['n_centrifuges'],
        'pump_flow_m3s': pump_flow_m3s, 'agitator_kW': agitator_kW,
        'feed_vol_m3': feed_vol_m3, 'base_vol_m3': base_vol_m3,
        'CIP_vol_m3': CIP_vol_m3, 'broth_vol_m3': broth_vol_m3,
        'ave_cool_m3s': ave_cool_m3s, 'steam_lb_hr': steam_lb_hr,
        'comp_kW': comp_kW, 'receiver_m3': receiver_m3,
        'hx_steri_m2': hx_steri_m2, 'hx_heatkill_m2': hx_heatkill_m2,
        'water_4hr_m3': water_4hr_m3, 'water_2hr_m3': water_2hr_m3,
        'water_m3s': water_m3s, 'ww_vol_m3': ww_vol_m3,
        'max_ferm_P_psig': max_ferm_P_psig,
        'ave_airflow_m3s': opex_results['ave_airflow_m3s'],
        'tank_total_vol_m3': tank_total_vol_m3,
        'steri_time_hr': steri_time_hr,
    }


def calculate_capex(sizing, DSP_CAPEX_frac):
    """
    Calculate total capital costs from equipment sizing.
    Implements Table 1 CAPEX rollup from Lynch 2021.
    """
    n  = sizing['n_tanks']
    nc = sizing['n_centrifuges']
    pf = sizing['pump_flow_m3s']

    a1_fermenters = equipment_TIC('fermenter', n_units=n)
    a1_agitators  = equipment_TIC('agitator',  sizing['agitator_kW'], n_units=n)
    a1_pumps      = equipment_TIC('main_pump', pf, n_units=n)
    a1_feed       = (equipment_TIC('feed_tank', sizing['feed_vol_m3'],  n_units=n)
                   + equipment_TIC('feed_pump', pf*0.02, n_units=n))
    a1_base       = (equipment_TIC('base_tank', sizing['base_vol_m3'],  n_units=n)
                   + equipment_TIC('base_pump', pf*0.005, n_units=n))
    a1_acid       = (equipment_TIC('acid_tank', sizing['base_vol_m3'],  n_units=n)
                   + equipment_TIC('acid_pump', pf*0.005, n_units=n))
    a1_media      = (equipment_TIC('dry_chem')
                   + equipment_TIC('media_prep', sizing['tank_total_vol_m3'])
                   + equipment_TIC('media_pump', pf))
    a1_CIP        = (equipment_TIC('CIP_tank', sizing['CIP_vol_m3'], n_units=3)
                   + equipment_TIC('CIP_pump', pf*0.01, n_units=3))

    area1_equip  = (a1_fermenters + a1_agitators + a1_pumps
                    + a1_feed + a1_base + a1_acid + a1_media + a1_CIP)
    area1_piping = 0.045 * area1_equip
    area1_total  = area1_equip + area1_piping

    area2_total = 0.27 * area1_total

    area3_equip  = (equipment_TIC('centrifuge', n_units=nc)
                  + equipment_TIC('broth_tank', sizing['broth_vol_m3'])
                  + equipment_TIC('broth_pump', pf * nc))
    area3_piping = 0.045 * area3_equip
    area3_total  = area3_equip + area3_piping

    a4_cooling = (equipment_TIC('cooling_tower', sizing['ave_cool_m3s'])
                + equipment_TIC('cooling_pump',  sizing['ave_cool_m3s']))
    a4_steam   =  equipment_TIC('boiler', sizing['steam_lb_hr'])
    a4_air     = (equipment_TIC('air_compressor', sizing['comp_kW'])
                + equipment_TIC('air_receiver',   max(abs(sizing['receiver_m3']), 1))
                + equipment_TIC('air_dryer',      sizing['ave_airflow_m3s']))
    a4_water   = (equipment_TIC('water_tank',     sizing['water_4hr_m3'])
                + equipment_TIC('water_softener', sizing['water_2hr_m3'])
                + equipment_TIC('water_pump',     sizing['water_m3s'])
                + equipment_TIC('potable_water'))
    a4_ww      = (equipment_TIC('ww_tank', sizing['ww_vol_m3'])
                + equipment_TIC('ww_pump', sizing['water_m3s'] * 0.1))
    a4_HX      = (equipment_TIC('heat_exchanger', sizing['hx_steri_m2'])
                + equipment_TIC('heat_exchanger', sizing['hx_heatkill_m2']))

    area4_equip   = a4_cooling + a4_steam + a4_air + a4_water + a4_ww + a4_HX
    area4_piping  = 0.045 * area4_equip
    TIC_elec_ctrl = 0.10 * (area1_equip + area3_equip + area4_equip)
    area4_total   = area4_equip + area4_piping + TIC_elec_ctrl

    TIC_upstream = area1_total + area2_total + area3_total + area4_total

    site_dev  = 0.09 * TIC_upstream
    warehouse = 0.04 * TIC_upstream
    admin     = 0.05 * TIC_upstream
    TDC       = TIC_upstream + site_dev + warehouse + admin

    indirect  = (0.10 + 0.10 + 0.20 + 0.10 + 0.10) * TDC
    FCI       = TDC + indirect
    WC        = 0.05 * FCI
    TCI       = FCI + WC

    TCI_total = TCI / (1 - DSP_CAPEX_frac)
    DSP_capex = DSP_CAPEX_frac * TCI_total

    return {
        'area1': area1_total, 'area2': area2_total,
        'area3': area3_total, 'area4': area4_total,
        'a1_fermenters': a1_fermenters, 'a1_agitators': a1_agitators,
        'a4_cooling': a4_cooling, 'a4_steam': a4_steam,
        'a4_air': a4_air, 'a4_water': a4_water,
        'TIC_elec_ctrl': TIC_elec_ctrl,
        'TIC_upstream': TIC_upstream, 'DSP_capex': DSP_capex,
        'TDC': TDC, 'indirect': indirect,
        'FCI': FCI, 'WC': WC,
        'TCI': TCI, 'TCI_total': TCI_total,
        'comp_kW': sizing['comp_kW'],
        'steam_lb_hr': sizing['steam_lb_hr'],
    }


# ════════════════════════════════════════════════════════════════════════════════
# SECTION 5 — FINANCIALS
# ════════════════════════════════════════════════════════════════════════════════

def build_loan_schedule(principal, annual_rate, term_yr, n_years):
    """
    Year-by-year loan amortisation schedule (equal annual payments).
    Returns (interest_payments, principal_payments), each a list of length n_years.
    """
    balance        = principal
    annual_payment = (principal * annual_rate
                      / (1 - (1 + annual_rate)**(-term_yr))
                      if annual_rate > 0 else principal / term_yr)
    interest_list  = []
    principal_list = []
    for yr in range(n_years):
        if balance > 0:
            interest = balance * annual_rate
            princ    = min(annual_payment - interest, balance)
            balance  = max(balance - princ, 0.0)
        else:
            interest = 0.0
            princ    = 0.0
        interest_list.append(interest)
        principal_list.append(princ)
    return interest_list, principal_list


def calculate_MSP(opex_total, capex, capacity_kg,
                  target_margin, tax_rate,
                  pct_debt, loan_interest, loan_term_yr,
                  construction_yr, depreciation_yr,
                  ongoing_capex_frac,
                  ramp_fractions=None):
    """
    Calculate the Minimum Selling Price ($/kg) at nameplate capacity.
    Evaluated at the first full-production year (after construction + ramp-up).
    """
    if ramp_fractions is None:
        ramp_fractions = RAMP_FRACTIONS

    TCI            = capex['TCI_total']
    debt_amount    = pct_debt * TCI
    annual_deprec  = TCI / depreciation_yr
    annual_ongoing = ongoing_capex_frac * TCI

    ramp_yr       = len(ramp_fractions)
    n_years_to_np = construction_yr + ramp_yr
    total_yrs     = n_years_to_np + 1

    interest_s, _ = build_loan_schedule(debt_amount, loan_interest,
                                         loan_term_yr, total_yrs)
    interest_at_nameplate = interest_s[n_years_to_np]

    numerator   = (opex_total + annual_deprec + annual_ongoing
                   + interest_at_nameplate) * (1 - tax_rate)
    denominator = (1 - tax_rate) - target_margin

    if denominator <= 0:
        return float('inf')

    MSP_revenue = numerator / denominator
    return MSP_revenue / capacity_kg


def calculate_DCF(opex_total, capex, capacity_kg, selling_price,
                   tax_rate, discount_rate, payback_period,
                   pct_debt, loan_interest, loan_term_yr,
                   construction_yr, ramp_fractions,
                   capex_yr1_frac, capex_yr2_frac,
                   ongoing_capex_frac, depreciation_yr):
    """
    Build a full discounted cash flow proforma.
    Returns dict with NPV, IRR, ROI, and year-by-year cash flow arrays.
    """
    TCI            = capex['TCI_total']
    debt_amount    = pct_debt * TCI
    equity_amount  = (1 - pct_debt) * TCI
    annual_deprec  = TCI / depreciation_yr
    annual_ongoing = ongoing_capex_frac * TCI
    ramp_yr        = len(ramp_fractions)
    total_yrs      = construction_yr + payback_period

    interest_s, principal_s = build_loan_schedule(
        debt_amount, loan_interest, loan_term_yr, total_yrs)

    cash_flows  = []
    revenues    = []
    net_incomes = []

    for yr in range(total_yrs):
        if yr < construction_yr:
            frac = capex_yr1_frac if yr == 0 else capex_yr2_frac
            ncf  = -(frac * equity_amount + interest_s[yr])
            cash_flows.append(ncf)
            revenues.append(0)
            net_incomes.append(0)
            continue

        prod_yr   = yr - construction_yr
        prod_frac = ramp_fractions[prod_yr] if prod_yr < ramp_yr else 1.0

        revenue = selling_price * capacity_kg * prod_frac
        opex_yr = opex_total * prod_frac
        deprec  = annual_deprec if prod_yr < depreciation_yr else 0.0

        EBIT  = revenue - opex_yr - deprec - annual_ongoing
        EBT   = EBIT - interest_s[yr]
        tax   = max(EBT * tax_rate, 0.0)
        NI    = EBT - tax

        ncf = NI + deprec - principal_s[yr]
        cash_flows.append(ncf)
        revenues.append(revenue)
        net_incomes.append(NI)

    NPV = sum(cf / (1 + discount_rate)**yr for yr, cf in enumerate(cash_flows))

    def npv_at(r):
        return sum(cf / (1+r)**yr for yr, cf in enumerate(cash_flows))

    try:
        lo, hi = -0.9999, 10.0
        for _ in range(400):
            mid = (lo + hi) / 2
            (lo := mid) if npv_at(mid) > 0 else (hi := mid)
        IRR = (lo + hi) / 2
    except Exception:
        IRR = float('nan')

    total_return = sum(cash_flows[construction_yr:])
    ROI = total_return / TCI * 100

    cum_flows = []
    running = 0
    for cf in cash_flows:
        running += cf
        cum_flows.append(running)

    return {
        'NPV':           NPV,
        'IRR':           IRR * 100,
        'ROI':           ROI,
        'cash_flows':    cash_flows,
        'cum_flows':     cum_flows,
        'revenues':      revenues,
        'net_incomes':   net_incomes,
        'total_yrs':     total_yrs,
        'construction_yr': construction_yr,
    }
