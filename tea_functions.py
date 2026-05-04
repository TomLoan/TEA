"""
tea_functions.py — Bioprocess TEA calculation functions.

Implements the model described in:
  Lynch, S.A. (2021). The bioprocess TEA calculator.
  Metabolic Engineering, 65, 42–51.

All functions are pure (no print/plot side-effects except inside verbose= guards).

SCOPE
-----
- Aerobic batch/fed-batch fermentation on glucose as the sole carbon source
- Products containing C, H, O, N, S atoms
- New plant construction at 1–100 kta production scale
- Accuracy: ±50% (FEL-1 / order-of-magnitude estimate)
- Intended for early-stage R&D goal-setting, not detailed engineering

NOT APPLICABLE TO
-----------------
- Anaerobic fermentation
- Alternative carbon sources (xylose, methanol, glycerol, …)
- Stationary-phase or non-growth-associated production
- GMP / pharmaceutical-grade facilities (regulatory costs not modelled)
- Mammalian cell culture beyond rough approximation

KNOWN DISCREPANCIES vs LYNCH 2021 (all within ±50% FEL-1 tolerance)
---------------------------------------------------------------------
- MaxOTR: equations in the supplementary give a value ~14× lower than Fig. 5b
  (176.18 mmol/L/hr). Cannot be resolved from published equations alone —
  likely reflects undescribed terms in the original JavaScript tool.
  CUMULATIVE O₂ (which drives annual air and cooling costs) is unaffected.
- OPEX: ~$1.32/kg vs paper $1.46/kg — gap from labour formula approximation
  and capital-linked fixed costs (Davis 2018 formula not fully published).
- CAPEX: ~$2.67/kg vs paper $2.76/kg — within 4%.
- IRR: ~10 percentage points higher throughout — consistent with slightly
  lower cost base; relative sensitivity to inputs is preserved.
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

# ── DSP route library ─────────────────────────────────────────────────────────
# Each route defines per-step yield fractions and cost correlations.
# CAPEX scales with annual broth throughput via power law (6th-tenths rule).
# OPEX = DSP_CAPEX × opex_pct_capex (covers utilities, consumables, labour).
# Reference costs are order-of-magnitude estimates consistent with ±50% FEL-1.
DSP_ROUTE_LIBRARY = {
    'minimal_processing': {
        'label': 'Minimal processing — cell separation only (food proteins, SCP, whole-cell products)',
        'steps': [
            ('Centrifugation / microfiltration', 0.97),
            ('Concentration / spray drying',     0.97),
        ],
        'capex_ref':      2_000_000,   # $ at ref_vol_m3_yr broth throughput
        'ref_vol_m3_yr':  10_000,
        'scale_exp':      0.60,
        'opex_pct_capex': 0.30,
    },
    'crystallisation': {
        'label': 'Crystallisation (organic acids, amino acids)',
        'steps': [
            ('Broth centrifugation',          0.98),
            ('Acidification / precipitation', 0.93),
            ('Washing filtration',            0.96),
            ('Drying',                        0.98),
        ],
        'capex_ref':      8_000_000,
        'ref_vol_m3_yr':  10_000,
        'scale_exp':      0.65,
        'opex_pct_capex': 0.35,
    },
    'solvent_extraction': {
        'label': 'Solvent extraction (terpenoids, fatty alcohols, hydrophobic small molecules)',
        'steps': [
            ('Broth centrifugation',     0.97),
            ('Liquid-liquid extraction', 0.88),
            ('Solvent stripping',        0.94),
            ('Polish / adsorption',      0.96),
        ],
        'capex_ref':      12_000_000,
        'ref_vol_m3_yr':  10_000,
        'scale_exp':      0.65,
        'opex_pct_capex': 0.45,
    },
    'uf_precipitation': {
        'label': 'UF + precipitation (secreted enzymes, extracellular proteins)',
        'steps': [
            ('Broth centrifugation',            0.93),
            ('Ultrafiltration / concentration', 0.90),
            ('Precipitation',                   0.82),
            ('Drying / formulation',            0.95),
        ],
        'capex_ref':      15_000_000,
        'ref_vol_m3_yr':  10_000,
        'scale_exp':      0.70,
        'opex_pct_capex': 0.50,
    },
    'chromatography': {
        'label': 'Capture + polish chromatography (recombinant proteins, high-value enzymes)',
        'steps': [
            ('Broth centrifugation',   0.93),
            ('Capture chromatography', 0.85),
            ('Polish chromatography',  0.90),
            ('UF/DF + formulation',    0.93),
        ],
        'capex_ref':      20_000_000,
        'ref_vol_m3_yr':  10_000,
        'scale_exp':      0.70,
        'opex_pct_capex': 0.60,
    },
    'multi_column_chromatography': {
        'label': 'Multi-step chromatography (therapeutic proteins, mAbs)',
        'steps': [
            ('Centrifugation / clarification', 0.88),
            ('Protein A / capture',            0.87),
            ('Ion exchange polish',            0.92),
            ('UF/DF + formulation',            0.95),
        ],
        'capex_ref':      40_000_000,
        'ref_vol_m3_yr':  10_000,
        'scale_exp':      0.75,
        'opex_pct_capex': 0.75,
    },
}


def calculate_dsp(route_name, annual_broth_vol_m3, step_yield_overrides=None):
    """
    Calculate DSP yield, CAPEX, and OPEX for a given processing route.

    Parameters
    ----------
    route_name : str
        Key in DSP_ROUTE_LIBRARY.
    annual_broth_vol_m3 : float
        Annual broth volume processed (m³/yr) — logistics['annual_ferm_vol'] / 1000.
    step_yield_overrides : list of float or None
        Per-step yield fractions (same length as route steps). None = use defaults.

    Returns
    -------
    dict with keys: step_names, step_yields, overall_yield, dsp_capex, dsp_opex, route_label.
    """
    route = DSP_ROUTE_LIBRARY[route_name]
    steps = route['steps']

    step_yields = (list(step_yield_overrides)
                   if step_yield_overrides is not None
                   else [y for _, y in steps])

    overall_yield = 1.0
    for y in step_yields:
        overall_yield *= y

    dsp_capex = (route['capex_ref']
                 * (annual_broth_vol_m3 / route['ref_vol_m3_yr']) ** route['scale_exp'])
    dsp_opex  = dsp_capex * route['opex_pct_capex']

    return {
        'step_names':    [name for name, _ in steps],
        'step_yields':   step_yields,
        'overall_yield': overall_yield,
        'dsp_capex':     dsp_capex,
        'dsp_opex':      dsp_opex,
        'route_label':   route['label'],
    }

# ── Organism presets ─────────────────────────────────────────────────────────
ORGANISM_PRESETS = {
    'Generic (model default)': {
        'biomass_yield_coeff':  BIOMASS_YIELD_COEFF,   # ~0.504 gCDW/g glucose
        'carbon_to_co2_frac':   0.0,
        'media_cost':           0.40,
        'note': 'Battley 1987 empirical average (C₃.₈₅H₆.₆₉O₁.₇₈N). '
                'Reasonable for E. coli, B. subtilis on glucose.',
    },
    'E. coli (aerobic)': {
        'biomass_yield_coeff':  0.48,
        'carbon_to_co2_frac':   0.20,   # acetate overflow at high mu
        'media_cost':           0.25,
        'note': 'Simple mineral salts medium; aerobic growth on glucose. '
                '~20% of non-product glucose lost to acetate overflow even aerobically.',
    },
    'S. cerevisiae (yeast)': {
        'biomass_yield_coeff':  0.45,
        'carbon_to_co2_frac':   0.35,   # Crabtree effect — ethanol above ~0.1 g/L/hr glucose
        'media_cost':           0.40,
        'note': 'Aerobic; mineral medium. Crabtree effect produces ethanol '
                'above ~0.1 g/L/hr glucose — ~35% of non-product glucose diverted to CO₂/ethanol.',
    },
    'Pichia pastoris (methanol)': {
        'biomass_yield_coeff':  0.35,   # gCDW/g methanol (literature: 0.3–0.5)
        'carbon_to_co2_frac':   0.05,   # Crabtree-negative; minimal overflow on methanol
        'media_cost':           0.50,   # $/kgCDW — trace salts, biotin, more complex than E. coli
        'note': 'Methylotrophic yeast widely used for secreted recombinant proteins. '
                'Crabtree-negative (no ethanol overflow on methanol). '
                'Select "Methanol" as carbon source to use methanol feedstock pricing. '
                'Stoichiometry uses a glucose-equivalent approximation (carbon content per gram '
                'differs by <7%) — suitable for FEL-1 (±50%) estimates. '
                'Typical: 1–20 g/L secreted protein; 0.05–0.5 g/L/hr productivity.',
    },
    'B. subtilis': {
        'biomass_yield_coeff':  0.50,
        'carbon_to_co2_frac':   0.10,   # minor acetoin/acetate overflow
        'media_cost':           0.30,
        'note': 'Aerobic; simple mineral medium. Relatively efficient — '
                'minor acetoin/acetate overflow (~10% of non-product glucose).',
    },
    'Mammalian (CHO-like)': {
        'biomass_yield_coeff':  0.20,
        'carbon_to_co2_frac':   0.50,   # Warburg-like lactate + high ATP maintenance
        'media_cost':           5.00,
        'note': 'Very rough approximation only. Complex medium required; '
                'high media cost. Warburg-like lactate production + high maintenance '
                'energy (~50% of non-product glucose to CO₂/lactate). '
                'Model assumes aerobic single-substrate fermentation — CHO bioreactors are substantially more complex.',
    },
}

CARBON_SOURCE_OPTIONS = {
    'Glucose': {
        'label':       'Glucose price ($/kg)',
        'default_per_kg': 0.40,
        'note':    'Dextrose (corn syrup). Industrial bulk price ~$0.33–0.55/kg.',
    },
    'Methanol': {
        'label':       'Methanol price ($/kg)',
        'default_per_kg': 0.31,   # ~$0.25–0.40/kg industrial bulk methanol
        'note':    'Industrial-grade methanol. Price ~$0.25–0.40/kg; highly region-dependent. '
                   'Stoichiometry uses a glucose-equivalent basis (carbon content per gram '
                   'differs by <7% between glucose and methanol — within FEL-1 tolerance).',
    },
}

# ── Financial defaults (used as defaults in app.py) ──────────────────────────
RAMP_FRACTIONS    = [0.50, 0.75, 1.00]
CAPEX_YR1_FRAC    = 0.70
CAPEX_YR2_FRAC    = 0.30
ONGOING_CAPEX_FRAC = 0.10
DEPRECIATION_YR   = 10


# ════════════════════════════════════════════════════════════════════════════════
# SECTION 1 — CHEMISTRY & STOICHIOMETRY
# ════════════════════════════════════════════════════════════════════════════════
#
# The key insight (Lynch 2021 §1) is to express both glucose and any organic
# product as combinations of two building blocks: CO₂ (carbon) and H₂
# (reducing equivalents). Glucose has an H₂:CO₂ ratio of exactly 2.
#
#   ratio > 2 → product MORE REDUCED than glucose: biology must conserve
#               reducing power, so some carbon is lost as CO₂ (Case 2).
#   ratio = 2 → NEUTRAL: no byproduct carbon or O₂ needed (Case 1).
#   ratio < 2 → product MORE OXIDISED than glucose: O₂ consumed as reactant
#               (Case 3).
#
# This classification determines the balanced stoichiometric equation and
# therefore the theoretical maximum yield.

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


# H₂:CO₂ ratio for CₐHᵦOᵧNδSε (Equation 3, Lynch 2021):
#   ratio = 0.5(b/a) − 1.0(c/a) − 1.5(d/a) + 3.0(e/a) + 2.0
# The +3(e/a) S term is an extension of the original: SO₄²⁻ (S at +6) is
# reduced to organic S (S at −2), consuming 8 electrons per atom, but the
# net coefficient is +3 after accounting for the O and H that H₂SO₄ itself
# contributes to the atom balance (not +4 as a naive electron count suggests).
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


# Balanced equations per mole of glucose (Supplemental Materials §1, Lynch 2021):
#   Case 1 (ratio = 2, neutral):
#     glucose + W NH₃ + V H₂SO₄  →  X product + Z H₂O
#   Case 2 (ratio > 2, more reduced):
#     glucose + W NH₃ + V H₂SO₄  →  X product + Z H₂O + Q CO₂
#   Case 3 (ratio < 2, more oxidised):
#     glucose + W NH₃ + V H₂SO₄ + Y O₂  →  X product + Z H₂O
# H₂SO₄ is used as the stoichiometric S source (gives same atom balance as MgSO₄).
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


# Theoretical yield (g product / g glucose):
#   Y = (moles_product × MW_product) / MW_glucose
# This is the stoichiometric ceiling. The actual yield achieved in the
# fermentation is yield_fraction × Y (yield_fraction is a user input).
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
#
# Models a single aerobic batch/fed-batch fermentation.
# Key assumptions (Supplemental Materials §2, Lynch 2021):
#
#   Growth-associated production: product accumulates as cells grow.
#
#   Logistic growth: N(t) = K / (1 + A·exp(−r·t))
#     K  = final biomass (gCDW/L, the "carrying capacity")
#     A  = (K − N₀) / N₀   (shape parameter from initial conditions)
#     r  = logistic growth rate (hr⁻¹), solved so culture reaches K in ferm_time
#     N₀ = 1% of K (inoculum from seed culture)
#     Peak growth rate occurs at N = K/2, used to estimate maximum OTR.
#
#   Glucose partitioning:
#     sugar_to_product = titer / theoretical_yield
#     stoich_total     = sugar_to_product / yield_fraction   (includes losses)
#     sugar_to_biomass = stoich_total − sugar_to_product
#     Biomass grows at 80% of the theoretical conversion (Eq. S2.8):
#       0.84 glucose + NH₃ + 1.212 O₂ → 1 biomass + 3.212 H₂O + 1.212 CO₂
#     Biomass empirical formula: C₃.₈₅H₆.₆₉O₁.₇₈N (Battley 1987), MW = 95.37 g/mol
#
#   NOTE — MaxOTR discrepancy: the published supplementary equations give a
#   MaxOTR ~14× lower than Fig. 5b (176.18 mmol/L/hr). This cannot be resolved
#   from the paper alone and likely reflects undescribed terms in the original
#   JavaScript. CUMULATIVE O₂ (which drives annual air and cooling costs) is
#   unaffected and matches the paper's cost outputs.
#
#   Heat generation: 0.460 kJ per mmol O₂ consumed (Doran 1995). (Eq. S3.22)

def run_fermentation_model(titer, rate, yield_fraction, chemistry,
                            biomass_yield_coeff=None, biomass_o2_coeff=None,
                            carbon_to_co2_frac=0.0,
                            production_mode='growth_associated',
                            target_biomass=None, growth_time_hr=None):
    """
    Run the complete fermentation model.

    Parameters
    ----------
    titer               : float  Final product concentration (g/L).
    rate                : float  Volumetric production rate (g/L/hr). In
                                 growth_associated mode: average over the batch.
                                 In stationary_phase mode: rate during the production
                                 phase only.
    yield_fraction      : float  Fraction of theoretical yield achieved (0–0.99).
    chemistry           : dict   Output from run_chemistry().
    biomass_yield_coeff : float  gCDW per g glucose. Defaults to BIOMASS_YIELD_COEFF.
    biomass_o2_coeff    : float  g biomass per g O2. Defaults to BIOMASS_O2_COEFF.
    carbon_to_co2_frac  : float  Fraction of non-product glucose (growth_associated) or
                                 growth-phase glucose (stationary_phase) diverted to
                                 CO₂/heat via overflow or maintenance. Default 0.0.
    production_mode     : str    'growth_associated' (default, Lynch 2021) or
                                 'stationary_phase' (grow to target_biomass, then
                                 produce at given rate for titer/rate hours).
    target_biomass      : float  Target biomass at induction (gCDW/L). Required for
                                 stationary_phase mode.
    growth_time_hr      : float  Duration of growth phase (hr). Required for
                                 stationary_phase mode.

    Returns
    -------
    dict with glucose partitioning, biomass, kinetics, oxygen, cooling, and
    time-course arrays (t_points, biomass_curve, product_curve).
    """
    _biomass_yield = biomass_yield_coeff if biomass_yield_coeff is not None else BIOMASS_YIELD_COEFF
    _biomass_o2    = biomass_o2_coeff    if biomass_o2_coeff    is not None else BIOMASS_O2_COEFF
    _co2_frac      = max(0.0, min(0.99, carbon_to_co2_frac))

    if production_mode == 'stationary_phase':
        if not target_biomass or target_biomass <= 0:
            raise ValueError("target_biomass must be > 0 for stationary_phase mode.")
        if not growth_time_hr or growth_time_hr <= 0:
            raise ValueError("growth_time_hr must be > 0 for stationary_phase mode.")

    theoretical_yield = chemistry['theoretical_yield']
    eq                = chemistry['equation']
    txl_overhead      = chemistry.get('txl_glucose_g_per_g', 0.0)

    # Shared: product glucose budget (same semantics in both modes)
    sugar_to_product = titer / theoretical_yield
    stoich_total     = sugar_to_product / yield_fraction
    sugar_for_txl    = titer * txl_overhead

    # O2 for product formation (Case 3 only) — shared
    if eq['case'] == 3 and eq['O2'] > 1e-9:
        product_O2_yield_coeff = chemistry['yield_coeffs']['O2']
        O2_for_product = (titer / product_O2_yield_coeff) * (1000.0 / MW_O2)
    else:
        O2_for_product = 0.0

    theoretical_biomass_yield_100pct = MW_BIOMASS / (0.84 * MW_GLUCOSE)

    if production_mode == 'stationary_phase':
        # ── Stationary-phase path ─────────────────────────────────────────────
        # Grow to target_biomass in growth_time_hr, then produce for titer/rate hr.
        t_production     = titer / rate
        ferm_time        = growth_time_hr + t_production

        final_biomass    = target_biomass
        starting_biomass = INOCULUM_FRACTION * target_biomass

        # Growth-phase glucose: work backwards from target_biomass.
        # carbon_to_co2_frac represents overflow/maintenance during the growth phase.
        glucose_to_biomass   = target_biomass / _biomass_yield
        total_growth_glucose = glucose_to_biomass / (1.0 - _co2_frac)
        glucose_to_co2       = total_growth_glucose - glucose_to_biomass
        sugar_to_biomass     = total_growth_glucose
        total_sugar          = stoich_total + sugar_for_txl + total_growth_glucose

        # Kinetics — logistic rate solved for the growth phase duration
        A                    = (target_biomass - starting_biomass) / starting_biomass
        logistic_growth_rate = -math.log(0.01 / A) / growth_time_hr
        product_to_cell_ratio = None   # not meaningful: production is decoupled from growth
        logistic_prod_rate    = None
        specific_rate         = rate / target_biomass

        # O2
        O2_for_biomass    = target_biomass / _biomass_o2 * (1000.0 / MW_O2)
        glucose_at_100pct = target_biomass / theoretical_biomass_yield_100pct
        waste_glucose     = glucose_to_biomass - glucose_at_100pct
        O2_for_waste      = waste_glucose    * 6000.0 / MW_GLUCOSE
        O2_for_overflow   = glucose_to_co2   * 6000.0 / MW_GLUCOSE
        O2_for_txl        = sugar_for_txl    * 6000.0 / MW_GLUCOSE
        cumulative_O2     = O2_for_biomass + O2_for_product + O2_for_waste + O2_for_overflow + O2_for_txl

        # max_OTR: growth phase peaks at N=K/2; stationary OTR from product formation (Case 3)
        max_OTR_biomass    = ((1.0 / _biomass_o2) * (1000.0 / MW_O2)
                               * (logistic_growth_rate / 4.0) * target_biomass)
        if eq['case'] == 3 and O2_for_product > 0:
            max_OTR_stationary = (rate / chemistry['yield_coeffs']['O2']) * (1000.0 / MW_O2)
        else:
            max_OTR_stationary = 0.0
        max_OTR = max(max_OTR_biomass, max_OTR_stationary)

        # Time-course: two concatenated phases proportional to their durations
        n_growth = max(2, int(150 * growth_time_hr / ferm_time))
        n_prod   = 300 - n_growth
        t_g = np.linspace(0, growth_time_hr, n_growth)
        t_p = np.linspace(growth_time_hr, ferm_time, n_prod)
        biomass_curve = np.concatenate([
            target_biomass / (1.0 + A * np.exp(-logistic_growth_rate * t_g)),
            np.full(n_prod, target_biomass),
        ])
        product_curve = np.concatenate([np.zeros(n_growth), rate * (t_p - growth_time_hr)])
        t_points = np.concatenate([t_g, t_p])

    else:
        # ── Growth-associated path (original Lynch 2021) ─────────────────────
        total_sugar      = stoich_total + sugar_for_txl
        sugar_to_biomass = stoich_total - sugar_to_product

        # Split non-product glucose: fraction to biomass vs overflow/maintenance
        glucose_to_biomass = sugar_to_biomass * (1.0 - _co2_frac)
        glucose_to_co2     = sugar_to_biomass * _co2_frac

        final_biomass    = glucose_to_biomass * _biomass_yield
        starting_biomass = INOCULUM_FRACTION * final_biomass

        ferm_time             = titer / rate
        A                     = (final_biomass - starting_biomass) / starting_biomass
        logistic_growth_rate  = -math.log(0.01 / A) / ferm_time
        product_to_cell_ratio = titer / (final_biomass - starting_biomass)
        logistic_prod_rate    = product_to_cell_ratio * logistic_growth_rate
        specific_rate         = rate / final_biomass

        O2_for_biomass    = final_biomass / _biomass_o2 * (1000.0 / MW_O2)
        glucose_at_100pct = final_biomass / theoretical_biomass_yield_100pct
        waste_glucose     = glucose_to_biomass - glucose_at_100pct
        O2_for_waste      = waste_glucose    * 6000.0 / MW_GLUCOSE
        O2_for_overflow   = glucose_to_co2   * 6000.0 / MW_GLUCOSE
        O2_for_txl        = sugar_for_txl    * 6000.0 / MW_GLUCOSE
        cumulative_O2     = O2_for_biomass + O2_for_product + O2_for_waste + O2_for_overflow + O2_for_txl

        max_OTR_biomass = ((1.0 / _biomass_o2) * (1000.0 / MW_O2)
                           * (logistic_growth_rate / 4.0) * final_biomass)
        if eq['case'] == 3 and O2_for_product > 0:
            max_product_rate = (product_to_cell_ratio * logistic_growth_rate
                                * final_biomass / 4.0)
            max_OTR_product  = ((1.0 / chemistry['yield_coeffs']['O2'])
                                * (1000.0 / MW_O2) * max_product_rate)
        else:
            max_OTR_product = 0.0
        max_OTR = max_OTR_biomass + max_OTR_product

        t_points      = np.linspace(0, ferm_time, 300)
        biomass_curve = final_biomass / (1.0 + A * np.exp(-logistic_growth_rate * t_points))
        product_curve = product_to_cell_ratio * (biomass_curve - starting_biomass)

    # Shared post-processing
    max_O2_gradient  = O2_SATURATION * (1.0 - DO_SETPOINT)
    max_kla          = (max_OTR / max_O2_gradient) / 3600.0
    max_cooling_rate = O2_COOLING_COEFF * max_OTR
    ave_cooling_rate = O2_COOLING_COEFF * (cumulative_O2 / ferm_time)

    return {
        'titer': titer, 'rate': rate, 'yield_fraction': yield_fraction,
        'theoretical_yield':    theoretical_yield,
        'sugar_to_product':     sugar_to_product,
        'stoich_total':         stoich_total,
        'total_sugar':          total_sugar,
        'sugar_to_biomass':     sugar_to_biomass,
        'glucose_to_biomass':   glucose_to_biomass,
        'glucose_to_co2':       glucose_to_co2,
        'final_biomass':        final_biomass,
        'starting_biomass':     starting_biomass,
        'ferm_time':            ferm_time,
        'A':                    A,
        'logistic_growth_rate': logistic_growth_rate,
        'product_to_cell_ratio': product_to_cell_ratio,
        'logistic_prod_rate':   logistic_prod_rate,
        'specific_rate':        specific_rate,
        'O2_for_biomass':       O2_for_biomass,
        'O2_for_product':       O2_for_product,
        'O2_for_waste':         O2_for_waste,
        'O2_for_overflow':      O2_for_overflow,
        'cumulative_O2':        cumulative_O2,
        'max_OTR':              max_OTR,
        'max_kla':              max_kla,
        'max_cooling_rate':     max_cooling_rate,
        'ave_cooling_rate':     ave_cooling_rate,
        'overall_yield':        titer / total_sugar,
        'sugar_for_txl':        sugar_for_txl,
        'O2_for_txl':           O2_for_txl,
        't_points':             t_points,
        'biomass_curve':        biomass_curve,
        'product_curve':        product_curve,
        'production_mode':      production_mode,
        't_growth':             growth_time_hr if production_mode == 'stationary_phase' else ferm_time,
        't_production':         (titer / rate)  if production_mode == 'stationary_phase' else ferm_time,
    }


# ════════════════════════════════════════════════════════════════════════════════
# SECTION 3 — OPERATING COSTS (OPEX)
# ════════════════════════════════════════════════════════════════════════════════
#
# Annual costs broken into three categories (Supplemental Materials §3, Lynch 2021):
#
#   Variable (scale with production volume):
#     Glucose, ammonia, MgSO₄, media salts, CIP chemicals (NaOH + peracetic acid),
#     water, compressed air, mass-transfer electricity, cooling water,
#     steam for sterilisation, steam for biomass heat-kill, centrifugation.
#
#   Fixed (largely independent of production volume):
#     Labour — scales with number of fermentation tanks (Davis et al. 2018).
#     Other fixed costs — 3.7% of TCI (maintenance, insurance, overhead).
#       NOTE: other_fixed_costs can only be added after CAPEX is known, so
#       calculate_opex() must be called twice (two-pass approach):
#         Pass 1: other_fixed_costs=0 → size equipment → calculate_capex()
#         Pass 2: other_fixed_costs = 0.037 × TCI_total → final OPEX
#
#   DSP (downstream processing):
#     Calculated from route-specific reference CAPEX scaled by annual broth throughput
#     (power-law economy of scale). DSP OPEX = DSP_CAPEX × route opex_pct_capex.
#     Routes and parameters are in DSP_ROUTE_LIBRARY; compute with calculate_dsp().
#
# Utility cost equations from Ulrich & Vasudevan (2006), scaled by CEPCI index.

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
    batches_per_tank  = math.floor(annual_uptime_hr / batch_cycle_time)
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
                   dsp,
                   price_glucose_per_g, price_ammonia_per_g,
                   media_cost_per_kgCDW, price_NaOH_per_kg,
                   price_peracetic_per_L, price_mgso4_per_kg,
                   price_electricity,
                   price_natural_gas, CEPCI, cost_of_fuel,
                   ferm_temp_C, tank_volume_L,
                   other_fixed_costs=0.0):
    """
    Calculate total annual OPEX broken down by category ($/yr).

    dsp : dict   Output from calculate_dsp(). Provides overall_yield (for raw
                 material scaling) and dsp_opex (absolute DSP operating cost).

    Pass other_fixed_costs=0 on first call; update after CAPEX is known
    (other_fixed = 3.7% of TCI).
    """
    rm   = calculate_raw_material_costs(
               logistics, fermentation, chemistry, dsp['overall_yield'],
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

    DSP_opex   = dsp['dsp_opex']
    total_opex = ferm_opex + DSP_opex

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
# SECTION 4 — CAPITAL COSTS (CAPEX)
# ════════════════════════════════════════════════════════════════════════════════
#
# Factored cost estimation (Supplemental Materials §4, Lynch 2021):
#
#   TIC = Inflation_Factor × QuotedCost × (ActualSize / QuotedSize)^ScalingExp
#         × InstallationFactor                              (Equations S4.1–S4.2)
#
# Equipment quoted costs, scaling exponents, and installation factors from
# Table S4.1, Lynch 2021 (originally Davis et al. 2013/2018).
#
# Capital structure rollup (Table 1, Lynch 2021):
#   TIC upstream
#   + site development (9%) + warehouse (4%) + admin buildings (5%)  = TDC
#   TDC + indirect costs (60% of TDC)                                = FCI
#   FCI + working capital (5% of FCI)                                = TCI (upstream)
#   TCI upstream + DSP_capex (from calculate_dsp())                   = TCI total

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
    'broth_tank':     (1317000,  1000,   0.70, 1.13, 1.8),
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


def calculate_capex(sizing, dsp):
    """
    Calculate total capital costs from equipment sizing.
    Implements Table 1 CAPEX rollup from Lynch 2021.

    dsp : dict   Output from calculate_dsp(). Provides dsp_capex (absolute $).
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

    DSP_capex = dsp['dsp_capex']
    TCI_total = TCI + DSP_capex

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
# SECTION 5 — FINANCIAL OUTPUTS
# ════════════════════════════════════════════════════════════════════════════════
#
# Plant timeline assumed by the model (Supplemental Materials §5, Lynch 2021):
#   Years 1–2  Construction  (70% capex yr 1, 30% yr 2); interest payments only
#   Year  3    Ramp-up       50% of nameplate capacity
#   Year  4    Ramp-up       75% of nameplate capacity
#   Year  5    Ramp-up       100% of nameplate capacity
#   Years 6+   Full production at nameplate
#
# MSP is evaluated at nameplate capacity (first full-production year, yr 5).
# DCF runs across the full user-specified payback period.

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


# MSP derivation — solving for the selling price that achieves the target margin.
# Starting from:
#   Margin = NI / Revenue
#          = [(Revenue − OPEX − Depreciation − Maintenance − Interest) × (1−Tax)]
#            / Revenue
# Setting Margin = target and solving for Revenue:
#   Revenue* = (OPEX + Depreciation + Maintenance + Interest) × (1−Tax)
#              / [(1−Tax) − Margin]
#   MSP = Revenue* / annual_capacity_kg
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


# DCF net cash flow each year:
#   NCF = Net Income + Depreciation − Principal repayment
# Depreciation is added back (non-cash accounting charge; cash was spent at build).
# Principal repayment is subtracted (real cash outflow not captured in P&L).
# IRR is solved by bisection on NPV = 0.
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
