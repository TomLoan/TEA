"""Headless smoke test: run the full TEA pipeline with app.py's default inputs.

Exercises every calculation in app.py's reactive block plus the dict keys the
tabs read, without starting a Streamlit server. Run with `python dryrun.py`.

The DEFAULTS block below must mirror the sidebar widget defaults in app.py.
Preset-derived values are read from tea_functions so they cannot drift; only
the literal widget values are repeated here.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from tea_functions import (
    run_chemistry, run_fermentation_model, calculate_plant_logistics,
    calculate_dsp, calculate_opex, size_equipment, calculate_capex,
    calculate_MSP, calculate_DCF,
    DSP_ROUTE_LIBRARY, ORGANISM_PRESETS, CARBON_SOURCE_OPTIONS,
    RAMP_FRACTIONS, CAPEX_YR1_FRAC, CAPEX_YR2_FRAC, ONGOING_CAPEX_FRAC,
    DEPRECIATION_YR,
)

# ── Defaults, mirroring app.py's sidebar ─────────────────────────────────────
# Product
formula          = "C7H12O4"
is_protein       = False
avg_residue_mw   = 110.0
atp_per_residue  = 5.0

# DSP - selectbox with no index, so the first route, at its default step yields
dsp_route_key    = list(DSP_ROUTE_LIBRARY)[0]
step_overrides   = [y for _name, y in DSP_ROUTE_LIBRARY[dsp_route_key]["steps"]]

# Organism - selectbox with no index, so the first preset
organism         = list(ORGANISM_PRESETS)[0]
_org             = ORGANISM_PRESETS[organism]
biomass_yield_override = float(_org["biomass_yield_coeff"])
mu_max_override        = float(_org["mu_max"])
# app.py renders this as an integer percent slider, so it round-trips
carbon_to_co2_override = round(_org["carbon_to_co2_frac"] * 100) / 100.0
media_cost_per_kgCDW   = _org["media_cost"]

is_stationary          = False        # "Growth-associated"
target_biomass_input   = 30.0
growth_time_input      = 24.0

# Fermentation performance
titer            = 150.0
rate             = 5.0
yield_fraction   = 0.90

# Plant configuration
capacity_kta     = 15.0
tank_volume_L    = 500_000
annual_uptime    = 0.90               # fraction of the year, not hours
ferm_temp_C      = 37.0
turnaround_time  = 16.0
batches_on_spec  = 0.95

# Raw material prices
carbon_source        = list(CARBON_SOURCE_OPTIONS)[0]
price_carbon_per_kg  = CARBON_SOURCE_OPTIONS[carbon_source]["default_per_kg"]
price_ammonia_per_kg = 0.26
price_mgso4_per_kg   = 0.30
price_NaOH_per_kg    = 0.15
price_peracetic_per_L= 5.00

# Utility prices
price_electricity    = 0.11
price_natural_gas    = 3.11
CEPCI                = 603

# Financial parameters
selling_price    = 2.50
target_margin    = 0.30
discount_rate    = 0.20
tax_rate         = 0.21
payback_period   = 20
pct_debt         = 0.60
loan_interest    = 0.08
loan_term_yr     = 10

# ── Derived inputs (mirrors app.py) ──────────────────────────────────────────
price_carbon_per_g  = price_carbon_per_kg / 1000
price_ammonia_per_g = price_ammonia_per_kg / 1000
cost_of_fuel        = price_natural_gas / 1.05505
capacity_kg         = capacity_kta * 1e6

OPEX_KWARGS = dict(
    price_feedstock_per_g=price_carbon_per_g,
    price_ammonia_per_g=price_ammonia_per_g,
    media_cost_per_kgCDW=media_cost_per_kgCDW,
    price_NaOH_per_kg=price_NaOH_per_kg,
    price_peracetic_per_L=price_peracetic_per_L,
    price_mgso4_per_kg=price_mgso4_per_kg,
    price_electricity=price_electricity,
    price_natural_gas=price_natural_gas,
    CEPCI=CEPCI,
    cost_of_fuel=cost_of_fuel,
    ferm_temp_C=ferm_temp_C,
    tank_volume_L=tank_volume_L,
)

# ── Pipeline ─────────────────────────────────────────────────────────────────
print("1. run_chemistry...")
chem = run_chemistry(formula, is_protein=is_protein,
                     avg_residue_mw=avg_residue_mw, atp_per_residue=atp_per_residue)

print("2. run_fermentation_model...")
ferm = run_fermentation_model(
    titer, rate, yield_fraction, chem,
    biomass_yield_coeff=biomass_yield_override,
    carbon_to_co2_frac=carbon_to_co2_override,
    mu_max=mu_max_override,
    production_mode="stationary_phase" if is_stationary else "growth_associated",
    target_biomass=target_biomass_input if is_stationary else None,
    growth_time_hr=growth_time_input    if is_stationary else None,
)

print("3. calculate_plant_logistics...")
logistics = calculate_plant_logistics(capacity_kta, annual_uptime, batches_on_spec,
                                      tank_volume_L, turnaround_time, ferm)

print("4. calculate_dsp...")
dsp = calculate_dsp(dsp_route_key, logistics["annual_ferm_vol"] / 1000,
                    step_yield_overrides=step_overrides)

print("5. calculate_opex pass 1...")
opex1 = calculate_opex(logistics, ferm, chem, dsp=dsp, **OPEX_KWARGS)

print("6. size_equipment...")
sizing = size_equipment(logistics, ferm, opex1, tank_volume_L, ferm_temp_C)

print("7. calculate_capex...")
capex = calculate_capex(sizing, dsp=dsp)

print("8. calculate_opex pass 2...")
opex = calculate_opex(logistics, ferm, chem, dsp=dsp,
                      other_fixed_costs=0.037 * capex["TCI_total"], **OPEX_KWARGS)

print("9. calculate_MSP...")
MSP = calculate_MSP(opex["total_opex"], capex, capacity_kg,
                    target_margin=target_margin, tax_rate=tax_rate,
                    pct_debt=pct_debt, loan_interest=loan_interest,
                    loan_term_yr=int(loan_term_yr), construction_yr=2,
                    depreciation_yr=DEPRECIATION_YR,
                    ongoing_capex_frac=ONGOING_CAPEX_FRAC)

print("10. calculate_DCF...")
dcf = calculate_DCF(opex["total_opex"], capex, capacity_kg, selling_price,
                    tax_rate=tax_rate, discount_rate=discount_rate,
                    payback_period=int(payback_period),
                    pct_debt=pct_debt, loan_interest=loan_interest,
                    loan_term_yr=int(loan_term_yr), construction_yr=2,
                    ramp_fractions=RAMP_FRACTIONS,
                    capex_yr1_frac=CAPEX_YR1_FRAC, capex_yr2_frac=CAPEX_YR2_FRAC,
                    ongoing_capex_frac=ONGOING_CAPEX_FRAC,
                    depreciation_yr=DEPRECIATION_YR)

print(f"    MSP=${MSP:.2f}/kg  IRR={dcf['IRR']:.1f}%  "
      f"NPV=${dcf['NPV']/1e6:.1f}M  TCI=${capex['TCI_total']/1e6:.0f}M  "
      f"tanks={logistics['n_tanks']}")

print("11. Checking display key accesses...")
_ = opex["feedstock"], opex["ammonia"], opex["total_opex"], opex["opex_per_kg"]
_ = capex["TCI_total"]
_ = dcf["IRR"], dcf["NPV"], dcf["cash_flows"], dcf["cum_flows"]
_ = ferm["max_OTR"], ferm["max_kla"], ferm["ferm_time"]

print("12. Test matplotlib figures...")
plt.close("all")
fig, ax = plt.subplots()
ax.plot([0, 1], [0, 1])
plt.tight_layout()
plt.close(fig)

# ── Regression baseline ──────────────────────────────────────────────────────
# Update these deliberately when a model change is expected to move the numbers.
print("13. Checking regression baseline...")
EXPECTED = {
    "MSP":         (MSP,                      3.0270),
    "IRR":         (dcf["IRR"],              28.3346),
    "NPV_M":       (dcf["NPV"] / 1e6,        11.9965),
    "TCI_M":       (capex["TCI_total"] / 1e6, 42.4617),
    "opex_per_kg": (opex["opex_per_kg"],      1.2305),
    "n_tanks":     (logistics["n_tanks"],     2),
}
failures = []
for name, (actual, expected) in EXPECTED.items():
    if abs(actual - expected) > max(abs(expected) * 0.001, 1e-9):
        failures.append(f"  {name}: got {actual:.4f}, expected {expected:.4f}")
if failures:
    print("BASELINE DRIFT:")
    print("\n".join(failures))
    raise SystemExit(1)

print("\nALL OK")
