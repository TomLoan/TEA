"""
app.py — Streamlit TEA calculator
Run with: streamlit run app.py
"""
import streamlit as st
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

from tea_functions import (
    run_chemistry, run_fermentation_model,
    calculate_plant_logistics, calculate_dsp, calculate_opex,
    size_equipment, calculate_capex,
    calculate_MSP, calculate_DCF,
    DSP_ROUTE_LIBRARY, ORGANISM_PRESETS, CARBON_SOURCE_OPTIONS,
    RAMP_FRACTIONS, CAPEX_YR1_FRAC, CAPEX_YR2_FRAC,
    ONGOING_CAPEX_FRAC, DEPRECIATION_YR,
)

st.set_page_config(page_title="Bioprocess TEA", layout="wide")
st.title("Bioprocess Techno-Economic Analysis")
st.caption("Based on Lynch 2021 FEL-1 model (±50% accuracy). All costs in 2020 USD.")

# ── Sidebar ──────────────────────────────────────────────────────────────────
with st.sidebar:
    st.header("Inputs")

    # Product
    st.subheader("Product")
    formula = st.text_input("Molecular formula", value="C7H12O4",
                             help="e.g. C7H12O4 for diethyl malonate. Supported atoms: C, H, O, N, S.")
    is_protein = st.checkbox("Protein target (add tx/tl energy cost)", value=False)
    col1, col2 = st.columns(2)
    with col1:
        avg_residue_mw = st.number_input("Avg residue MW (g/mol)", value=110.0, min_value=50.0, max_value=200.0,
                                         help="Average amino acid residue MW. Typical: 110 g/mol.",
                                         disabled=not is_protein)
    with col2:
        atp_per_residue = st.number_input("ATP per residue", value=5.0, min_value=1.0, max_value=10.0,
                                          help="Energetic cost of translation. Default 5 ATP-equiv/residue "
                                               "(4 translation + 1 mRNA). Reasonable range: 4–6.",
                                          disabled=not is_protein)

    # DSP
    st.subheader("Downstream Processing")
    dsp_route_key = st.selectbox(
        "DSP route",
        options=list(DSP_ROUTE_LIBRARY.keys()),
        format_func=lambda k: DSP_ROUTE_LIBRARY[k]['label'],
        help="Select the processing route that best matches your product. "
             "Costs scale with annual broth throughput via the 6th-tenths rule."
    )
    _route_steps = DSP_ROUTE_LIBRARY[dsp_route_key]['steps']
    with st.expander("Advanced DSP — per-step yield overrides"):
        step_overrides = []
        for _step_name, _step_default in _route_steps:
            _y = st.slider(
                _step_name, min_value=0.50, max_value=1.00,
                value=_step_default, step=0.01, format="%.2f",
                key=f"dsp_{dsp_route_key}_{_step_name}"
            )
            step_overrides.append(_y)

    # Organism
    st.subheader("Organism")
    organism = st.selectbox("Organism preset", list(ORGANISM_PRESETS.keys()))
    _org = ORGANISM_PRESETS[organism]
    st.caption(_org['note'])
    if organism == 'Mammalian (CHO-like)':
        st.warning("CHO-like cells: this model is aerobic single-substrate only. "
                   "Treat outputs as order-of-magnitude estimates.")
    _preset_media_cost = _org['media_cost']

    with st.expander("Advanced biology settings"):
        biomass_yield_override = st.number_input(
            "Biomass yield (gCDW/g carbon source)",
            value=float(_org['biomass_yield_coeff']),
            min_value=0.05, max_value=0.80, step=0.01,
            help="Cell mass produced per gram of carbon source consumed for growth. "
                 "E. coli (glucose): ~0.48. S. cerevisiae (glucose): ~0.45. "
                 "Pichia pastoris (methanol): ~0.35. Mammalian: ~0.20. "
                 "Override the preset if you have measured data."
        )
        carbon_to_co2_override = st.slider(
            "Carbon to CO₂/heat (%)",
            min_value=0, max_value=80,
            value=int(round(_org['carbon_to_co2_frac'] * 100)),
            help="Fraction of non-product glucose diverted to CO₂ and heat via overflow "
                 "metabolism or maintenance energy, rather than to biomass. "
                 "E. coli: ~20% (acetate overflow). Yeast: ~35% (Crabtree effect). "
                 "CHO: ~50% (Warburg-like lactate + high maintenance). "
                 "Genetic engineering to reduce overflow metabolism can lower this significantly."
        ) / 100.0

        st.markdown("---")
        production_mode_label = st.radio(
            "Production mode",
            ["Growth-associated", "Stationary phase"],
            help="Growth-associated: product accumulates proportionally to biomass "
                 "(logistic model, Lynch 2021 default). "
                 "Stationary phase: cells grow to a target density then production is "
                 "induced — decouples biomass from titer but lengthens batch time."
        )
        is_stationary = (production_mode_label == "Stationary phase")

        target_biomass_input = st.number_input(
            "Target biomass at induction (gCDW/L)",
            value=30.0, min_value=1.0, max_value=300.0, step=1.0,
            disabled=not is_stationary,
            help="Cell density at end of growth phase / induction point. "
                 "Typical fed-batch: 20–80 gCDW/L."
        )
        growth_time_input = st.number_input(
            "Growth phase duration (hr)",
            value=24.0, min_value=1.0, max_value=200.0, step=1.0,
            disabled=not is_stationary,
            help="Time from inoculation to induction. "
                 "Production phase duration = titer / rate."
        )

    # Fermentation
    st.subheader("Fermentation Performance")
    titer = st.number_input("Titer (g/L)", value=150.0, min_value=0.01, max_value=500.0,
                             step=0.1,
                             help="Final product concentration at harvest (g/L). Higher titer "
                                  "means fewer tanks for the same annual output, directly "
                                  "reducing capital and fixed costs per kg. The single biggest "
                                  "fermentation lever for most processes. "
                                  "Values below ~1 g/L are typical for membrane proteins and "
                                  "high-value recombinant proteins — the model handles these "
                                  "correctly but will predict high tank counts or high MSP.")
    rate = st.number_input("Rate (g/L/hr)", value=5.0, min_value=0.005, max_value=50.0,
                            step=0.01, format="%.3f",
                            help="Volumetric production rate. Growth-associated mode: average "
                                 "across the whole batch. Stationary phase mode: rate during "
                                 "the production phase only. Determines batch duration and "
                                 "therefore how many batches fit per tank per year. "
                                 "Low-expression proteins (membrane proteins, inclusion body "
                                 "products) typically 0.005–0.1 g/L/hr.")
    yield_fraction = st.slider("Yield fraction (% of theoretical)", min_value=1, max_value=99,
                                value=90,
                                help="Fraction of the thermodynamic maximum carbon yield "
                                     "actually achieved (e.g. 90 = 90% of max). Lower yield "
                                     "= more glucose consumed per kg of product. Driven by "
                                     "competing pathways, side reactions, and product "
                                     "degradation.") / 100.0

    # Plant
    st.subheader("Plant Configuration")
    capacity_kta = st.number_input("Capacity (kta)", value=15.0, min_value=0.001, max_value=500.0,
                                    step=0.001, format="%.3f",
                                    help="Target annual nameplate production in kilotonnes per year. "
                                         "Larger plants benefit from economy of scale on capital and "
                                         "labour, but require more upfront investment. This model is "
                                         "most reliable in the 1–100 kta range. Below 0.1 kta "
                                         "(100 t/yr) the capital cost correlations are less reliable "
                                         "but outputs remain useful for FEL-1 order-of-magnitude "
                                         "estimates of specialty proteins and high-value products.")
    tank_volume_L = st.selectbox("Tank volume (L)", [250_000, 500_000, 1_000_000], index=1,
                                  help="Total fermentation vessel volume. Working volume is 85% of "
                                       "total. Larger tanks reduce unit count and simplify operations "
                                       "but cost more per unit. 500,000 L is typical for large-scale "
                                       "commodity fermentation.")
    annual_uptime = st.slider("Annual uptime (%)", min_value=50, max_value=99, value=90,
                               help="Fraction of the year the plant is in productive operation, "
                                    "accounting for planned maintenance shutdowns, cleaning cycles, "
                                    "and unplanned downtime. 90% ≈ 7,884 hr/yr.") / 100.0
    ferm_temp_C = st.number_input("Fermentation temp (°C)", value=37.0, min_value=4.0, max_value=70.0,
                                   help="Fermentation temperature affects cooling water demand: higher "
                                        "temperature increases the temperature differential to cooling "
                                        "water (~4°C supply), reducing the flow rate needed to remove "
                                        "the same heat load. Also affects sterilisation energy.")
    turnaround_time = st.number_input("Turnaround time (hr)", value=16.0, min_value=1.0, max_value=72.0,
                                       help="Time between the end of one batch and the start of the "
                                            "next — draining, CIP cleaning, sterilisation, and filling. "
                                            "Directly limits how many batches fit per tank per year. "
                                            "Typical range: 12–24 hr for large vessels.")
    batches_on_spec = st.slider("Batches on-spec (%)", min_value=50, max_value=100, value=95,
                                 help="Fraction of batches that meet product quality specification "
                                      "and are sold. Failed batches consume raw materials and "
                                      "utilities but generate no revenue. 95% is a reasonable "
                                      "mature-process assumption.") / 100.0

    # Raw material prices
    st.subheader("Raw Material Prices")
    carbon_source = st.selectbox(
        "Carbon source",
        options=list(CARBON_SOURCE_OPTIONS.keys()),
        index=0,
        help="Primary fermentation feedstock. Choose 'Methanol' together with the "
             "'Pichia pastoris (methanol)' organism preset. Stoichiometry uses a "
             "glucose-equivalent basis (carbon content per gram differs by <7%).",
    )
    _cs = CARBON_SOURCE_OPTIONS[carbon_source]
    price_carbon_per_kg = st.number_input(
        _cs['label'],
        value=_cs['default_per_kg'],
        min_value=0.01, max_value=10.0,
        help=_cs['note'] + " The dominant raw material cost for most processes — "
             "see the Sensitivity tab to quantify the impact.",
    )
    price_ammonia_per_kg = st.number_input("Ammonia ($/kg)", value=0.26, min_value=0.01, max_value=10.0,
                                            help="Anhydrous ammonia or equivalent nitrogen source. "
                                                 "Used for biomass growth (all products) and as a "
                                                 "stoichiometric reactant for N-containing products. "
                                                 "Industrial price: ~$0.22–0.44/kg.")
    price_mgso4_per_kg = st.number_input("MgSO4 ($/kg)", value=0.30, min_value=0.0, max_value=10.0,
                                          help="Industrial-grade magnesium sulfate — the sulfur source "
                                               "for products containing S atoms. Contributes negligibly "
                                               "to cost for non-S products; set to 0 if not applicable.")
    media_cost_per_kgCDW = st.number_input("Media cost ($/kgCDW)", value=_preset_media_cost,
                                            min_value=0.0, max_value=50.0,
                                            help="Mineral salts and micronutrients cost per kg of cell "
                                                 "dry weight produced. Default $0.40/kgCDW (Lynch 2021) "
                                                 "for simple mineral medium. Typical range: "
                                                 "$0.25–$0.80/kgCDW for bacteria; $5+/kgCDW for "
                                                 "mammalian cells requiring complex media.")
    price_NaOH_per_kg = st.number_input("NaOH ($/kg)", value=0.15, min_value=0.01, max_value=5.0,
                                         help="Caustic soda used for pH control during fermentation "
                                              "and CIP cleaning between batches. Industrial grade: "
                                              "~$0.10–0.25/kg.")
    price_peracetic_per_L = st.number_input("Peracetic acid ($/L)", value=5.00, min_value=0.1, max_value=50.0,
                                             help="Peracetic acid solution used for CIP sterilisation "
                                                  "of tanks and pipework between batches. Typical "
                                                  "industrial price: $3–8/L for ~15% solution.")

    # Utility prices
    st.subheader("Utility Prices")
    price_electricity = st.number_input("Electricity ($/kWh)", value=0.11, min_value=0.01, max_value=1.0,
                                         help="Grid electricity cost. Drives aeration (mass transfer "
                                              "blowers), agitation, and cooling system pumps. US "
                                              "industrial average: ~$0.07–0.12/kWh. Higher electricity "
                                              "cost significantly impacts O₂-intensive Case 3 products.")
    price_natural_gas = st.number_input("Natural gas ($/MMBtu)", value=3.11, min_value=0.1, max_value=30.0,
                                         help="Natural gas for steam generation — used for sterilisation "
                                              "of fermentation media and biomass heat-kill. US industrial "
                                              "price: ~$2–6/MMBtu. Henry Hub spot price as of 2020: "
                                              "$2.03/MMBtu (Lynch 2021 default: $3.11/MMBtu).")
    CEPCI = st.number_input("CEPCI", value=603, min_value=100, max_value=1500,
                             help="Chemical Engineering Plant Cost Index — used to scale equipment "
                                  "purchase costs to current year. 2020 = 603 (paper baseline). "
                                  "2024 ≈ 800 (+33%). Update this to get current-dollar estimates. "
                                  "Source: Chemical Engineering magazine.")

    # Financial
    st.subheader("Financial Parameters")
    selling_price = st.number_input("Selling price ($/kg)", value=2.50, min_value=0.01, max_value=100_000.0,
                                     help="Your product's target or estimated market price. Used to "
                                          "calculate revenue, NPV, and IRR in the Financials tab. "
                                          "Does not affect MSP — MSP is the break-even price at the "
                                          "target margin, independent of what you choose to charge. "
                                          "Specialty proteins and pharmaceutical-grade products can "
                                          "range from $1,000 to $100,000+/kg.")
    target_margin = st.slider("Target margin (%)", min_value=0, max_value=80, value=30,
                               help="Net profit margin target used to calculate MSP. E.g. 30% means "
                                    "30 cents of every revenue dollar is net income after tax, "
                                    "depreciation, and interest. Higher margin = higher MSP. "
                                    "Setting 0% gives the true break-even price.") / 100.0
    discount_rate = st.slider("Discount rate (%)", min_value=1, max_value=50, value=20,
                               help="Hurdle rate — the minimum acceptable return on invested capital. "
                                    "Used to discount future cash flows when calculating NPV. "
                                    "20% is typical for early-stage bioprocess ventures (Lynch 2021); "
                                    "established commodity chemical plants may use 10–12%.") / 100.0
    tax_rate = st.slider("Tax rate (%)", min_value=0, max_value=50, value=21,
                          help="Corporate income tax rate applied to pre-tax earnings. "
                               "US federal statutory rate: 21%. Combined federal + state "
                               "effective rates typically 25–28%.") / 100.0
    payback_period = st.number_input("Payback period (years)", value=20, min_value=5, max_value=40,
                                      help="Total horizon over which the DCF model runs, including "
                                           "construction years. Does not affect MSP. Longer periods "
                                           "improve NPV and IRR by capturing more production years, "
                                           "but increase uncertainty.")
    pct_debt = st.slider("Debt fraction (%)", min_value=0, max_value=90, value=60,
                          help="Fraction of Total Capital Investment financed by debt (bank loan or "
                               "bonds). The remainder is equity. 60% debt / 40% equity is a common "
                               "assumption for mid-scale industrial bioprocess projects.") / 100.0
    loan_interest = st.slider("Loan interest (%)", min_value=1, max_value=20, value=8,
                               help="Annual interest rate on the debt portion of capital. Affects "
                                    "interest expense in the P&L and cash flows. Typical project "
                                    "finance range: 6–12% depending on credit quality and market "
                                    "conditions.") / 100.0
    loan_term_yr = st.number_input("Loan term (years)", value=10, min_value=1, max_value=30,
                                    help="Number of years over which the debt is repaid via equal "
                                         "annual payments (annuity). Longer terms reduce annual "
                                         "principal payments but increase total interest paid.")

# ── Derived inputs ────────────────────────────────────────────────────────────
price_carbon_per_g  = price_carbon_per_kg / 1000
price_ammonia_per_g = price_ammonia_per_kg / 1000
cost_of_fuel = price_natural_gas / 1.05505
capacity_kg = capacity_kta * 1e6

OPEX_KWARGS = dict(
    price_glucose_per_g=price_carbon_per_g,
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

# ── Main calculation (no Run button — reactive) ───────────────────────────────
try:
    # S-warning
    chem = run_chemistry(formula, is_protein=is_protein,
                         avg_residue_mw=avg_residue_mw, atp_per_residue=atp_per_residue)
    if chem['atoms']['S'] > 0:
        st.info("S detected in formula. Verify the MgSO₄ price and DSP route are appropriate "
                "for a sulfur-containing product.")

    ferm = run_fermentation_model(
        titer, rate, yield_fraction, chem,
        biomass_yield_coeff=biomass_yield_override,
        carbon_to_co2_frac=carbon_to_co2_override,
        production_mode='stationary_phase' if is_stationary else 'growth_associated',
        target_biomass=target_biomass_input if is_stationary else None,
        growth_time_hr=growth_time_input    if is_stationary else None,
    )
    logistics = calculate_plant_logistics(
        capacity_kta, annual_uptime, batches_on_spec,
        tank_volume_L, turnaround_time, ferm
    )

    if rate >= titer:
        st.warning(f"Rate ({rate} g/L/hr) ≥ titer ({titer} g/L): fermentation time is "
                   f"{ferm['ferm_time']:.1f} hr — likely unrealistic.")
    if logistics['n_tanks'] > 30:
        st.warning(f"{logistics['n_tanks']} tanks required. Consider increasing titer, "
                   f"tank volume, or reducing capacity target.")
    if CEPCI < 700:
        st.info(f"CEPCI = {CEPCI} (2020 baseline). Current value (2024) ≈ 800. "
                "Equipment cost estimates may be understated by ~30%.")

    # DSP: compute route CAPEX/OPEX (scales with broth throughput)
    annual_broth_m3 = logistics['annual_ferm_vol'] / 1000
    dsp = calculate_dsp(dsp_route_key, annual_broth_m3,
                        step_yield_overrides=step_overrides)

    # Two-pass OPEX (first pass without other_fixed to get CAPEX, then re-run)
    opex_pass1 = calculate_opex(logistics, ferm, chem, dsp=dsp, **OPEX_KWARGS)
    sizing = size_equipment(logistics, ferm, opex_pass1, tank_volume_L, ferm_temp_C)
    capex = calculate_capex(sizing, dsp=dsp)
    other_fixed = 0.037 * capex['TCI_total']
    opex = calculate_opex(logistics, ferm, chem, dsp=dsp,
                          other_fixed_costs=other_fixed, **OPEX_KWARGS)

    MSP = calculate_MSP(opex['total_opex'], capex, capacity_kg,
                        target_margin=target_margin, tax_rate=tax_rate,
                        pct_debt=pct_debt, loan_interest=loan_interest,
                        loan_term_yr=int(loan_term_yr), construction_yr=2,
                        depreciation_yr=DEPRECIATION_YR,
                        ongoing_capex_frac=ONGOING_CAPEX_FRAC)
    dcf = calculate_DCF(opex['total_opex'], capex, capacity_kg, selling_price,
                        tax_rate=tax_rate, discount_rate=discount_rate,
                        payback_period=int(payback_period),
                        pct_debt=pct_debt, loan_interest=loan_interest,
                        loan_term_yr=int(loan_term_yr), construction_yr=2,
                        ramp_fractions=RAMP_FRACTIONS,
                        capex_yr1_frac=CAPEX_YR1_FRAC, capex_yr2_frac=CAPEX_YR2_FRAC,
                        ongoing_capex_frac=ONGOING_CAPEX_FRAC,
                        depreciation_yr=DEPRECIATION_YR)

    error_msg = None

except ValueError as e:
    error_msg = str(e)
except Exception as e:
    error_msg = f"Unexpected error: {e}"

# ── Output area ───────────────────────────────────────────────────────────────
if error_msg:
    st.error(f"Calculation error: {error_msg}")
    st.stop()

if MSP == float('inf'):
    st.error("MSP is undefined: (1 − tax rate) − target margin ≤ 0. "
             "Reduce target margin or tax rate.")
    st.stop()

# Headline metrics
mcol1, mcol2, mcol3, mcol4 = st.columns(4)
mcol1.metric("MSP", f"${MSP:.2f}/kg",
    help="Minimum Selling Price — the product price at which the project exactly meets the "
         "target profit margin, assuming full nameplate production. The primary output for "
         "R&D goal-setting: if the current MSP exceeds your market price target, use the "
         "Sensitivity tab to identify which parameters to improve.")
mcol2.metric("IRR", f"{dcf['IRR']:.1f}%",
    help="Internal Rate of Return — the annualised return on invested capital over the "
         "payback period. Compared against the hurdle rate (discount rate): IRR > hurdle "
         "rate means the project creates value. Does not depend on selling price — it is "
         "evaluated at the MSP.")
mcol3.metric("NPV (20-yr)", f"${dcf['NPV']/1e6:.1f}M",
    help="Net Present Value over the payback period, discounted at the hurdle rate. Positive "
         "NPV = the project returns more than the cost of capital. Calculated at the user-set "
         "selling price (not MSP) — use the Financials tab to see the full cash flow profile.")
mcol4.metric("TCI (total)", f"${capex['TCI_total']/1e6:.0f}M",
    help="Total Capital Investment — all capital required to build the plant, including "
         "fermentation equipment, utilities, site development, indirect costs (engineering, "
         "contingency, startup), working capital, and the DSP section. See the CAPEX tab "
         "for a full breakdown.")

with st.expander("Model scope & assumptions"):
    st.markdown("""
**Applicable to:** Aerobic batch/fed-batch fermentation, glucose or methanol carbon source,
products containing C/H/O/N/S atoms, scale ~1–100 kta. **Accuracy: ±50% (FEL-1 level).**

**Methanol (Pichia pastoris):** Supported via the organism preset and carbon source selector.
Stoichiometry uses a glucose-equivalent basis (carbon content per gram differs by <7%).

**Not modelled:** Anaerobic fermentation · Other alternative feedstocks (xylose, glycerol, acetate) ·
GMP/regulatory costs · Multi-substrate media · Fed-batch feeding strategies.

| Hard-coded assumption | Value | Source |
|---|---|---|
| Seed train cost | 27% of main fermentation area | Lynch 2021 |
| Labour | 2.5 FTE/tank × $60k loaded × 1.5× overhead | Lynch 2021 |
| DO setpoint | 25% of air saturation | Typical aerobic |
| Tank working volume | 85% of total volume | Engineering rule of thumb |
| Depreciation | 10-yr straight-line | Lynch 2021 |
| Ongoing maintenance capex | 10% TCI/yr | Lynch 2021 |
| Ramp-up schedule | 50% → 75% → 100% over 3 years | Lynch 2021 |
| Inoculum | 1% of final biomass | Typical practice |

*Reference: Lynch et al. 2021, J. Cleaner Production*
""")

st.divider()

# Tabs
tab_chem, tab_ferm, tab_opex, tab_capex, tab_fin, tab_sens = st.tabs(
    ["Chemistry", "Fermentation", "OPEX", "CAPEX", "Financials", "Sensitivity"]
)

# ── Chemistry tab ─────────────────────────────────────────────────────────────
with tab_chem:
    st.subheader("Stoichiometry & Yields")
    eq = chem['equation']
    ccol1, ccol2, ccol3 = st.columns(3)
    ccol1.metric("Molecular weight", f"{chem['MW']:.2f} g/mol",
        help="Molecular weight of the target product in g/mol, calculated from the formula.")
    ccol2.metric("H₂:CO₂ ratio", f"{eq['ratio']:.3f}",
        help="Ratio of reducing equivalents to carbon in the product, relative to glucose "
             "(ratio = 2). Above 2: product is more reduced than glucose — some carbon "
             "must leave as CO₂ (Case 2). Below 2: product is more oxidised — O₂ is "
             "consumed as a stoichiometric reactant (Case 3). Equal to 2: neutral (Case 1).")
    ccol3.metric("Redox class", f"Case {eq['case']}",
        help="Case 1 (neutral): no O₂ consumed or CO₂ byproduct from the product reaction. "
             "Case 2 (reduced): product more reduced than glucose — CO₂ is a byproduct. "
             "Case 3 (oxidised): product more oxidised than glucose — O₂ is a reactant, "
             "increasing aeration cost significantly.")

    ccol4, ccol5 = st.columns(2)
    ccol4.metric("Theoretical yield", f"{chem['theoretical_yield']:.4f} g/g glucose",
        help="Maximum possible yield from atom balance alone (g product per g glucose), "
             "assuming 100% of carbon goes to product with no biomass or waste. The actual "
             "yield achieved is this × yield fraction × (1 - overflow fraction).")
    if carbon_source == 'Methanol':
        st.info(
            "**Methanol selected.** Stoichiometry above is shown on a glucose basis. "
            "Glucose and methanol have nearly identical carbon content per gram "
            "(40.0% vs 37.5%), so the glucose-based theoretical yield approximates "
            "the methanol-based value to within ~7% — within FEL-1 (±50%) tolerance. "
            "Set your yield fraction and titer/rate based on actual methanol fermentation data."
        )
    if is_protein and chem['txl_glucose_g_per_g'] > 0:
        eff = 1.0 / (1.0 / chem['theoretical_yield'] + chem['txl_glucose_g_per_g'])
        ccol5.metric("Tx/tl glucose overhead", f"{chem['txl_glucose_g_per_g']:.4f} g/g protein",
                     delta=f"Effective yield: {eff:.4f} g/g", delta_color="inverse")

    # Stoichiometry table
    st.markdown("**Balanced equation (per mol glucose):**")
    rows = [
        {"Species": "Glucose", "Role": "Reactant", "Moles": f"{eq['glucose']:.4f}"},
        {"Species": "NH₃", "Role": "Reactant", "Moles": f"{eq['NH3']:.4f}"},
    ]
    if eq['H2SO4'] > 1e-9:
        rows.append({"Species": "H₂SO₄", "Role": "Reactant", "Moles": f"{eq['H2SO4']:.4f}"})
    if eq['O2'] > 1e-9:
        rows.append({"Species": "O₂", "Role": "Reactant", "Moles": f"{eq['O2']:.4f}"})
    rows.append({"Species": formula, "Role": "Product", "Moles": f"{eq['product']:.4f}"})
    rows.append({"Species": "H₂O", "Role": "Product", "Moles": f"{eq['H2O']:.4f}"})
    if eq['CO2'] > 1e-9:
        rows.append({"Species": "CO₂", "Role": "Byproduct", "Moles": f"{eq['CO2']:.4f}"})
    st.dataframe(rows, use_container_width=True)

# ── Fermentation tab ──────────────────────────────────────────────────────────
with tab_ferm:
    st.subheader("Fermentation Model")
    fcol1, fcol2, fcol3, fcol4 = st.columns(4)
    fcol1.metric("Fermentation time", f"{ferm['ferm_time']:.1f} hr",
        help="Total batch duration = titer ÷ rate (growth-associated), or growth phase + "
             "production phase (stationary phase). Determines how many batches fit per tank "
             "per year and therefore the number of tanks required.")
    fcol2.metric("Final biomass", f"{ferm['final_biomass']:.1f} g/L CDW",
        help="Predicted cell dry weight at end of fermentation (g/L). Drives oxygen demand, "
             "media cost, and the biomass heat-kill step. In stationary-phase mode this is "
             "the user-set target biomass at induction.")
    fcol3.metric("Overall yield", f"{ferm['overall_yield']:.4f} g/g glucose",
        help="Total product mass per unit of glucose consumed (g product / g glucose), "
             "accounting for all glucose used by product formation, biomass growth, and "
             "overflow metabolism. Determines how much glucose the plant buys per year.")
    fcol4.metric("Cumulative O₂", f"{ferm['cumulative_O2']:.1f} mmol/L",
        help="Total oxygen consumed per litre of broth per batch (mmol O₂/L). The primary "
             "driver of annual compressed air cost and cooling water demand — both scale "
             "directly with this number.")

    fcol5, fcol6, fcol7, fcol8 = st.columns(4)
    fcol5.metric("Max OTR", f"{ferm['max_OTR']:.2f} mmol/L/hr",
        help="Peak oxygen transfer rate at the logistic growth midpoint (mmol O₂/L/hr). "
             "Display metric only — does not affect any cost calculation. Note: the published "
             "equations yield a value ~14× below Lynch Fig. 5b, likely due to undescribed "
             "terms in the original tool; cumulative O₂ (which drives costs) is unaffected.")
    fcol6.metric("Max kLa", f"{ferm['max_kla']:.4f} s⁻¹",
        help="Minimum volumetric mass-transfer coefficient the bioreactor must achieve to "
             "supply oxygen at the growth peak (s⁻¹). Useful for bioreactor design scoping. "
             "Display metric only — does not affect cost calculations.")
    fcol7.metric("Max cooling", f"{ferm['max_cooling_rate']:.1f} kJ/L/hr",
        help="Peak heat removal rate at the logistic growth midpoint (kJ/L/hr), derived from "
             "max OTR × 0.46 kJ/mmol O₂ (Doran 1995). Display metric only — equipment and "
             "costs are sized from cumulative O₂, not from the instantaneous peak.")
    fcol8.metric("Tanks required", f"{logistics['n_tanks']}",
        help="Number of fermentation tanks needed to meet the annual capacity target at this "
             "configuration. Drives capital cost (Area 1) and labour. Reduce by increasing "
             "titer, tank volume, or annual uptime.")

    # Phase duration metrics (stationary phase only)
    if ferm['production_mode'] == 'stationary_phase':
        pcol1, pcol2, pcol3 = st.columns(3)
        pcol1.metric("Growth phase",     f"{ferm['t_growth']:.1f} hr")
        pcol2.metric("Production phase", f"{ferm['t_production']:.1f} hr")
        pcol3.metric("Total batch",      f"{ferm['ferm_time']:.1f} hr")

    # Sugar partitioning
    st.markdown("**Glucose partitioning (g/L):**")
    if ferm['production_mode'] == 'stationary_phase':
        stoich_losses = ferm['stoich_total'] - ferm['sugar_to_product']
        sugar_rows = [
            {"Category": "→ Product (stoichiometric)",            "g/L": f"{ferm['sugar_to_product']:.2f}"},
            {"Category": "→ Yield losses (product pathway)",      "g/L": f"{stoich_losses:.2f}"},
            {"Category": "→ Biomass growth (independent budget)", "g/L": f"{ferm['glucose_to_biomass']:.2f}"},
            {"Category": "→ Overflow / maintenance (CO₂ + heat)", "g/L": f"{ferm['glucose_to_co2']:.2f}"},
        ]
    else:
        sugar_rows = [
            {"Category": "→ Product (stoichiometric)",            "g/L": f"{ferm['sugar_to_product']:.2f}"},
            {"Category": "→ Biomass formation",                   "g/L": f"{ferm['glucose_to_biomass']:.2f}"},
            {"Category": "→ Overflow / maintenance (CO₂ + heat)", "g/L": f"{ferm['glucose_to_co2']:.2f}"},
        ]
    if ferm['sugar_for_txl'] > 0:
        sugar_rows.append({"Category": "→ Tx/tl energy (protein only)", "g/L": f"{ferm['sugar_for_txl']:.2f}"})
    sugar_rows.append({"Category": "TOTAL glucose consumed", "g/L": f"{ferm['total_sugar']:.2f}"})
    st.dataframe(sugar_rows, use_container_width=True)

    # Time course plot
    st.markdown("**Fermentation time course:**")
    fig_ferm, ax1 = plt.subplots(figsize=(8, 4))
    ax2 = ax1.twinx()
    ax1.plot(ferm['t_points'], ferm['biomass_curve'], 'b-', linewidth=2, label='Biomass (gCDW/L)')
    ax2.plot(ferm['t_points'], ferm['product_curve'], 'r--', linewidth=2, label=f'{formula} (g/L)')
    ax1.set_xlabel('Time (hr)')
    ax1.set_ylabel('Biomass (gCDW/L)', color='blue')
    ax2.set_ylabel(f'{formula} (g/L)', color='red')
    ax1.tick_params(axis='y', labelcolor='blue')
    ax2.tick_params(axis='y', labelcolor='red')
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper left')
    plt.tight_layout()
    st.pyplot(fig_ferm)
    plt.close(fig_ferm)

# ── OPEX tab ──────────────────────────────────────────────────────────────────
with tab_opex:
    st.subheader("Operating Cost Breakdown")
    ocol1, ocol2, ocol3 = st.columns(3)
    ocol1.metric("Total OPEX", f"${opex['total_opex']/1e6:.2f}M/yr",
        help="Total annual operating cost at nameplate capacity: fermentation (raw materials, "
             "utilities, labour, other fixed costs) plus DSP operating cost.")
    ocol2.metric("OPEX per kg", f"${opex['opex_per_kg']:.3f}/kg",
        help="Specific OPEX — total annual operating cost divided by annual nameplate "
             "production. The main driver of MSP alongside capital charges.")
    ocol3.metric("Fermentation OPEX", f"${opex['ferm_opex']/1e6:.2f}M/yr",
        help="OPEX from the fermentation section only (raw materials, utilities, labour, "
             "and other fixed costs), before adding DSP operating costs.")

    opex_items = [
        ("Glucose",          opex['glucose']),
        ("Ammonia",          opex['ammonia']),
        ("Sulfate (MgSO4)",  opex['sulfate']),
        ("Media (salts)",    opex['media']),
        ("CIP chemicals",    opex['CIP']),
        ("Water",            opex['water']),
        ("Compressed air",   opex['compressed_air']),
        ("Mass transfer",    opex['mass_transfer']),
        ("Cooling water",    opex['cooling_water']),
        ("Sterilisation",    opex['sterilisation']),
        ("Heat kill",        opex['heat_kill']),
        ("Centrifugation",   opex['centrifugation']),
        ("Labour",           opex['labour']),
        ("Other fixed",      opex['other_fixed']),
        ("DSP (allocated)",  opex['DSP_opex']),
    ]
    opex_table = [
        {"Item": name, "$/yr": f"${cost:,.0f}", "$/kg": f"${cost/capacity_kg:.3f}"}
        for name, cost in opex_items if cost > 0
    ]
    st.dataframe(opex_table, use_container_width=True)

    # Bar chart
    labels = [r[0] for r in opex_items if r[1] > 0]
    values = [r[1]/1e6 for r in opex_items if r[1] > 0]
    fig_opex, ax = plt.subplots(figsize=(8, 4))
    ax.barh(labels, values, color='steelblue')
    ax.set_xlabel("$/yr (M)")
    ax.set_title("OPEX breakdown")
    plt.tight_layout()
    st.pyplot(fig_opex)
    plt.close(fig_opex)

    st.divider()
    st.subheader("DSP Recovery Funnel")
    st.caption(f"Route: {dsp['route_label']}")
    dsp_col1, dsp_col2, dsp_col3 = st.columns(3)
    dsp_col1.metric("Overall DSP recovery", f"{dsp['overall_yield']*100:.1f}%",
        help="Product of all per-step yield fractions — the fraction of fermented product "
             "that reaches final specification. Affects raw material demand (more glucose "
             "needed to compensate for losses) and DSP cost.")
    dsp_col2.metric("DSP CAPEX", f"${dsp['dsp_capex']/1e6:.1f}M",
        help="Capital cost of the downstream processing section, scaled from a route-specific "
             "reference cost using annual broth throughput and the 6th-tenths rule. Larger "
             "plants get a lower specific DSP capital cost per tonne.")
    dsp_col3.metric("DSP OPEX", f"${dsp['dsp_opex']/1e6:.2f}M/yr",
        help="Annual DSP operating cost — covers utilities, consumables (resins, membranes, "
             "solvents), and allocated labour. Estimated as a route-specific fraction of DSP "
             "CAPEX per year.")

    # Cumulative recovery at each step
    _cum = [1.0]
    for _y in dsp['step_yields']:
        _cum.append(_cum[-1] * _y)
    _n_steps = len(dsp['step_names'])
    fig_dsp, ax_dsp = plt.subplots(figsize=(7, max(2.2, _n_steps * 0.6 + 0.8)))
    _y_pos = list(range(_n_steps))
    _bars = ax_dsp.barh(_y_pos, [c * 100 for c in _cum[1:]], color='steelblue', height=0.5)
    ax_dsp.set_xlim(0, 105)
    ax_dsp.set_yticks(_y_pos)
    ax_dsp.set_yticklabels(dsp['step_names'])
    ax_dsp.invert_yaxis()
    ax_dsp.set_xlabel("Cumulative recovery (%)")
    ax_dsp.set_title(f"Overall DSP recovery: {dsp['overall_yield']*100:.1f}%")
    for _bar, _val in zip(_bars, _cum[1:]):
        ax_dsp.text(_bar.get_width() + 0.5, _bar.get_y() + _bar.get_height() / 2,
                    f"{_val*100:.1f}%", va='center', fontsize=9)
    plt.tight_layout()
    st.pyplot(fig_dsp)
    plt.close(fig_dsp)

# ── CAPEX tab ─────────────────────────────────────────────────────────────────
with tab_capex:
    st.subheader("Capital Cost Breakdown")
    kacol1, kacol2, kacol3 = st.columns(3)
    kacol1.metric("TCI (upstream)", f"${capex['TCI']/1e6:.1f}M",
        help="Total Capital Investment for the fermentation + utilities sections only "
             "(Areas 1–4): equipment, installation, site development, indirect costs "
             "(engineering, contingency, startup), and working capital. Excludes DSP.")
    kacol2.metric("TCI (incl. DSP)", f"${capex['TCI_total']/1e6:.1f}M",
        help="Total Capital Investment including both fermentation and DSP sections. "
             "This is the number used to calculate MSP, IRR, and NPV.")
    kacol3.metric("TCI per kg/yr", f"${capex['TCI_total']/capacity_kg:.2f}/kg",
        help="Capital intensity — total TCI divided by annual nameplate production. "
             "Useful for benchmarking against published process economics. Typical range: "
             "$1–10/kg for commodity bioprocesses; $50–500+/kg for therapeutics.")

    capex_area_items = [
        ("Area 1: Main fermentation",      capex['area1']),
        ("Area 2: Seed fermentation (27%)", capex['area2']),
        ("Area 3: Cell removal",           capex['area3']),
        ("Area 4: Process utilities",      capex['area4']),
        ("DSP capital",                    capex['DSP_capex']),
    ]
    capex_struct_items = [
        ("TIC (installed)",   capex['TIC_upstream']),
        ("Site dev + admin",  capex['TDC'] - capex['TIC_upstream']),
        ("Indirect costs",    capex['indirect']),
        ("Working capital",   capex['WC']),
        ("DSP capital",       capex['DSP_capex']),
    ]

    kacol4, kacol5 = st.columns(2)
    with kacol4:
        st.markdown("**By area:**")
        area_table = [{"Area": name, "$M": f"${cost/1e6:.1f}M"}
                      for name, cost in capex_area_items]
        st.dataframe(area_table, use_container_width=True)
    with kacol5:
        st.markdown("**Capital structure:**")
        struct_table = [{"Item": name, "$M": f"${cost/1e6:.1f}M"}
                        for name, cost in capex_struct_items]
        st.dataframe(struct_table, use_container_width=True)

    # Bar chart
    fig_capex, ax = plt.subplots(figsize=(7, 3))
    ax.bar([n for n, _ in capex_area_items],
           [v/1e6 for _, v in capex_area_items], color='darkorange')
    ax.set_ylabel("$M")
    ax.set_title("CAPEX by area")
    plt.xticks(rotation=30, ha='right')
    plt.tight_layout()
    st.pyplot(fig_capex)
    plt.close(fig_capex)

# ── Financials tab ────────────────────────────────────────────────────────────
with tab_fin:
    st.subheader("DCF Financial Analysis")
    fincol1, fincol2, fincol3, fincol4 = st.columns(4)
    fincol1.metric("MSP", f"${MSP:.2f}/kg")
    fincol2.metric("IRR", f"{dcf['IRR']:.1f}%",
                   delta=f"vs hurdle {discount_rate*100:.0f}%",
                   delta_color="normal" if dcf['IRR'] > discount_rate * 100 else "inverse")
    fincol3.metric("NPV (20-yr)", f"${dcf['NPV']/1e6:.1f}M",
                   delta_color="normal" if dcf['NPV'] > 0 else "inverse")
    fincol4.metric("ROI", f"{dcf['ROI']:.1f}%")

    # Cash flow table
    years = list(range(dcf['total_yrs']))
    cf_table = []
    for yr in years:
        label = (f"Yr {yr+1} (const.)" if yr < dcf['construction_yr']
                 else f"Yr {yr+1} (prod.)")
        cf_table.append({
            "Year": label,
            "Revenue ($M)": f"${dcf['revenues'][yr]/1e6:.2f}",
            "Net Income ($M)": f"${dcf['net_incomes'][yr]/1e6:.2f}",
            "Cash Flow ($M)": f"${dcf['cash_flows'][yr]/1e6:.2f}",
            "Cumulative ($M)": f"${dcf['cum_flows'][yr]/1e6:.2f}",
        })
    st.dataframe(cf_table, use_container_width=True)

    # Cumulative cash flow plot
    fig_dcf, ax = plt.subplots(figsize=(9, 4))
    yr_labels = list(range(1, dcf['total_yrs'] + 1))
    cum = [v / 1e6 for v in dcf['cum_flows']]
    ax.bar(yr_labels, [v / 1e6 for v in dcf['cash_flows']],
           color=['#d62728' if v < 0 else '#2ca02c' for v in dcf['cash_flows']],
           alpha=0.6, label='Annual cash flow')
    ax.plot(yr_labels, cum, 'k-o', markersize=4, linewidth=1.5, label='Cumulative')
    ax.axhline(0, color='black', linewidth=0.8)
    ax.axvline(dcf['construction_yr'] + 0.5, color='grey', linestyle='--', linewidth=1,
               label='Production start')
    ax.set_xlabel('Year')
    ax.set_ylabel('$M')
    ax.set_title('Cash flow proforma')
    ax.legend()
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'${x:.0f}M'))
    plt.tight_layout()
    st.pyplot(fig_dcf)
    plt.close(fig_dcf)

    # MSP vs selling price comparison
    if selling_price > 0:
        delta_pct = (selling_price - MSP) / MSP * 100
        label = "above MSP" if selling_price >= MSP else "below MSP"
        st.info(f"Selling price ${selling_price:.2f}/kg is "
                f"{abs(delta_pct):.1f}% {label} (MSP ${MSP:.2f}/kg)")

# ── Sensitivity tab ───────────────────────────────────────────────────────────
with tab_sens:
    st.subheader("One-at-a-time Sensitivity Analysis")
    st.caption(
        "Each parameter is varied independently over a plausible range while all "
        "others stay at their current values. Longer bars = higher economic leverage. "
        "Parameters are ranked by their total MSP swing (best case → worst case)."
    )

    target_msp = st.number_input(
        "Target MSP ($/kg)", value=round(MSP * 0.7, 2), min_value=0.01,
        help="Set a cost target. The table below shows what each parameter would need "
             "to reach alone, and whether it is achievable within the modelled range."
    )

    # ── Helper: rerun full pipeline with one parameter overridden ──────────────
    def _msp(**ov):
        _titer         = ov.get('titer',              titer)
        _rate          = ov.get('rate',               rate)
        _yield_frac    = ov.get('yield_fraction',     yield_fraction)
        _dsp_yield_ov  = ov.get('dsp_yield_override', None)
        _co2_frac      = ov.get('carbon_to_co2_frac', carbon_to_co2_override)
        _cs_per_kg     = ov.get('price_carbon_per_kg', price_carbon_per_kg)
        _cap_kta       = ov.get('capacity_kta',       capacity_kta)
        _cap_kg        = _cap_kta * 1e6
        _opex_kw       = dict(OPEX_KWARGS, price_glucose_per_g=_cs_per_kg / 1000)
        try:
            _ferm = run_fermentation_model(
                _titer, _rate, _yield_frac, chem,
                biomass_yield_coeff=biomass_yield_override,
                carbon_to_co2_frac=_co2_frac,
                production_mode='stationary_phase' if is_stationary else 'growth_associated',
                target_biomass=target_biomass_input if is_stationary else None,
                growth_time_hr=growth_time_input    if is_stationary else None,
            )
            _log = calculate_plant_logistics(
                _cap_kta, annual_uptime, batches_on_spec,
                tank_volume_L, turnaround_time, _ferm
            )
            _broth_m3 = _log['annual_ferm_vol'] / 1000
            if _dsp_yield_ov is not None:
                # scale all step yields proportionally to reach the target overall yield
                _base = max(dsp['overall_yield'], 1e-9)
                _scaled = [min(1.0, y * _dsp_yield_ov / _base) for y in dsp['step_yields']]
                _dsp = calculate_dsp(dsp_route_key, _broth_m3, step_yield_overrides=_scaled)
            else:
                _dsp = calculate_dsp(dsp_route_key, _broth_m3, step_yield_overrides=step_overrides)
            _op1   = calculate_opex(_log, _ferm, chem, dsp=_dsp, **_opex_kw)
            _sz    = size_equipment(_log, _ferm, _op1, tank_volume_L, ferm_temp_C)
            _cx    = calculate_capex(_sz, dsp=_dsp)
            _op    = calculate_opex(_log, _ferm, chem, dsp=_dsp,
                         other_fixed_costs=0.037 * _cx['TCI_total'], **_opex_kw)
            _m     = calculate_MSP(_op['total_opex'], _cx, _cap_kg,
                         target_margin=target_margin, tax_rate=tax_rate,
                         pct_debt=pct_debt, loan_interest=loan_interest,
                         loan_term_yr=int(loan_term_yr), construction_yr=2,
                         depreciation_yr=DEPRECIATION_YR, ongoing_capex_frac=ONGOING_CAPEX_FRAC)
            return _m if _m != float('inf') else None
        except Exception:
            return None

    # ── Parameter sweep definitions ────────────────────────────────────────────
    # (label, kwarg_key, low_val, high_val, format_fn)
    def _pct(v):  return f"{v*100:.0f}%"
    def _dol(v):  return f"${v:.2f}"
    def _gL(v):   return f"{v:.0f} g/L"
    def _gLhr(v): return f"{v:.1f} g/L/hr"
    def _kta(v):  return f"{v:.0f} kta"

    _base_dsp_yield = dsp['overall_yield']
    sweep_defs = [
        ("Titer",            "titer",              titer * 0.5,                            titer * 2.0,                            _gL),
        ("Prod. rate",       "rate",               rate * 0.5,                             rate * 2.0,                             _gLhr),
        ("Yield fraction",   "yield_fraction",     max(0.05, yield_fraction - 0.30),       min(0.99, yield_fraction + 0.30),       _pct),
        ("DSP recovery",     "dsp_yield_override", max(0.10, _base_dsp_yield - 0.20),      min(0.99, _base_dsp_yield + 0.15),      _pct),
        ("Carbon overflow loss","carbon_to_co2_frac", 0.0,                                 min(0.80, carbon_to_co2_override + 0.30), _pct),
        (f"{carbon_source} price", "price_carbon_per_kg", price_carbon_per_kg * 0.5,       price_carbon_per_kg * 2.0,              _dol),
        ("Scale",            "capacity_kta",       max(0.5, capacity_kta * 0.5),           capacity_kta * 2.0,                     _kta),
    ]

    # ── Compute tornado ────────────────────────────────────────────────────────
    tornado = []
    for label, key, lo, hi, fmt in sweep_defs:
        m_lo = _msp(**{key: lo})
        m_hi = _msp(**{key: hi})
        if m_lo is None or m_hi is None:
            continue
        best  = min(m_lo, m_hi)
        worst = max(m_lo, m_hi)
        # val_best = whichever parameter value gives the lower MSP
        val_best  = lo if m_lo <= m_hi else hi
        val_worst = hi if m_lo <= m_hi else lo
        tornado.append({
            'label':     label,
            'label_ext': f"{label}\n{fmt(lo)} → {fmt(hi)}",
            'key': key, 'fmt': fmt,
            'lo': lo, 'hi': hi,
            'msp_best': best, 'msp_worst': worst,
            'val_best': val_best, 'val_worst': val_worst,
            'swing': worst - best,
        })
    tornado.sort(key=lambda d: d['swing'], reverse=True)

    # ── Tornado plot ───────────────────────────────────────────────────────────
    baseline = MSP
    fig_t, ax_t = plt.subplots(figsize=(9, max(3, len(tornado) * 0.65 + 1.2)))
    for i, d in enumerate(tornado):
        # green bar: improvement region (best → baseline)
        ax_t.barh(i, d['msp_best'] - baseline, left=baseline,
                  color='#2ca02c', alpha=0.75)
        # red bar: deterioration region (baseline → worst)
        ax_t.barh(i, d['msp_worst'] - baseline, left=baseline,
                  color='#d62728', alpha=0.75)
    ax_t.axvline(baseline, color='black', linewidth=1.5,
                 label=f'Current  ${baseline:.2f}/kg')
    ax_t.axvline(target_msp, color='steelblue', linewidth=1.2, linestyle='--',
                 label=f'Target  ${target_msp:.2f}/kg')
    ax_t.set_yticks(range(len(tornado)))
    ax_t.set_yticklabels([d['label_ext'] for d in tornado])
    ax_t.invert_yaxis()
    ax_t.set_xlabel('MSP ($/kg)')
    ax_t.set_title('Sensitivity tornado — one parameter at a time')
    ax_t.legend(loc='lower right')
    plt.tight_layout()
    st.pyplot(fig_t)
    plt.close(fig_t)

    st.caption(
        "Green = improvement from baseline; red = deterioration. "
        "Ranges: titer/rate/glucose price/scale ±50% of current; "
        "yield fraction/DSP yield ±20–30 pp; overflow fraction from 0% to current+30 pp. "
        "Note: parameters interact — the tornado treats each as independent."
    )

    # ── Numeric table ──────────────────────────────────────────────────────────
    st.markdown("**Parameter detail:**")
    sens_rows = []
    for d in tornado:
        # target achievability: linear interpolation
        if d['msp_best'] <= target_msp <= d['msp_worst']:
            frac = (target_msp - d['msp_worst']) / (d['msp_best'] - d['msp_worst'])
            val_needed = d['val_worst'] + frac * (d['val_best'] - d['val_worst'])
            target_str = d['fmt'](val_needed)
        elif target_msp < d['msp_best']:
            target_str = f"< best case ({d['fmt'](d['val_best'])})"
        else:
            target_str = "Not required"
        sens_rows.append({
            "Parameter":      d['label'],
            "Best case":      f"${d['msp_best']:.2f}/kg  @ {d['fmt'](d['val_best'])}",
            "Worst case":     f"${d['msp_worst']:.2f}/kg  @ {d['fmt'](d['val_worst'])}",
            "Swing ($/kg)":   f"${d['swing']:.2f}",
            f"To reach ${target_msp:.2f}/kg": target_str,
        })
    st.dataframe(sens_rows, use_container_width=True)

    with st.expander("What does each parameter mean — and how can it be improved?"):
        st.caption(
            "Parameters are varied one at a time — real improvements often involve "
            "several parameters moving together (e.g. a higher-titre strain also tends "
            "to have lower rate). Use the main sidebar inputs to explore combined scenarios."
        )
        st.markdown("""
| Parameter | What it measures | R&D levers to improve it |
|---|---|---|
| **Titer** | Final product concentration at harvest (g/L). Directly reduces tanks needed and spreads fixed costs over more product per batch. | Strain engineering for higher product tolerance or volumetric titre; fed-batch optimisation; stationary-phase mode (decouple growth from production and target high cell density). |
| **Prod. rate** | Volumetric production rate (g/L/hr) — determines how long the production phase runs. | Stronger or tighter-regulated promoters; improved enzyme kinetics or gene copy number; higher cell density at induction (stationary-phase mode); better dissolved-oxygen control. |
| **Yield fraction** | Fraction of the stoichiometric maximum carbon yield to product actually achieved. Primary driver of glucose (raw material) cost. | Pathway engineering to redirect carbon flux; deletion of competing biosynthetic pathways; cofactor (NADH/NADPH) rebalancing; minimise product degradation during fermentation. |
| **DSP recovery** | Overall fraction of fermentation product recovered after all downstream processing steps. Drives both DSP OPEX (consumables, utilities) and CAPEX (equipment sizing). Use the Advanced DSP expander to see and adjust per-step losses. | Optimise the highest-loss step first: improve extraction selectivity; switch to a secreted product to avoid lysis losses; reduce product degradation during processing; consider a simpler route archetype (e.g. crystallisation vs chromatography). |
| **Carbon overflow loss** | Fraction of non-product glucose diverted to CO₂ and heat via overflow metabolism — acetate in *E. coli*, ethanol in yeast (Crabtree effect) — or cell maintenance, rather than channelled to biomass. **Currently set by the organism preset; use the Advanced biology settings slider to model engineering improvements.** | Knockout overflow pathways (e.g. *pta-ackA* deletion in *E. coli*; replace pyruvate decarboxylase in yeast); use Crabtree-negative strains; operate at lower growth rate; switch to stationary-phase production (removes growth-phase overflow entirely). *Tip: adjust the slider in Advanced settings to model the economic value of a proposed strain improvement before committing lab resource.* |
| **Feedstock price** | Carbon source cost ($/lb) — glucose or methanol — which scales all raw material costs proportionally. | Glucose: bulk supply contracts; crude sugar or molasses; co-location with a sugar producer. Methanol: industrial spot price, highly region-dependent; consider on-site generation or proximity to a petrochemical supply chain. |
| **Scale** | Plant nameplate capacity (kta). Larger plants spread fixed capital and labour over more product — the classic economy-of-scale effect. | Phased capacity expansion; partnership or licensing to reach minimum efficient scale faster; toll manufacturing during early commercialisation. |
""")
