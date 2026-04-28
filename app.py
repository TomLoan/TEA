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
    calculate_plant_logistics, calculate_opex, size_equipment, calculate_capex,
    calculate_MSP, calculate_DCF,
    DSP_PRESETS, ORGANISM_PRESETS,
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
    product_type = st.selectbox("Product type", list(DSP_PRESETS.keys()),
                                 format_func=lambda x: x.replace("_", " ").title())
    preset = DSP_PRESETS[product_type]
    dsp_yield = st.number_input("DSP yield (fraction)", value=preset['yield'],
                                 min_value=0.01, max_value=0.99, step=0.01,
                                 help="Product recovered after downstream processing. "
                                      "Small molecules: 0.90. Enzymes: 0.60–0.80. Therapeutics: 0.30–0.60.")
    dsp_opex_frac = st.number_input("DSP OPEX fraction", value=preset['opex_frac'],
                                     min_value=0.01, max_value=0.90, step=0.01,
                                     help="DSP costs as share of total OPEX. "
                                          "Small molecules: 0.20. Enzymes: 0.20–0.40. Therapeutics: 0.50–0.80.")
    dsp_capex_frac = st.number_input("DSP CAPEX fraction", value=preset['capex_frac'],
                                      min_value=0.01, max_value=0.90, step=0.01,
                                      help="DSP capital costs as fraction of total TCI. "
                                           "Small molecules: 0.20. Enzymes: 0.35. Therapeutics: 0.50.")

    # Organism
    st.subheader("Organism")
    organism = st.selectbox("Organism preset", list(ORGANISM_PRESETS.keys()))
    _org = ORGANISM_PRESETS[organism]
    st.caption(_org['note'])
    if organism == 'Mammalian (CHO-like)':
        st.warning("CHO-like cells: this model is aerobic single-substrate only. "
                   "Treat outputs as order-of-magnitude estimates.")
    biomass_yield_override = st.number_input(
        "Biomass yield (gCDW/g glucose)",
        value=float(_org['biomass_yield_coeff']),
        min_value=0.05, max_value=0.80, step=0.01,
        help="Cell mass produced per gram of glucose consumed. "
             "E. coli aerobic: ~0.48. Yeast: ~0.45. Mammalian: ~0.20. "
             "Override the preset if you have measured data."
    )
    _preset_media_cost = _org['media_cost']

    # Fermentation
    st.subheader("Fermentation Performance")
    titer = st.number_input("Titer (g/L)", value=150.0, min_value=1.0, max_value=500.0)
    rate = st.number_input("Rate (g/L/hr)", value=5.0, min_value=0.1, max_value=50.0)
    yield_fraction = st.slider("Yield fraction (% of theoretical)", min_value=1, max_value=99,
                                value=90) / 100.0

    # Plant
    st.subheader("Plant Configuration")
    capacity_kta = st.number_input("Capacity (kta)", value=15.0, min_value=0.1, max_value=500.0)
    tank_volume_L = st.selectbox("Tank volume (L)", [250_000, 500_000, 1_000_000], index=1)
    annual_uptime = st.slider("Annual uptime (%)", min_value=50, max_value=99, value=90) / 100.0
    ferm_temp_C = st.number_input("Fermentation temp (°C)", value=37.0, min_value=4.0, max_value=70.0)
    turnaround_time = st.number_input("Turnaround time (hr)", value=16.0, min_value=1.0, max_value=72.0)
    batches_on_spec = st.slider("Batches on-spec (%)", min_value=50, max_value=100, value=95) / 100.0

    # Raw material prices
    st.subheader("Raw Material Prices")
    price_glucose_per_lb = st.number_input("Glucose ($/lb)", value=0.18, min_value=0.01, max_value=5.0)
    price_ammonia_per_lb = st.number_input("Ammonia ($/lb)", value=0.12, min_value=0.01, max_value=5.0)
    price_mgso4_per_kg = st.number_input("MgSO4 ($/kg)", value=0.30, min_value=0.0, max_value=10.0,
                                          help="Industrial-grade magnesium sulfate (sulfur source). "
                                               "Only relevant for S-containing products.")
    media_cost_per_kgCDW = st.number_input("Media cost ($/kgCDW)", value=_preset_media_cost,
                                            min_value=0.0, max_value=50.0,
                                            help="Mineral salts media. Default $0.40/kgCDW (Lynch 2021). "
                                                 "Typical range: $0.30–$0.80/kgCDW for bacteria; "
                                                 "$5+/kgCDW for mammalian cells.")
    price_NaOH_per_kg = st.number_input("NaOH ($/kg)", value=0.15, min_value=0.01, max_value=5.0)
    price_peracetic_per_L = st.number_input("Peracetic acid ($/L)", value=5.00, min_value=0.1, max_value=50.0)

    # Utility prices
    st.subheader("Utility Prices")
    price_electricity = st.number_input("Electricity ($/kWh)", value=0.11, min_value=0.01, max_value=1.0)
    price_natural_gas = st.number_input("Natural gas ($/MMBtu)", value=3.11, min_value=0.1, max_value=30.0)
    CEPCI = st.number_input("CEPCI", value=603, min_value=100, max_value=1500,
                             help="Equipment cost index for scaling. 2020 = 603 (paper default). "
                                  "2024 ≈ 800. Update for current estimates.")

    # Financial
    st.subheader("Financial Parameters")
    selling_price = st.number_input("Selling price ($/kg)", value=2.50, min_value=0.01, max_value=1000.0)
    target_margin = st.slider("Target margin (%)", min_value=0, max_value=80, value=30) / 100.0
    discount_rate = st.slider("Discount rate (%)", min_value=1, max_value=50, value=20,
                               help="Hurdle rate for NPV. 20% is typical early-stage bioprocess assumption "
                                    "(Lynch 2021).") / 100.0
    tax_rate = st.slider("Tax rate (%)", min_value=0, max_value=50, value=21) / 100.0
    payback_period = st.number_input("Payback period (years)", value=20, min_value=5, max_value=40)
    pct_debt = st.slider("Debt fraction (%)", min_value=0, max_value=90, value=60) / 100.0
    loan_interest = st.slider("Loan interest (%)", min_value=1, max_value=20, value=8) / 100.0
    loan_term_yr = st.number_input("Loan term (years)", value=10, min_value=1, max_value=30)

# ── Derived inputs ────────────────────────────────────────────────────────────
price_glucose_per_g = price_glucose_per_lb / 453.592
price_ammonia_per_g = price_ammonia_per_lb / 453.592
cost_of_fuel = price_natural_gas / 1.05505
capacity_kg = capacity_kta * 1e6

OPEX_KWARGS = dict(
    price_glucose_per_g=price_glucose_per_g,
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
    if chem['atoms']['S'] > 0 and product_type == 'small_molecule':
        st.warning("S detected in formula but Product type is 'Small Molecule'. "
                   "Consider switching to Industrial Enzyme or Therapeutic Protein.")

    ferm = run_fermentation_model(titer, rate, yield_fraction, chem,
                                   biomass_yield_coeff=biomass_yield_override)
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

    # Two-pass OPEX (first pass without other_fixed to get CAPEX, then re-run)
    opex_pass1 = calculate_opex(logistics, ferm, chem,
                                 DSP_yield=dsp_yield, DSP_OPEX_frac=dsp_opex_frac,
                                 **OPEX_KWARGS)
    sizing = size_equipment(logistics, ferm, opex_pass1, tank_volume_L, ferm_temp_C)
    capex = calculate_capex(sizing, dsp_capex_frac)
    other_fixed = 0.037 * capex['TCI_total']
    opex = calculate_opex(logistics, ferm, chem,
                          DSP_yield=dsp_yield, DSP_OPEX_frac=dsp_opex_frac,
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
mcol1.metric("MSP", f"${MSP:.2f}/kg")
mcol2.metric("IRR", f"{dcf['IRR']:.1f}%")
mcol3.metric("NPV (20-yr)", f"${dcf['NPV']/1e6:.1f}M")
mcol4.metric("TCI (total)", f"${capex['TCI_total']/1e6:.0f}M")

with st.expander("Model scope & assumptions"):
    st.markdown("""
**Applicable to:** Aerobic batch/fed-batch fermentation, glucose as sole carbon source,
products containing C/H/O/N/S atoms, scale ~1–100 kta. **Accuracy: ±50% (FEL-1 level).**

**Not modelled:** Anaerobic fermentation · Alternative feedstocks (xylose, methanol, etc.) ·
Stationary-phase production · GMP/regulatory costs · Multi-substrate media · Fed-batch feeding strategies.

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
tab_chem, tab_ferm, tab_opex, tab_capex, tab_fin = st.tabs(
    ["Chemistry", "Fermentation", "OPEX", "CAPEX", "Financials"]
)

# ── Chemistry tab ─────────────────────────────────────────────────────────────
with tab_chem:
    st.subheader("Stoichiometry & Yields")
    eq = chem['equation']
    ccol1, ccol2, ccol3 = st.columns(3)
    ccol1.metric("Molecular weight", f"{chem['MW']:.2f} g/mol")
    ccol2.metric("H₂:CO₂ ratio", f"{eq['ratio']:.3f}")
    ccol3.metric("Redox class", f"Case {eq['case']}")

    ccol4, ccol5 = st.columns(2)
    ccol4.metric("Theoretical yield", f"{chem['theoretical_yield']:.4f} g/g glucose")
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
    fcol1.metric("Fermentation time", f"{ferm['ferm_time']:.1f} hr")
    fcol2.metric("Final biomass", f"{ferm['final_biomass']:.1f} g/L CDW")
    fcol3.metric("Overall yield", f"{ferm['overall_yield']:.4f} g/g glucose")
    fcol4.metric("Cumulative O₂", f"{ferm['cumulative_O2']:.1f} mmol/L")

    fcol5, fcol6, fcol7, fcol8 = st.columns(4)
    fcol5.metric("Max OTR", f"{ferm['max_OTR']:.2f} mmol/L/hr")
    fcol6.metric("Max kLa", f"{ferm['max_kla']:.4f} s⁻¹")
    fcol7.metric("Max cooling", f"{ferm['max_cooling_rate']:.1f} kJ/L/hr")
    fcol8.metric("Tanks required", f"{logistics['n_tanks']}")

    # Sugar partitioning
    st.markdown("**Glucose partitioning (g/L):**")
    sugar_rows = [
        {"Category": "→ Product (stoichiometric)", "g/L": f"{ferm['sugar_to_product']:.2f}"},
        {"Category": "→ Yield losses (CO₂/byproducts)", "g/L": f"{ferm['sugar_to_biomass']:.2f}"},
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
    ocol1.metric("Total OPEX", f"${opex['total_opex']/1e6:.2f}M/yr")
    ocol2.metric("OPEX per kg", f"${opex['opex_per_kg']:.3f}/kg")
    ocol3.metric("Fermentation OPEX", f"${opex['ferm_opex']/1e6:.2f}M/yr")

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

# ── CAPEX tab ─────────────────────────────────────────────────────────────────
with tab_capex:
    st.subheader("Capital Cost Breakdown")
    kacol1, kacol2, kacol3 = st.columns(3)
    kacol1.metric("TCI (upstream)", f"${capex['TCI']/1e6:.1f}M")
    kacol2.metric("TCI (incl. DSP)", f"${capex['TCI_total']/1e6:.1f}M")
    kacol3.metric("TCI per kg/yr", f"${capex['TCI_total']/capacity_kg:.2f}/kg")

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
