# Python TEA vs the original Lynch 2021 calculator - measured comparison

Date: 2026-09-25. Python model at commit 5a66664 (`tea_functions.py`). Original: github.com/mdlynch3/bioprocess-tea-calculator, commit 8d0a36e ("Release 3").

## 1. Summary

1. **MSP: ours $3.17/kg, original $2.38/kg (+33%).** Of the +$0.79/kg gap, $0.60 comes from the two tools defining MSP differently. Our 30% margin is net profit after tax and depreciation. Lynch's is a pre-tax cash margin that includes loan repayments. Another $0.40 is our ongoing capital spend of 10% of TCI per year against Lynch's 1% (`ONGOING_CAPEX_FRAC`, tea_functions.py:293). That looks like a transcription slip: app.py:481 attributes the 10% to "Lynch 2021". Our cheaper cost base pulls the gap back by $0.21/kg.
2. **The original's headline NPV (+$147.7M) is not a real NPV.** BOO_DCF.js divides the 20% discount rate by 100 twice, so it discounts at 0.2%. With 20% applied, the original gives -$19.4M. Its IRR (13.8%) is also pulled down because debt-financed capital is counted twice: once as capital spent and again as loan repayments. Our NPV is +$11.5M and our IRR 29.3%. On the same cost base, our cash-flow method adds about 11 IRR points for the debt fix and 9 for scaling OPEX during ramp-up, and our 10% ongoing capital takes about 10 back.
3. **Real bugs in ours, each worth $0.1-0.25/kg of MSP:**
   - The fermenter cost is not scaled with vessel size. It is priced as a 76 m3 vessel, and fixing it adds $7.6M to TCI.
   - A biomass unit error (/1e6 instead of /1e3) makes media and ammonia costs 1000x too small, which drops $1.0M/yr.
   - The quote sizes for the cooling tower, broth tank and boiler look mis-transcribed, and probably media prep as well. The first three together overstate TCI by about $14M (-$0.43/kg MSP if fixed).
   - Other smaller slips: the tank radius formula, the log base and units of the air-pressure term, and the basis for "other fixed costs".
4. **The original's mass balance double-counts non-product glucose.** It bills 8% more glucose than its own reported yield implies, and it uses 3.2x more O2 than an electron balance allows. Ours matches the electron balance to within 2%, so keep ours.
5. **Plant sizing logic differs.**
   - The original counts turnaround time twice, giving 106 batches per tank per year against our 171. It therefore needs 3 fermenters, and sells 18.25 kt rather than 15 kt.
   - Ours needs 2 fermenters but is internally inconsistent. Utilities, CIP and media are charged on 145 ML of broth, glucose on about 111 ML, and sales on 15 kt.
6. **Completeness check.** Applying the 17 itemised differences to our model one after another reproduces the original's OPEX to within 1.0% and its TCI to within 0.1%. The attribution below therefore accounts for essentially the whole gap.

## 2. How the comparison was done

### 2.1 Running the original (executed, not traced)

The original's calculations run in the browser (`public/js/*.js`), reading inputs from and writing results to the page. I ran them unmodified in headless Chrome:

- `views/calculator.ejs` was rendered to static HTML:
  - The EJS tags (`<%= user... %>`, `<%- include nav2 %>`) were stripped.
  - Script paths were pointed at the repo's `public/js` files.
  - The CDN jQuery was replaced with a 30-line stub covering ready/click/change/remove/append. Chart.js and FontAwesome were removed.
  - `Chart` was stubbed as a no-op, and `window.Color` too (Chart.js 2.x defines that global and the chart code references it).
  - The only load error was a missing image.
- The page's own `$(document).ready` handler calls `UpdateTEA()` exactly as in production. Every numeric global (the model uses implicit globals throughout), the `bioprocessOutputs` and `DCFOutput` objects, and every output span were then dumped.
- Fidelity checks:
  - The chart functions do not feed back into the numbers. `ProFormaChart` rescales `cumCashFlow` only after NPV is computed. The pie-chart functions overwrite some globals with `toFixed` strings, but only after `bioprocessOutputs` is built, so values were read from `bioprocessOutputs`.
  - An independent Python port of BOO_DCF.js + calculateLoanPayments.js reproduces the executed NPV, IRR and MSP to machine precision.
- Caveats:
  - The GitHub code may differ from the version used for the paper or the live site. The executed code gives OPEX $1.44/kg and CAPEX $2.74/kg against the paper's $1.46 and $2.76, so they are close but not identical.
  - `ProcessConstants.js` is not loaded by calculator.ejs; it is dead code with different values.
- Scripts and raw output are in the scratchpad `lynch_compare/` folder:
  - `harness.py` runs the original.
  - `ours_run.py` and `compare.py` build the side-by-side.
  - `attrib.py` does the cost attribution.
  - `fin_attrib.py` and `dcf_bridge.py` do the finance attribution.
  - `js_variants.py` re-runs the original with individual bugs fixed.

### 2.2 Inputs - original (defaults from calculator.ejs and advancedVariables*.js)

| Input | Value | Source |
|---|---|---|
| Product | Diethyl malonate, C7H12O4 (N0) | calculator.ejs:104-116 |
| Selling price | $2.50/kg | calculator.ejs:132 |
| MSP margin | 30% | calculator.ejs:138 |
| Payback period (DCF length) | 20 y, including the 2 construction years | calculator.ejs:142 |
| Discount rate | 20% | calculator.ejs:147 |
| Tax rate | 21% | calculator.ejs:152 |
| Debt financed | 60% | calculator.ejs:156 |
| Debt interest | 8% | calculator.ejs:160 |
| Loan term | 15 y (includes 2 construction years, so 13 y of repayment) | calculator.ejs:165, BOO_DCF.js:128 |
| Capacity target | 15 kt/yr | calculator.ejs:178 |
| Annual uptime | 90% | calculator.ejs:182 |
| Batches on spec | 95% | calculator.ejs:187 |
| Main fermenter | 500,000 L | calculator.ejs:195 |
| Glucose | $0.18/lb ($0.3968/kg) | calculator.ejs:213 |
| Ammonia | $0.12/lb ($0.2646/kg) | calculator.ejs:217 |
| Natural gas | $3.10/MMBtu | calculator.ejs:222 |
| Electricity | $0.11/kWh | calculator.ejs:227 |
| CEPCI | 603 | calculator.ejs:231 |
| Rate / titer / yield | 5 g/L/h, 150 g/L, 90% of theoretical | calculator.ejs:243-252 |
| Turnaround | 16 h | calculator.ejs:257 |
| Media | $0.40/kg CDW | calculator.ejs:262 |
| Temperature | 37 C | calculator.ejs:267 |
| DSP yield / share of OPEX / share of installed cost | 90% / 20% / 20% | calculator.ejs:282-292 |
| Biomass | C3.85H6.69O1.78N, MW 95.37; 0.84 glucose + NH3 + 1.212 O2 -> biomass; achieved biomass yield 80% of that | advancedVariablesProcess.js:11-16 |
| Inoculum, working volume, aspect ratio | 1%, 0.85, 3 | advancedVariablesProcess.js:24-26 |
| NaOH, peracetic acid | $0.15/kg, $5/L | advancedVariablesProcess.js:27-28 |
| Centrifuge | sigma 200,000 m2, settling velocity 6.81e-9 m/s (9,806 L/h), 80% uptime | advancedVariablesProcess.js:33-38 |
| Labour | 1.13 x $203,923 per fermenter per year | advancedVariablesProcess.js:41 |
| Other fixed costs | 1.13 x 3.7% of TDC | advancedVariablesProcess.js:42, UpdateTEA.js:794 |
| Construction | 2 y, 70/30 split | advancedVariablesFinancial.js:5-7 |
| Depreciation | 10 y straight line | advancedVariablesFinancial.js:8 |
| Ongoing capital | 1% of TCI per year | advancedVariablesFinancial.js:9 |
| Ramp-up | 50%, 75%, then 100% | BOO_DCF.js:147-163 |
| Receivable/inventory/payable days | 45/60/30 (computed but never used in cash flow) | advancedVariablesFinancial.js:11-13, BOO_DCF.js:196-213 |
| Capital factors | piping 4.5%, controls 10%, seed 27%, site 9%, warehouse 4%, admin 5%, indirect 70% of TDC, WC 5% of FCI, real estate 6% of FCI | advancedVariablesCapital.js:34, 151-165 |

### 2.3 Inputs - ours (matched)

These are in `ours_run.py` (`LYNCH` dict):

- **Product and organism:** formula C7H12O4, not a protein. Organism "Generic (model default)": biomass yield 0.504 gCDW/g (= 0.8 x 95.37/(0.84 x 180.156)), carbon_to_co2_frac 0, mu_max 0.40 1/h. Growth-associated mode.
- **Fermentation:** titer 150 g/L, rate 5 g/L/h, yield fraction 0.90.
- **Plant:** capacity 15 kt/yr, tank 500,000 L, uptime 0.90, 37 C, turnaround 16 h, on-spec 0.95.
- **Raw materials:** glucose $0.3968/kg, ammonia $0.2646/kg, media $0.40/kgCDW, NaOH $0.15/kg, peracetic $5.00/L, MgSO4 $0.30/kg (unused, since the product has no S).
- **Utilities and cost basis:** electricity $0.11/kWh, natural gas $3.10/MMBtu. CEPCI 603, so the capital scaling ratio is exactly 1. CEPCI also enters the Ulrich and Vasudevan utility-cost formulas in both models, so 603 matches those too.
- **Location and escalation:** location_factor 1.0, escalation 0, years_to_construction 0.
- **Finance:** price $2.50/kg, margin 0.30, discount 0.20, tax 0.21, payback 20, debt 0.60, interest 0.08, loan term 15.
- **Constants left as coded:** ongoing capex 10%, depreciation 10 y, ramp 50/75/100, construction 70/30, other fixed 3.7% of TCI.

**DSP choice.** The original has no unit operations for downstream processing. It sets:

- DSP yield to 90%;
- DSP OPEX so that it is 20% of total OPEX, i.e. 0.25 x fermentation OPEX (UpdateTEA.js:799);
- DSP installed cost to 20% of total installed cost (UpdateTEA.js:766-767). That DSP cost then carries the same site, indirect and working-capital mark-ups as the rest of the plant.

Our route library gives very different numbers for this product:

| Route | DSP share of OPEX | DSP share of TCI | Overall yield |
|---|---|---|---|
| solvent_extraction | 61% | 68% | 0.77 |
| crystallisation | 47% | 58% | 0.86 |
| minimal_processing | 16% | 24% | 0.94 |

At those shares the comparison would be about DSP assumptions, not about the model. I therefore emulated Lynch's rule by passing `calculate_opex`/`calculate_capex` a hand-built `dsp` dict:

- overall_yield 0.90;
- dsp_opex = 0.25 x fermentation OPEX, iterated because "other fixed" depends on TCI;
- dsp_capex = 0.25 x upstream TCI, which puts it through the same mark-ups as in Lynch.

The table also shows the fermentation-only subtotals.

A third column, **"ours @ Lynch logistics"**, forces our model onto Lynch's plant schedule: 3 tanks x 106 batches, 135.15 ML of broth and 18.25 kt sold. This separates formula differences from sizing differences.

## 3. Side-by-side table

Money is in $/yr for OPEX and $ for CAPEX unless stated. "diff" is ours relative to the original.

| Quantity | Original (executed) | Ours | diff | Ours @ Lynch logistics | diff | Classification |
|---|---|---|---|---|---|---|
| Theoretical yield (g/g) | 0.667 | 0.667 | 0.0% | 0.667 | 0.0% | match |
| Product MW (g/mol) | 160.16 | 160.17 | 0.0% | 160.17 | 0.0% | match |
| Glucose used (g/L broth) | 269.9 | 250.0 | -7.4% | 250.0 | -7.4% | ORIGINAL QUESTIONABLE (4.7) |
| Final biomass (gCDW/L) | 12.60 | 12.60 | 0.0% | 12.60 | 0.0% | match |
| Fermentation time (h) | 30.0 | 30.0 | 0.0% | 30.0 | 0.0% | match |
| Cumulative O2 (mmol/L) | 1,032 | 326.6 | -68.4% | 326.6 | -68.4% | ORIGINAL QUESTIONABLE (4.8) |
| Max OTR (mmol/L/h) | 176.2 | 116.4 | -33.9% | 116.4 | -33.9% | BUG (ours) + INTENTIONAL (4.10) |
| Max kLa (1/s) | 0.326 | 0.216 | -33.9% | 0.216 | -33.9% | follows max OTR |
| Number of fermenters | 3 | 2 | -33% | 3 | 0% | ORIGINAL QUESTIONABLE (4.6) |
| Batches per tank per yr | 106 | 171 | +61% | 106 | 0% | ORIGINAL QUESTIONABLE (4.6) |
| Annual broth volume (ML) | 135.15 | 145.35 | +7.5% | 135.15 | 0% | BUG (ours) (4.6) |
| Product sold (kt/yr) | 18.25 | 15.00 | -17.8% | 18.25 | 0% | design difference (4.6) |
| OPEX glucose | 14,473,347 | 11,021,012 | -23.9% | 13,405,408 | -7.4% | ORIGINAL QUESTIONABLE (4.7) |
| OPEX ammonia | 100,546 | 87 | -99.9% | 80 | -99.9% | BUG (ours) (4.11) |
| OPEX media salts | 681,040 | 733 | -99.9% | 681 | -99.9% | BUG (ours) (4.11) |
| OPEX CIP chemicals | 2,067,000 | 1,162,800 | -43.7% | 1,081,200 | -47.7% | UNEXPLAINED (4.12) |
| OPEX process water | 111,785 | 132,157 | +18.2% | 125,466 | +12.2% | ORIGINAL QUESTIONABLE, minor (4.13) |
| OPEX compressed air | 265,644 | 498,633 | +87.7% | 473,321 | +78.2% | mixed (4.8, 4.9) |
| OPEX mass-transfer power | 841,038 | 381,600 | -54.6% | 354,821 | -57.8% | ORIGINAL QUESTIONABLE (4.8) |
| OPEX cooling water | 599,888 | 603,628 | +0.6% | 597,299 | -0.4% | match by coincidence (4.8) |
| OPEX sterilisation | 31,884 | 35,648 | +11.8% | 33,146 | +4.0% | minor (4.13) |
| OPEX heat kill | 234 | 0.28 | -99.9% | 0.26 | -99.9% | BUG (ours), negligible $ (4.11) |
| OPEX centrifugation | 4,548 | 4,797 | +5.5% | 4,460 | -1.9% | volume only |
| OPEX labour | 691,299 | 450,000 | -34.9% | 675,000 | -2.4% | tank count; INTENTIONAL approximation (4.14) |
| OPEX other fixed | 1,109,287 | 1,502,973 | +35.5% | 1,690,113 | +52.4% | BUG (ours) (4.15) |
| OPEX fermentation subtotal | 20,977,541 | 15,794,067 | -24.7% | 18,440,996 | -12.1% | sum |
| OPEX DSP (emulated, see 2.3) | 5,244,385 | 3,948,517 | -24.7% | 4,610,249 | -12.1% | follows subtotal |
| OPEX total | 26,221,926 | 19,742,583 | -24.7% | 23,051,245 | -12.1% | sum |
| OPEX per kg sold ($/kg) | 1.437 | 1.316 | -8.4% | 1.263 | -12.1% | sum |
| Main fermentation area (installed incl. piping) | 6,466,428 | 4,105,424 | -36.5% | 5,350,220 | -17.3% | BUG (ours) mainly fermenter (4.16) |
| Seed train (27% of main area) | 1,745,936 | 1,108,464 | -36.5% | 1,444,559 | -17.3% | follows main area |
| Primary cell removal | 3,745,507 | 3,804,789 | +1.6% | 4,370,341 | +16.7% | offsetting errors (4.17) |
| Utilities + controls | 6,029,665 | 7,373,929 | +22.3% | 7,268,586 | +20.5% | BUG (ours) quote sizes (4.17) |
| TIC, fermentation side | 17,987,536 | 16,392,607 | -8.9% | 18,433,706 | +2.5% | sum |
| TDC excl. DSP | 21,225,292 | 19,343,276 | -8.9% | 21,751,774 | +2.5% | follows TIC |
| Indirect costs excl. DSP | 14,857,704 | 11,605,966 | -21.9% | 13,051,064 | -12.2% | 60% vs 70% of TDC (4.18) |
| FCI excl. DSP | 36,082,997 | 30,949,242 | -14.2% | 34,802,838 | -3.5% | sum |
| Working capital excl. DSP | 1,804,150 | 1,547,462 | -14.2% | 1,740,142 | -3.5% | 5% of FCI in both |
| Real estate excl. DSP | 2,164,980 | 0 | -100% | 0 | -100% | BUG (ours), omission (4.18) |
| TCI excl. DSP | 40,052,126 | 32,496,704 | -18.9% | 36,542,980 | -8.8% | sum |
| DSP share of TCI | 10,013,032 | 8,124,176 | -18.9% | 9,135,745 | -8.8% | emulated as 20% |
| TCI total | 50,065,158 | 40,620,880 | -18.9% | 45,678,725 | -8.8% | sum |
| TCI per kg of annual output ($/kg) | 2.744 | 2.708 | -1.3% | 2.504 | -8.8% | sum |
| MSP ($/kg) | 2.384 | 3.172 | +33.1% | 3.008 | +26.2% | definition + ongoing capex (4.1, 4.2) |
| NPV ($M) | 147.67 | 11.52 | -92.2% | 20.43 | -86.2% | ORIGINAL QUESTIONABLE (4.3) |
| IRR (%) | 13.76 | 29.28 | +15.5 pts | 34.12 | +20.4 pts | mostly ORIGINAL QUESTIONABLE (4.4) |

The original reports "TDC", "FCI" and so on including the 20% DSP share, so the headline figures are TDC $26.53M, indirect $18.57M, FCI $45.10M, WC $2.26M and real estate $2.71M. Our model adds DSP after the rollup, so the rows above compare the non-DSP parts.

ROI is not compared. The original divides the sum of all cash flows by the sum of capital spent (258.7%), while ours divides operating cash flows by TCI (434.6%). They are different definitions.

## 4. Differences above about 2%

Each item's effect on our MSP was measured in `attrib.py`. The change was patched into our code in memory, one item at a time, starting from our matched run: MSP $3.172, OPEX $19.74M, TCI $40.62M, NPV $11.52M, IRR 29.28%. The tables below give "if ours adopted Lynch's version" deltas (dMSP in $/kg, dIRR in points).

### 4.1 MSP definition (largest single MSP item)

- **Ours** (tea_functions.py:1364-1400): MSP is the price at which net income (after tax, depreciation, ongoing capital and interest) is 30% of revenue in the first full-production year: `(OPEX + depreciation + ongoing + interest) x (1-tax) / ((1-tax) - margin)`.
- **Original** (BOO_DCF.js:247-255): `(ongoing capital + OPEX + loan payment) / (production x (1 - margin))` in project year 4, which is also the first full-production year, so the timing matches. Taxes and depreciation are ignored, and principal repayment is treated as a cost.
- **In plain terms:** the original asks "what price leaves 30% of sales after paying the cash bills, including the bank?". Ours asks "what price leaves 30% of sales as profit after tax?". Ours is the more usual definition of a profit margin.
- **Measured** on the original's own cost base (OPEX $26.22M, TCI $50.07M, 18.25 kt), using fin_attrib.py:
  - original formula: $2.384;
  - our formula with 1% ongoing: $2.980 (+$0.596, +25%);
  - our formula with our 10% ongoing: $3.378.
- **Classification:** ORIGINAL QUESTIONABLE (it ignores tax), but this is a definition choice and ours is not documented as different. Keep ours and label it in the app.

### 4.2 Ongoing capital 10% vs 1% of TCI per year

- **Ours:** `ONGOING_CAPEX_FRAC = 0.10` (tea_functions.py:293). It is charged every year in both MSP and DCF (tea_functions.py:1380, 1423, 1450). app.py:481 labels it "10% TCI/yr - Lynch 2021".
- **Original:** `ongoingCapitalReinvestmentRate = 0.01` (advancedVariablesFinancial.js:9, BOO_DCF.js:69).
- **Effect:** switching ours to 1% gives MSP -$0.393/kg, NPV +$11.9M, IRR +9.0 points. Ten per cent of the whole investment re-spent every year is also very high next to the 3.7% "other fixed" charge, which already covers maintenance.
- **Classification:** BUG (ours).

### 4.3 NPV in the original

- **Original:**
  - BOO_DCF.js:16 sets `discountRate = Input[10]/100` (0.20).
  - BOO_DCF.js:240 then discounts by `Math.pow(discountRate/100 + 1, projectYear + 1)`, i.e. at 0.2% per year.
  - The reported $147.7M is therefore almost the undiscounted sum of cash flows ($152.8M).
  - Year 0 is also discounted one year (exponent `projectYear + 1`), whereas ours leaves year 0 undiscounted (tea_functions.py:1462).
- **Measured:**
  - fixing the rate turns the original's NPV into -$16.2M;
  - also using our year-0 convention gives -$19.4M.
  - The JS-variant run (js_variants.py, "npv") confirms -$19.39M.
- **Classification:** ORIGINAL QUESTIONABLE (a bug). The NPV headline gap (-92%) is almost entirely this.

### 4.4 Cash-flow structure and IRR

`dcf_bridge.py` rebuilds both cash-flow methods from one parametrised routine. It reproduces each tool's NPV and IRR exactly. It then switches from the original's method to ours one item at a time, on the original's cost base:

| Step (cumulative) | NPV $M | IRR % | dIRR |
|---|---|---|---|
| Original as coded | 147.67 | 13.75 | |
| Discount at 20%, not 0.2% | -16.16 | 13.75 | 0 |
| Year 0 undiscounted | -19.39 | 13.75 | 0 |
| Loan proceeds netted: only the equity share of capital is a cash outflow | 9.14 | 24.99 | +11.24 |
| Ongoing capital not also depreciated | 8.83 | 24.83 | -0.16 |
| Working capital not depreciated | 8.69 | 24.76 | -0.07 |
| Working capital recovered in final year | 8.76 | 24.78 | +0.02 |
| OPEX scaled with ramp-up output | 20.54 | 33.69 | +8.91 |
| 22-year horizon (2 build + 20 operating) | 21.23 | 33.77 | +0.08 |
| Loan: annual, 15 y from first operating year | 21.79 | 34.31 | +0.54 |
| Ongoing capital 10% | 6.98 | 24.72 | -9.58 |

The last row equals our `calculate_DCF` on the original's cost base.

**Debt counted twice.** The original's net cash flow (BOO_DCF.js:219) subtracts the full capital spend in years 0-1 (BOO_DCF.js:47-71) and also subtracts every loan payment. The 60% borrowed is effectively paid for twice. Ours treats the loan as money in: only 40% of capital is an outflow, and the loan is then repaid (tea_functions.py:1437, 1455), which is the standard equity view. Classification: ORIGINAL QUESTIONABLE. This is worth +11.2 IRR points.

**OPEX during ramp-up.**

- The original charges full OPEX in years 2-3 while producing 50% and 75% (BOO_DCF.js:105-113). That is too harsh: glucose and utilities should fall with output.
- Ours scales all OPEX with output (tea_functions.py:1447), including labour, "other fixed" and DSP. That is too kind: fixed costs do not fall.
- Charging fixed costs in full during ramp-up would cost ours about -$0.76M NPV and -0.76 IRR points.
- Classification: ORIGINAL QUESTIONABLE for the full-OPEX assumption; BUG (ours), small, for scaling fixed costs.

**Smaller items.**

- The original both expenses and depreciates ongoing capital (BOO_DCF.js:89-96 and 171), which double-counts it for tax. ORIGINAL QUESTIONABLE, -0.16 points.
- It depreciates working capital and real estate (BOO_DCF.js:82), and never uses the receivable, inventory and payable days it calculates (BOO_DCF.js:196-213).
- "Payback period 20" means 20 years in total for the original but 2 + 20 in ours (tea_functions.py:1425). Document this.
- "Loan term 15" means 13 years of repayment after construction in the original (BOO_DCF.js:128), but 15 years in ours. Document this too. Using term 13 in ours moves MSP by -$0.004 and IRR by -0.6 points.

### 4.5 Intentional changes, quantified

- **CAPEX follows CEPCI:** zero effect here (CEPCI 603, ratio 1.000). At our default CEPCI 800, capital is x1.327, and the utility-cost formulas also move, as they do in the original.
- **Loan drawn during construction and amortised from the first operating year; working capital not depreciated and recovered; MSP interest taken in the first full-production year** (the 2026-09-24 change). Measured with the pre-change functions from git (5ae7062^) on our matched cost base:
  - MSP $3.159 -> $3.172 (+$0.013);
  - NPV $11.25M -> $11.52M;
  - IRR 28.75% -> 29.28%.
  - On the original's cost base: $3.368 -> $3.378.
  - These are small at a 15-year loan. The "14% never repaid" problem was mostly a short-loan-term effect.
- **MaxOTR uses organism mu_max:** see 4.10. It raises max OTR by 30% relative to the batch-fitted rate and has no cost effect in either model.
- **Location factor 1.0 and escalation 0:** no effect.

### 4.6 Plant logistics: batches, fermenters, broth volume, product sold

- **Original** (UpdateTEA.js):
  - UpdateTEA.js:364 multiplies uptime (7,884 h) by the fermenting fraction 30/46 to get 5,142 h;
  - UpdateTEA.js:365 then divides by the full 46 h cycle and applies on-spec: floor(5,142/46 x 0.95) = 106 batches.
  - Turnaround is thereby subtracted twice; the correct count is floor(7,884/46 x 0.95) = 162.
  - Fermentation output is grossed up for the 90% DSP yield (UpdateTEA.js:359), giving ceil(16.67 kt / 6.76 kt per tank) = 3 tanks.
  - All 3 tanks run full, and the plant sells what they make: 18.25 kt, reported as "Optimal Annual Capacity" (UpdateTEA.js:371).
- **Ours** (tea_functions.py:775-790):
  - floor(7,884/46) = 171 batches;
  - n_tanks = ceil(15/0.95 kt / 10.9 kt per tank) = 2. This is not grossed up for DSP yield, although the answer is still 2 here.
  - Revenue and per-kg figures use the 15 kt nameplate.
  - Utilities, CIP, media and labour are charged on 2 x 171 x 425 m3 = 145.35 ML of broth, which is enough for about 19.6 kt of product.
  - Glucose is charged on 15 kt, equivalent to 111.1 ML of broth, which also ignores off-spec batches (tea_functions.py:819-820).
  - The broth actually needed for 15 kt sold is 116.96 ML.
- **Classification:** ORIGINAL QUESTIONABLE (turnaround double count); BUG (ours) for the inconsistent volumes and the missing DSP gross-up.
- **Effects:**
  - Running ours on the original's logistics: OPEX +$3.31M/yr, TCI +$5.06M, MSP -$0.164/kg (more product sold), IRR +4.8 points.
  - Scaling our volume-driven costs to the 116.96 ML actually needed: OPEX -$0.61M, TCI -$1.97M, MSP -$0.117/kg.
  - Fixing only the original's batch count (js_variants.py): 2 tanks, TCI $46.2M, OPEX $1.413/kg.

### 4.7 Glucose per litre (-7.4%)

- **Original** (UpdateTEA.js:422-434):
  - glucose = product glucose + biomass glucose + "by-product" glucose;
  - "by-product" is set to product glucose / 0.9 - product glucose (line 433). That is **all** the non-product glucose, which already contains the biomass glucose of line 423.
  - Result: 224.9 + 20.0 + 25.0 = 269.9 g/L.
  - Its own reported fermentation yield (UpdateTEA.js:398, 0.600 g/g) implies 249.9 g/L.
- **Ours** (tea_functions.py:570-571, 643, 819-820): 249.95 g/L.
- **Classification:** ORIGINAL QUESTIONABLE (double count).
- **Effect if ours adopted it:** glucose +$1.10M/yr, MSP +$0.118/kg.

### 4.8 Cumulative O2 (-68%) and the utilities driven by it

- **Original** (UpdateTEA.js:445-448):
  - biomass O2 = X / Y(x/O2) / 0.8, giving 200 mmol/L. Dividing by 0.8 also charges O2 for biomass that never forms.
  - by-product O2 = the double-counted 25 g/L of glucose from 4.7, fully burnt at 6 mol O2/mol: 832 mmol/L.
  - Total 1,032 mmol/L.
- **Ours** (tea_functions.py:660-666): biomass 160 + waste 166.5 = 326.6 mmol/L.
- **Electron-balance check:** glucose supplies 33.30 mol e-/L (1.387 mol x 24). The product takes 29.97 (0.937 mol x 32) and biomass 2.05 (0.132 mol x 15.5). That leaves 1.28 mol e-, i.e. **319 mmol O2/L**. Ours is within 2%; the original is 3.2x too high.
- **Classification:** ORIGINAL QUESTIONABLE.
- **Effect if ours adopted it:** OPEX +$2.43M, TCI +$7.69M (bigger air and cooling plant), MSP +$0.460/kg, IRR -11.6 points.
- **Knock-on utility lines:**
  - Mass-transfer power: same formula in both (1.8 x 0.233 kWh per kg air; UpdateTEA.js:473-475, tea_functions.py:881-882). The difference is O2 (x3.16) partly offset by the air utilisation in 4.9.
  - Cooling water: the Ulrich and Vasudevan cooling-water price (UpdateTEA.js:489, tea_functions.py:891) contains a term `0.00003/flow x CEPCI` per m3. Multiplied by annual m3, it becomes a flow-independent $0.00003 x CEPCI per second of operation: $513k/yr for ours (7,884 h) and $335k/yr for the original (5,142 h). It dominates, so a 3.2x difference in O2 shows up as a 0.6% difference in cost. The match is a coincidence.

### 4.9 Compressed-air pricing: utilisation, pressure term, tank geometry

Compressed air is +88% in ours despite 3.2x less O2. Three factors:

1. **Air volume.**
   - Ours bills O2 / (9.375 mol/m3 x 0.75 utilisation) (tea_functions.py:868).
   - The original divides by 9.375 only, and applies the 1/0.75 to the flow used in the unit-price formula but not to the volume billed (UpdateTEA.js:450-453).
   - Ours is right. Classification: ORIGINAL QUESTIONABLE.
   - Adopting the original's version: -$0.28M/yr, MSP -$0.030.
2. **Pressure term.**
   - The original uses ln(1.742), the hydrostatic head at working fill plus 0.25 bar (UpdateTEA.js:375-379, 452).
   - Ours uses log10(32.1 psig) (tea_functions.py:874-878). That is 1.506 against 0.555, a 2.7x higher air price per m3.
   - The correlation, as coded in the original, takes the natural log of pressure in bar. Using a base-10 log of psig mixes both the log base and the units. Classification: BUG (ours).
   - Adopting the original's version: -$0.39M/yr, MSP -$0.042.
3. **Tank geometry.**
   - Ours computes radius = (V/(pi x AR))^(1/3) and height = 2 x AR x r (tea_functions.py:872-873, 1092-1093). The implied volume pi r^2 h = 2V, so the tank is 22.5 m tall instead of 17.9 m.
   - The original: r = (V/(2 pi AR))^(1/3) (UpdateTEA.js:375).
   - Classification: BUG (ours).
   - Fixing it: OPEX -$0.06M, TCI -$0.44M, MSP -$0.018.

### 4.10 Max OTR and kLa (-34%)

- **Original** (UpdateTEA.js:458-470):
  - The peak rate is evaluated at the **batch-fitted logistic rate** k = -ln(0.01/A)/t = 0.307 1/h (UpdateTEA.js:397), not a biological mu_max. The current docstring says the reverse; it is wrong.
  - The biomass term is divided by BiomassYieldFraction = 0.8 and the by-product term by (1 - 0.8) = 0.2 (advancedVariablesProcess.js:16).
- **Ours** (tea_functions.py:671-674, same at 620-621):
  - uses mu_max = 0.40 (INTENTIONAL, +30%);
  - but divides by `_biomass_yield` (0.504 gCDW/g glucose) and (1 - 0.504). That mis-transcribes Lynch's dimensionless 0.8 fraction as the biomass yield coefficient. Classification: BUG (ours).
- **Check:** our formula with 0.8 and k reproduces the original exactly (176.24). With 0.8 and mu_max 0.40 it gives 229.9.
- **Both are inconsistent with cumulative O2.** Lynch's by-product term implies 166 mmol O2 per g of biomass formed. Over the batch that is about 2,080 mmol/L, twice even the original's own cumulative figure. If O2 demand tracks growth, as the growth-associated model assumes, the peak is cumulative O2 x (peak growth rate) / (biomass formed) = 326.6 x (0.307 x 12.6/4) / 12.47 = **25 mmol/L/h** for ours (80 with the original's O2). Both tools report 116-176, so they overstate the aeration challenge 5-7x. Classification: ORIGINAL QUESTIONABLE.
- **No cost effect:** max OTR, kLa and max cooling feed no cost in either model.

### 4.11 Media, ammonia and heat kill (BUG, ours)

- **Biomass units.** tea_functions.py:789 sets `annual_biomass_kg = annual_ferm_vol x final_biomass / 1e6`. L x g/L gives grams, so /1e6 gives tonnes, not kg.
  - Media (tea_functions.py:837) and biomass ammonia (tea_functions.py:827) are therefore 1000x too small: $733 instead of $733k, and $87 instead of about $87k.
  - The original: UpdateTEA.js:432, 440.
  - Fixing it: OPEX +$1.02M/yr, MSP +$0.110/kg, IRR -2.3 points.
  - After fixing, our biomass NH3 is still 20% below the original's, because the original divides by 0.8 again (UpdateTEA.js:432) and charges NH3 for biomass that is never made. Keep ours.
- **Heat kill.** tea_functions.py:897 uses /1e6 where the original uses /1000 (UpdateTEA.js:519). The dollar effect is negligible either way ($234/yr in the original), but the same volume sizes the waste-water tank (tea_functions.py:1142, clamped at 10 m3).

### 4.12 CIP chemicals (-44%)

- **Original** (UpdateTEA.js:500-503): per batch, gross vessel volume x (2% NaOH x $0.15 + 0.002 L/L peracetic x $5) = $6,500. Its own comment says peracetic is at "200 ppm or 0.02%", which the 0.002 does not match.
- **Ours** (tea_functions.py:839-842): working volume x (2% x $0.15 + 0.001 L/L x $5) = $3,400.
- **Classification:** UNEXPLAINED. Neither matches the stated 200 ppm, and ours is not documented.
- **Effect if ours adopted the original's:** +$1.33M/yr, MSP +$0.142. This is a surprisingly large cost line (8% of the original's OPEX) for an assumption with no cited source.

### 4.13 Water, sterilisation, centrifugation (minor)

- **Water:** +12% at the same volume.
  - The original averages water demand over the 5,142 fermenting hours (UpdateTEA.js:417); ours over 7,884 h (tea_functions.py:863).
  - The lower flow in ours raises the correlation's unit price: $0.93 against $0.83/m3.
  - The original's basis inherits the double-count of 4.6. Classification: ORIGINAL QUESTIONABLE, minor.
  - The same hours basis gives most of the small residual in cooling water (4.8).
- **Sterilisation:** ours applies a 1.05 seed-volume allowance against 1.01 (tea_functions.py:894 vs UpdateTEA.js:495), +4%, about $1k/yr.
- **Centrifugation:** our flow is 10,000 L/h against 9,806 and our uptime basis differs (tea_functions.py:901-907 vs advancedVariablesProcess.js:33-38). The count is 2 in both; the cost difference is under $0.3k.

### 4.14 Labour (-35%; -2.4% at the same tank count)

- **Ours:** 2.5 operators x $60k x 1.5 = $225k per fermenter (tea_functions.py:927-932).
- **Original:** 1.13 x $203,923 = $230.4k (advancedVariablesProcess.js:41).
- Nearly all of the headline gap is the tank count (4.6).
- **Classification:** INTENTIONAL (documented as an approximation in the docstring). Effect +$0.001/kg.

### 4.15 Other fixed costs (+36%)

- **Ours:** 3.7% of TCI_total (app.py:397, dryrun.py:147).
- **Original:** 1.13 x 3.7% of TDC (UpdateTEA.js:794). TDC excludes indirect costs, working capital and land; TCI is about 1.9x TDC.
- **Classification:** BUG (ours), wrong basis.
- **Effect of adopting the original's basis:** -$0.62M/yr, MSP -$0.066.

### 4.16 Fermenter purchase cost (-82%)

- **Ours:** EQUIP_DB `'fermenter': (176000, None, None, 1.13, 2.0)` (tea_functions.py:1035), used unscaled at tea_functions.py:1208. Every fermenter costs $398k installed whatever its size.
- **Original:** $176k x 1.12 is the NREL quote for a **76 m3** vessel, scaled by (V/76,000 L)^0.7 (advancedVariablesCapital.js:22-25, UpdateTEA.js:534-544). That gives $737k purchased, $1.47M installed per 500 m3 vessel.
- **Classification:** BUG (ours). With it, our fermenter cost does not respond to the tank-size choice at all.
- **Effect of fixing:** TCI +$7.61M, MSP +$0.235/kg, IRR -7.1 points.

### 4.17 Other equipment items

Installed cost before piping, at matched inputs (2 tanks in ours):

| Item | Original | Ours | diff | Main reason | Classification |
|---|---|---|---|---|---|
| Agitators | 466,280 | 524,786 | +13% | Ours sized on 2 kW/m3 power, $36k at 36 kW (tea_functions.py:1036, 1099); original on vessel volume, $36k x 1.12 at 75.7 m3 (advancedVariablesCapital.js:28-31) | UNEXPLAINED (different basis) |
| Fermenter transfer pumps | 38,000 | 151,244 | +298% | Ours sized on 1 h drain flow; original 2 per tank on vessel volume | UNEXPLAINED, <$0.2M |
| Glucose feed tank + pump | 117,126 | 438,439 | +274% | Ours one tank per fermenter at 2% x 12 h of volume, quote size 100 m3; original one tank for 12 h of glucose, quote size 265 m3 (advancedVariablesCapital.js:38-42) | UNEXPLAINED |
| Base + acid tanks + pumps | 22,816 | 472,396 | +1970% | Ours 0.5% x 12 h of volume per fermenter each; original sized on NH3 consumption, so tiny | UNEXPLAINED |
| Media prep | 639,541 | 1,131,542 | +77% | Ours quote size 100 m3 and inflation 1.17 (tea_functions.py:1045); original 264,978 L and 1.12 (advancedVariablesCapital.js:56-60) | BUG (ours), quote size |
| CIP | 491,095 | 414,709 | -16% | Original `Math.pow(x),SF` typo (UpdateTEA.js:595, 597) makes each tank quote x 0.7 = $77k instead of the scaled $13k; ours quote size 10 m3 against 106 m3 (tea_functions.py:1047) | both wrong; ORIGINAL QUESTIONABLE + BUG (ours) |
| Centrifuges | 1,717,200 | 1,860,300 | +8% | Quote $325k against $300k (tea_functions.py:1049 vs advancedVariablesCapital.js:86) | UNEXPLAINED, small |
| Broth tank + pump | 1,867,017 | 1,780,647 | -5% | Same `Math.pow` typo in the original (UpdateTEA.js:630-632) gives quote x 0.7 = $1.03M (scaled would be $0.42M); ours quote size 1,000 m3 against 4,542 m3 (tea_functions.py:1050) overstates 2.9x. Fixing ours: TCI -$3.07M, MSP -$0.095 | both wrong; BUG (ours) |
| Cooling tower + pumps | 761,283 | 1,978,191 | +160% | Ours quote sizes 0.1 m3/s (tea_functions.py:1052-1053); the NREL quotes are 44,200 and 16,120 gpm = 2.79 and 1.02 m3/s (advancedVariablesCapital.js:96-103). Fixing ours: TCI -$4.82M, MSP -$0.149 | BUG (ours) |
| Boiler | 2,666,910 | 2,708,599 | +2% | Coincidence, see below. Fixing ours: TCI -$5.94M, MSP -$0.183 | BUG (ours) + ORIGINAL QUESTIONABLE |
| Air compressor | 319,885 | 80,587 | -75% | Original sized on scfm (2,280 scfm, about 3.2x our O2 flow); ours on shaft power at $1,000/kW (tea_functions.py:1055, 1113-1121) | follows 4.8; basis UNEXPLAINED |
| Air receiver + dryer | 65,303 | 249,533 | +282% | Different quotes and factors (ours $17k at 10 m3, installation factor 3.1; original $104.6k at 25,000 gal) | UNEXPLAINED, <$0.25M |
| Water handling | 275,237 | 356,692 | +30% | Different sizing bases | UNEXPLAINED, small |
| Waste water | 79,448 | 151,339 | +90% | Different tank size basis (4.11) | UNEXPLAINED, small |
| Heat exchangers | 145,642 | 254,073 | +74% | Ours quote size 10 m2 against 140 m2 | UNEXPLAINED, small |

Why the boiler matches by coincidence:

- The original's `AVC.BoilerPackageCosts_SF = 1,8;` (advancedVariablesCapital.js:107) is a comma typo. It sets the scaling exponent to 1, and the installation factor it reuses (UpdateTEA.js:671) is also 1.
- The original sizes steam without the 0.2 heat-recovery factor used in its own OPEX, giving 166,700 lb/h (UpdateTEA.js:658-668).
- Ours uses 35,500 lb/h but a quote size of 1,000 lb/h against 10,000 lb/h, and installation factor 2.0 against 1.8 (tea_functions.py:1054).
- The errors happen to cancel.

Several of our quote sizes look like unit slips when set against the original's NREL values. The docstring cites Lynch's Table S4.1, which I could not check here. Confirm against the paper before changing them.

### 4.18 Capital rollup

- **Indirect costs.** Ours uses 60% of TDC (tea_functions.py:1263: 10+10+20+10+10). The original uses 70% (advancedVariablesCapital.js:156-160, UpdateTEA.js:779-784). The difference is project contingency, 10% in ours against 20% in the original. "Contingency" is a reserve for the costs an early estimate cannot see yet. Ours matches the NREL biorefinery reports; for a plus or minus 50% early estimate, 20% is the more defensible figure.
  - Classification: BUG (ours), an undocumented deviation.
  - Effect: TCI +$2.54M, MSP +$0.078.
- **Real estate (land).** The original adds 6% of FCI (UpdateTEA.js:787-788). Ours has none.
  - Classification: BUG (ours), an omission.
  - Effect: TCI +$2.32M, MSP +$0.072.
  - The original then depreciates the land, which is not allowed for tax. ORIGINAL QUESTIONABLE, small.
- **Controls.** Ours takes 10% of area 1+3+4 equipment (tea_functions.py:1253). The original takes 10% of all areas including piping and the seed train (UpdateTEA.js:762). Effect TCI +$0.42M, MSP +$0.013. Minor.
- **Seed train, piping, site, warehouse, admin, working capital:** same factors in both.

### 4.19 Original bugs that do not show at the defaults

Re-running the original with individual fixes applied (js_variants.py; copies in the scratchpad, repo untouched):

| Original with... | OPEX $M | TCI $M | Tanks | MSP | NPV $M | IRR % |
|---|---|---|---|---|---|---|
| as coded | 26.22 | 50.07 | 3 | 2.384 | 147.7 | 13.76 |
| `Math.pow` comma typos fixed (CIP, broth, waste-water tanks/pumps) | 26.08 | 44.80 | 3 | 2.338 | 159.7 | 15.74 |
| boiler comma typo fixed | 26.13 | 46.67 | 3 | 2.354 | 155.4 | 15.01 |
| batch double count fixed | 26.27 | 46.17 | 2 | 2.318 | 166.0 | 15.89 |
| glucose/O2 double count fixed | 23.69 | 47.98 | 3 | 2.172 | 187.3 | 17.27 |
| NPV discount fixed | 26.22 | 50.07 | 3 | 2.384 | -19.4 | 13.76 |
| all of the above | 23.43 | 34.79 | 2 | 2.026 | 11.4 | 24.60 |

Also:

- `else if (vesselSize = 500000)` at UpdateTEA.js:536 is an assignment, not a comparison. With the 200 m3 option selected, the original prices 7 tanks as 500 m3 vessels and resets the global vessel size to 500,000 L for all later calculations (verified: TCI $76.7M). Validation against the original at 200 m3 is therefore meaningless.
- CIP filter and heater are bought but left out of the installed total (UpdateTEA.js:599-605).

### 4.20 Cumulative check (completeness)

Applying all the "Lynch" versions to ours one after another (attrib.py) moves:

| Stage | OPEX | TCI | MSP (ours definition) |
|---|---|---|---|
| Start | $19.74M | $40.62M | $3.172 |
| After all 17 items | $26.48M | $50.14M | $3.008 |
| Original | $26.22M | $50.07M | |

The remaining OPEX residual of +1.0% is mostly the hours basis in the cooling and water formulas (4.8, 4.13). The remaining MSP gap, $3.008 against $2.384, is the MSP definition (4.1).

## 5. Recommendations

### 5.1 Fix (BUG, ours) - in order of MSP impact

1. **Ongoing capital: `ONGOING_CAPEX_FRAC = 0.01`** (tea_functions.py:293), and correct app.py:481. If 10% is deliberate, stop attributing it to Lynch. Effect MSP -$0.39, IRR +9.
2. **Fermenter cost:** set EQUIP_DB `'fermenter': (176000, 76.0, 0.70, 1.12, 2.0)` and call `tic('fermenter', sizing['tank_total_vol_m3'], n_units=n)` (tea_functions.py:1035, 1208). Effect MSP +$0.24.
3. **Quote sizes** (confirm against the paper's Table S4.1 first). Cooling tower, broth tank and boiler together: MSP -$0.43, TCI -$13.8M; the others are small.
   - Cooling tower 2.79 m3/s and cooling pump 1.02 m3/s (tea_functions.py:1052-1053).
   - Broth tank 4,542 m3 with inflation 1.12 (1050).
   - Boiler 10,000 lb/h with installation factor 1.8 (1054).
   - Media prep and feed tank 265 m3 (1038, 1045).
   - CIP tank 106 m3 (1047).
   - Heat exchanger 140 m2 (1064).
4. **Biomass units:** `/ 1e6` -> `/ 1e3` at tea_functions.py:789 (and check 897, 1142). Media and NH3 become real. Effect MSP +$0.11.
5. **Plant volume consistency** (tea_functions.py:775-790, 819-820). Either sell what the tanks make (the original's approach), or scale every volume-driven cost to the broth actually needed. Also gross tank sizing up for DSP yield, and charge glucose for off-spec batches. Effect about -$0.12 on MSP for the second option.
6. **Other fixed costs:** use 1.13 x 3.7% of TDC (including DSP), or document why TCI (app.py:397). Effect MSP -$0.07.
7. **Compressed air:** natural log of pressure in bar, not log10 of psig (tea_functions.py:877-878). Effect MSP -$0.04.
8. **Tank radius:** `(V / (2 * math.pi * ASPECT_RATIO)) ** (1/3)` at tea_functions.py:872 and 1092. Effect MSP -$0.02.
9. **Max OTR:** replace `_biomass_yield` with the biomass yield *fraction* (0.8) at tea_functions.py:620-621 and 673-674, and correct the docstrings (lines 26-31, 510-515). Consider replacing the formula with the mass-balance-consistent peak in 4.10.
10. **Ramp-up:** scale only variable costs with ramp output (tea_functions.py:1447). Effect about -0.8 IRR points.
11. **Decide and document** indirect-cost contingency (10% vs Lynch 20%), land (6% of FCI), and CIP dosing. Together these are worth +$0.15 (contingency and land) and +$0.14 (CIP) on MSP if you adopt Lynch's values.

### 5.2 Keep ours (ORIGINAL QUESTIONABLE) - but say so in the app

- Glucose and O2 accounting (4.7, 4.8): ours closes the electron balance.
- Batch count without the turnaround double-count (4.6).
- Air volume including the 75% O2 utilisation (4.9).
- NPV discounted at the stated rate (4.3).
- Equity cash flow with loan proceeds netted (4.4).
- Ongoing capital expensed once; working capital not depreciated and recovered; land not depreciated.
- MSP as an after-tax net margin (4.1). Label this clearly, because it is the biggest single reason our MSP sits above the original's.
- Document that "payback period" and "loan term" are counted from the first operating year in ours but include construction in the original.

## 6. Proposed replacement for "KNOWN DISCREPANCIES vs LYNCH 2021"

```text
KNOWN DISCREPANCIES vs LYNCH 2021
---------------------------------
Measured 2026-09-25 by executing the original JavaScript (github.com/mdlynch3/
bioprocess-tea-calculator, commit 8d0a36e) at its default inputs (diethyl malonate
C7H12O4, 150 g/L, 5 g/L/h, 90% yield, 15 kta, 500 m3, CEPCI 603) and running this
model at matched inputs (Generic organism, CEPCI 603, location 1, escalation 0,
DSP emulated as Lynch's 90% yield / 20% of OPEX / 20% of installed cost).
Full report: docs/lynch_comparison.md.

Headline (original -> this model): OPEX $26.2M -> $19.7M/yr ($1.44 -> $1.32/kg),
TCI $50.1M -> $40.6M, MSP $2.38 -> $3.17/kg, NPV +$147.7M -> +$11.5M,
IRR 13.8% -> 29.3%.

Deliberate differences (keep):
- Mass balance: glucose 250 vs 270 g/L and cumulative O2 327 vs 1032 mmol/L.
  The original counts non-product glucose twice; this model closes the electron
  balance (319 mmol O2/L) to within 2%.
- Batches: floor(uptime / cycle) = 171 per tank; the original subtracts
  turnaround twice (106 per tank) and so needs 3 fermenters instead of 2.
- Air billed at O2 / (9.375 x 0.75); the original omits the 0.75 from the volume.
- Finance: NPV at the stated discount rate (the original divides it by 100
  twice: 0.2%); equity cash flow with loan proceeds netted (the original pays
  for debt-financed capital twice); loan drawn during construction and amortised
  from the first operating year; WC not depreciated and recovered at the end;
  CAPEX scaled by CEPCI / 603. Together worth about +11 IRR points.
- MSP is an after-tax net-income margin; the original's is a pre-tax cash
  margin including loan principal. Worth +$0.60/kg on the same cost base.
- Peak OTR uses organism mu_max (0.40/h). The original uses the batch-fitted
  logistic rate (0.31/h), not a biological mu_max. Neither max OTR feeds a cost.

Open issues found by the comparison (see report section 5):
- ONGOING_CAPEX_FRAC is 0.10; the original uses 0.01 (+$0.39/kg MSP, -9 IRR pts).
- Fermenter cost is not scaled with vessel volume (the $176k quote is for 76 m3).
- annual_biomass_kgCDW is in tonnes (/1e6), so media and NH3 are 1000x too low.
- Quote sizes for cooling tower/pumps, broth tank, boiler, media prep, CIP tank
  and heat exchangers differ from the original's NREL values by 3-28x.
- Max OTR divides by the biomass yield coefficient (0.504 g/g) where the
  original divides by the biomass yield fraction (0.8).
- Other fixed costs use 3.7% of TCI; the original uses 1.13 x 3.7% of TDC.
- Indirect costs 60% of TDC (original 70%: contingency 20%); no land (6% FCI).
- Tank radius formula gives twice the volume; air-pressure term uses log10(psig).
- Volume-driven costs use full tank capacity (145 ML) while glucose and sales
  use 15 kt (about 111 ML of broth).
With all of these switched to the original's choices, this model reproduces the
original's OPEX within 1% and TCI within 0.1%.
```
