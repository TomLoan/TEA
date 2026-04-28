# Bioprocess TEA Calculator

A browser-based tool for early-stage techno-economic analysis of aerobic fermentation processes. Enter a product formula, fermentation targets, and cost inputs — get MSP, IRR, and a full OPEX/CAPEX breakdown in return.

Based on the FEL-1 model from [Lynch et al. 2021](https://doi.org/10.1016/j.ymben.2021.03.004). Accuracy is ±50% — intended for R&D goal-setting, not detailed engineering.

---

## Requirements

- Python ≥ 3.12
- Git

---

## Setup & run

```bash
git clone https://github.com/TomLoan/TEA.git
cd TEA
pip install -r requirements.txt
streamlit run app.py
```

The app opens automatically at `http://localhost:8501`.

If you use [uv](https://docs.astral.sh/uv/):

```bash
git clone https://github.com/TomLoan/TEA.git
cd TEA
uv pip install -r requirements.txt
streamlit run app.py
```

---

## Files

| File | Purpose |
|---|---|
| `app.py` | Streamlit UI — run this |
| `tea_functions.py` | All calculation functions (importable) |
| `bioprocess_tea_calculator.ipynb` | Source notebook with full derivations and validation against Lynch 2021 |
| `requirements.txt` | Python dependencies |
