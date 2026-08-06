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
| `tea_functions.py` | All calculation functions (importable, no Streamlit dependency) |
| `dryrun.py` | Headless smoke test — runs the full pipeline at default inputs |
| `requirements.txt` | Python dependencies |

Derivations and the sourcing of each constant are documented in the docstrings and
section comments of `tea_functions.py`.

---

## Checking a change

`dryrun.py` runs the whole calculation chain outside Streamlit and compares the
headline outputs against a stored baseline:

```bash
python dryrun.py
```

It exits non-zero if a result moves by more than 0.1%. When a model change is
*meant* to move the numbers, update the `EXPECTED` dict at the bottom of the file.
