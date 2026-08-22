# Contributing to CARMEN-Q

Contributions are welcome when they sharpen the access model, test a stated theorem, add an independently falsifiable control, or improve reproducibility and usability.

Before opening a pull request, state whether each scientific claim is known, derived here, conjectural, or interpretive, and cite primary sources for claims about prior work. New behaviour should include a focused test. Run the complete local check with:

```bash
python -m pip install -e ".[dev,reproducibility]"
python scripts/run_all.py
python scripts/build_order_pdf.py
python scratch/originality_gate/validate_streaming_parity.py
```

Update `CLAIMS.md` whenever a central evidential status changes. Reports of endpoint-equivalent shortcuts, classical strategies that exceed the bound, or failures of the causal interface are especially valuable: the project is designed to narrow or abandon claims when their anti-shortcut conditions fail.

Code should target Python 3.10 or later, keep the deterministic seed controls intact, and avoid committing local environments or generated build directories. By contributing code you agree that it is distributed under the MIT License; contributions to the manuscript, figures, and data are distributed under CC BY 4.0 unless explicitly noted.
