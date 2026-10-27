# Display reproduction environment

Tested separately with Python 3.12.14 and the complete pip lock in requirements-display-lock.txt. This is a new display-only verification environment, not a reconstruction of the historical model-fitting environment. No model is loaded or fitted.

Historical analysis evidence records Python 3.9.6, NumPy 1.26.4, pandas 2.3.3, SciPy 1.13.1, scikit-learn 1.6.1, scikit-image 0.24.0, tifffile 2024.8.30, Matplotlib 3.9.4, Pillow 11.3.0, torch 2.8.0 and Transformers 4.48.1. The older environment.yml describes an earlier preprocessing stage and is not a complete model-fitting lock. Historical full-pipeline reconstruction has not been validated.

Arial is used where installed, otherwise Helvetica / DejaVu Sans. Numeric values are identical across fonts; raster/font geometry may differ across operating systems.
