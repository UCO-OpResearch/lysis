"""RandomDrawSource protocol + NativeDrawSource / FortranDrawSource.

STATUS: placeholder -- this is the first extraction to actually do,
since nearly every strategy depends on it and it's the most
tangled duplicate_fortran cluster in macroscale.py.

Source material in macroscale.py: every site with
`self.random_numbers[RandomDraw.X]` vs `self.rng.random(...)` /
`self.rng.integers(...)`, plus the KissRandomGenerator vs default_rng
construction near __init__.
"""
