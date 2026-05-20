# Truncated Fortran Test Fixture

Truncated subset of real simulation data from `data/2026-02-18-1723/`.

## Source

- **Full dataset**: ~665 MB, 50,000 microscale simulations, 10 macroscale simulations
- **Created by**: `scripts/create_test_fixture.py`

## Truncation

- **Microscale**: All files copied as-is (full 50,000 values). Kept intact so
  `generate_macroscale_in()` can reproduce the macroscale_in data.
- **Macroscale input**: All files copied as-is (aggregated data, not per-simulation).
- **Macroscale output**: Only simulation 00 (of 10); truncated to 3 snapshots (of 299):
  - m_loc, m_bound: first 3 x 21,105 int32 values each
  - m_bind_t: all rows with timestamp < tsave[2] (2,127 rows)
  - f_deg_list: all rows with timestamp < tsave[2] (97 rows)
  - tsave: first 3 values
  - Nsave: scalar int32 = 2
  - mfpt: kept as-is (21,105 float64 values)
  - macro log, params.json: copied as-is

## File Codes

- Microscale + macroscale_in: `_PLG2_tPA01_TB-xiii`
- Macroscale output: `_TB-xiii__21_105`

## Required Overrides

When reading this data with `read_data_collection()`, these parameter
overrides are needed (they are not present in the log files):

```
fibrinogen_length=45nm
fibrinogen_radius=1.2nm
micro_log_lvl=40
micro_version=micro_rates
snap_proportion=0.66666667
```

## Size

~4.4 MB total (suitable for git).
