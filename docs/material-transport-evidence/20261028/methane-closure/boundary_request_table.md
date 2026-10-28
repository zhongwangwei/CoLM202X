# 15 missing CRUJRA boundary files (production reader requests)

Source: audit job codex-ad6b48e0cb0044a9bbf0ce2903195234, `boundary-requests.csv` (23,407 request rows) reduced with the same logic as `audit.py`.

| # | variable | file_year | record (1-based) | bound | model time requested | stencil role |
|---|----------|-----------|------------------|-------|----------------------|--------------|
| 1 | dlwrf | 1994 | 1460 | LB | 1995 day 1 sec 0 | LB (first step 1995-01-01 00:00) |
| 2 | dlwrf | 1997 | 1 | UB | 1996 day 366 sec 75600 | UB (last step 1996-12-31) |
| 3 | dswrf | 1997 | 1 | UB | 1996 day 366 sec 64800 | UB (last step 1996-12-31) |
| 4 | pre | 1994 | 1460 | LB | 1995 day 1 sec 0 | LB (first step 1995-01-01 00:00) |
| 5 | pre | 1997 | 1 | UB | 1996 day 366 sec 75600 | UB (last step 1996-12-31) |
| 6 | pres | 1994 | 1460 | LB | 1995 day 1 sec 0 | LB (first step 1995-01-01 00:00) |
| 7 | pres | 1997 | 1 | UB | 1996 day 366 sec 75600 | UB (last step 1996-12-31) |
| 8 | spfh | 1994 | 1460 | LB | 1995 day 1 sec 0 | LB (first step 1995-01-01 00:00) |
| 9 | spfh | 1997 | 1 | UB | 1996 day 366 sec 75600 | UB (last step 1996-12-31) |
| 10 | tmp | 1994 | 1460 | LB | 1995 day 1 sec 0 | LB (first step 1995-01-01 00:00) |
| 11 | tmp | 1997 | 1 | UB | 1996 day 366 sec 75600 | UB (last step 1996-12-31) |
| 12 | ugrd | 1994 | 1460 | LB | 1995 day 1 sec 0 | LB (first step 1995-01-01 00:00) |
| 13 | ugrd | 1997 | 1 | UB | 1996 day 366 sec 75600 | UB (last step 1996-12-31) |
| 14 | vgrd | 1994 | 1460 | LB | 1995 day 1 sec 0 | LB (first step 1995-01-01 00:00) |
| 15 | vgrd | 1997 | 1 | UB | 1996 day 366 sec 75600 | UB (last step 1996-12-31) |

Record-time meaning: record 1460 of a 365-day, 6-hourly CRUJRA v2.5 file is the last record (1994-12-31 18:00 raw time); record 1 is 1997-01-01 00:00 raw time.
The model reads only one record per missing file but addresses it by absolute index, so full-year files are the safe fix.
