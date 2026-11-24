# blend-ceiling data

| File | What | Made by |
|---|---|---|
| `toi_2026-11-24.csv` | ExoFOP-TESS TOI table, full CSV. Downloaded 2026-11-24 from `https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv` (public, no sign-in); 8,175 TOI rows. Not edited. | curl |
| `lc_search.csv` | MAST light-curve search results (SPOC and QLP) per TIC: sector, author, exptime. Metadata only. sector -1 = nothing found, -2 = search failed. | `build_sample.py` |
| `sample.csv` | Every TOI that passes the label and period cuts. A blank `drop_reason` means it is in the sample. | `build_sample.py` |
| `windows.csv` | Each usable predicted transit time for the kept TOIs: sector, light-curve author, T (BTJD), sigma_T (h). | `build_sample.py` |
| `blend_comment_matches.csv` | Every comment-rule match, with its segment and any negation that removed it. | `tag_blends.py` |
| `tic_gaia.csv`, `gaia_dr3_target.csv`, `gaia_dr3_neighbours.csv` | TIC v8 to Gaia DR3 target match, and the DR3 sources within 42″ of each target. | `tag_blends.py` |
| `blends.csv` | The kept sample with the blend tags `blend_C`, `blend_G` and `blend_CG`. | `tag_blends.py` |

No light curves are stored in git.
