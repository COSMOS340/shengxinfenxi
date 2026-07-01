# Literature Data Downloads

Use this reference when a reproduction project needs to cache article files,
supplementary tables, supplementary figures, or public repository files before
analysis starts.

## PMC Supplementary Files With CloudPMC Proof-of-Work

Some public PMC article assets return a short proof-of-work HTML page before
the real PDF or XLSX file. The reusable script is:

```bash
python3 bioinformatics-reproduction/scripts/download_pmc_supplementary_pow.py \
  --url-file pmc_supplementary_urls.tsv \
  --output-dir 03_data_raw/paper/supplementary \
  --manifest 99_logs/pmc_supplementary_download_manifest.tsv
```

The input TSV must contain:

```text
url	filename
https://pmc.ncbi.nlm.nih.gov/articles/instance/.../bin/file.pdf	file.pdf
```

The `filename` column is optional. If it is absent or blank, the script derives
the file name from the URL path.

## What The Script Does

1. Requests the public PMC file URL with a browser-like user agent.
2. Accepts the payload directly when the response is already a real file.
3. Parses `POW_CHALLENGE`, `POW_DIFFICULTY`, and `POW_COOKIE_NAME` from the
   CloudPMC proof-of-work page when present.
4. Solves the SHA-256 nonce challenge and retries the same URL with the required
   cookie.
5. Validates common file signatures, including `%PDF` for PDF and `PK` for XLSX.
6. Writes a TSV manifest with URL, output path, status, byte size, SHA-256 hash,
   proof-of-work difficulty, nonce, timestamps, and error text.

## Boundaries

- Use only for public PMC file URLs that the user is allowed to access.
- Do not use it for publisher paywalls, institutional authentication, private
  repositories, or access controls outside public PMC article assets.
- Keep downloaded files in the project data cache, not inside manuscript notes.
- Record every downloaded file in `00_metadata/data_manifest.tsv` or a project
  download manifest before using it as figure or table evidence.

## QC Checks

- For each PDF, verify the file starts with `%PDF`.
- For each XLSX, verify the file starts with `PK`.
- Keep blocked HTML files when validation fails; they are evidence for the
  failure mode.
- Record SHA-256 hashes so the same supplementary files can be checked later.
