# PC binary figure text handoff

Status: complete
Created: 2026-07-05T17:43:27+08:00

required_found: 4
required_total: 4
optional_found: 0

This package converts PNG figures into base64 text chunks so the GitHub connector can upload them as text files.

After this directory is uploaded to GitHub, run:

```bash
python3 reassemble_binary_figures.py
```

from this handoff directory to reconstruct the PNG files under `reassembled/`.
