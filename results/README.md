# results/

Captured stdout/stderr from every measurement script in `code/`.
Filenames pair with the scripts: `code/<name>.py` writes
`results/<name>_result.txt` and `<name>_err.txt`.

Scripts are run from `code/`, so the redirect is one level up:

```bash
cd code
PYTHONIOENCODING=utf-8 python <name>.py > ../results/<name>_result.txt 2> ../results/<name>_err.txt
```

`PYTHONIOENCODING=utf-8` is required on Windows; without it a redirected
stdout falls back to cp949 and non-ASCII output is mangled.
