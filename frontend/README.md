# LignoForge web interface

React + TypeScript (Vite). The production bundle is written to
`../lignoforge/web/static` and served by `lignoforge gui`; end users never need Node.

```bash
npm ci
npm run dev          # http://localhost:5173, proxies /api to 127.0.0.1:8765
lignoforge gui --dev --no-browser --port 8765     # in another terminal (enables CORS)
npm run build        # type-check + bundle into the Python package
npx playwright install chromium && npm run test:e2e   # browser tests against the real backend
```

* `src/api/` typed REST client mirroring `lignoforge/web/app.py`
* `src/pages/` one file per screen; `src/components/` shared UI, charts, 3-D viewer
* Design tokens (colours from the logo, light/dark) are in `src/styles.css`
