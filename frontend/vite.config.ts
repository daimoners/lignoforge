import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

// The production bundle is written straight into the Python package so that
// `pip install lignoforge` ships a ready-to-serve interface (no Node needed).
export default defineConfig({
  plugins: [react()],
  build: { outDir: "../lignoforge/web/static", emptyOutDir: true, chunkSizeWarningLimit: 1500 },
  server: { port: 5173, proxy: { "/api": "http://127.0.0.1:8765" } },
});
