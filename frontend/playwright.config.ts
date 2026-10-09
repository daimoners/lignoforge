import { defineConfig } from "@playwright/test";
import { mkdtempSync } from "node:fs";
import { tmpdir } from "node:os";
import { join } from "node:path";

// Starts the real backend (with the built bundle) in a throw-away workspace.
const workspace = mkdtempSync(join(tmpdir(), "lignoforge-e2e-"));
const PORT = 8798;

export default defineConfig({
  testDir: "./e2e",
  timeout: 240_000,
  workers: 1,
  reporter: [["list"]],
  use: {
    baseURL: `http://127.0.0.1:${PORT}`,
    viewport: { width: 1440, height: 900 },
    launchOptions: { args: ["--use-gl=angle", "--use-angle=swiftshader", "--enable-unsafe-swiftshader", "--ignore-gpu-blocklist"] },
    trace: "retain-on-failure",
    screenshot: "only-on-failure",
  },
  webServer: {
    command: `python -m lignoforge.cli.main gui --workspace "${workspace}" --no-browser --port ${PORT}`,
    cwd: "..",
    env: { PYTHONPATH: ".." },
    url: `http://127.0.0.1:${PORT}/api/v1/system`,
    reuseExistingServer: false,
    timeout: 60_000,
  },
});
