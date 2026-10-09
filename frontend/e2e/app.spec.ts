import { expect, test, type Page } from "@playwright/test";

const PROJECT = "e2e demo";
test.describe.configure({ mode: "serial" });

async function waitBadge(page: Page, text: string, timeout = 120_000) {
  await expect(page.locator(".page-head .badge", { hasText: text })).toBeVisible({ timeout });
}

test("home page, branding and theme toggle", async ({ page }) => {
  await page.goto("/");
  await expect(page.getByRole("heading", { name: "Projects", level: 1 })).toBeVisible();
  await expect(page.getByAltText(/LignoForge/)).toBeVisible();
  await page.getByRole("button", { name: /Toggle colour theme/ }).click();
  await expect(page.locator("html")).toHaveAttribute("data-theme", /dark|light/);
  const first = await page.locator("html").getAttribute("data-theme");
  await page.getByRole("button", { name: /Toggle colour theme/ }).click();
  expect(await page.locator("html").getAttribute("data-theme")).not.toBe(first);
});

test("system page reports GROMACS and model status", async ({ page }) => {
  await page.goto("/system");
  await expect(page.getByRole("heading", { name: "System", level: 1 })).toBeVisible();
  await expect(page.getByText("provisional").first()).toBeVisible();
  await expect(page.getByText(/GROMACS/).first()).toBeVisible();
});

test("create a project", async ({ page }) => {
  await page.goto("/");
  await page.getByLabel("Name").fill(PROJECT);
  await page.getByRole("button", { name: "Create project" }).click();
  await expect(page).toHaveURL(/\/p\/e2e%20demo$/);
  await expect(page.getByRole("heading", { name: PROJECT })).toBeVisible();
});

test("build a chain, inspect structure, force field, charges, linkages and files", async ({ page }) => {
  await page.goto(`/p/${encodeURIComponent(PROJECT)}/chains/new`);
  await expect(page.getByRole("heading", { name: "Build a chain" })).toBeVisible();
  // live preview of the resolved specification
  await expect(page.getByText(/monomers · seed/)).toBeVisible();
  await page.getByLabel("Name (optional)").fill("test-5mer");
  await page.getByLabel("Number of monomers").fill("5");
  await page.getByLabel("H weight value").fill("0");
  await page.getByRole("button", { name: "Build chain" }).click();
  await expect(page).toHaveURL(/chains\/c0001/);
  await waitBadge(page, "ready");

  // structure: WebGL canvas and a legend
  const viewer = page.getByTestId("viewer3d");
  await expect(viewer.locator("canvas")).toBeVisible({ timeout: 20_000 });
  await expect(viewer.getByText(/Guaiacyl|Syringyl/).first()).toBeVisible();
  await page.getByRole("button", { name: "OPLS type" }).click();
  await expect(viewer.getByText("opls_145")).toBeVisible();
  await page.getByRole("button", { name: "Charge" }).click();
  await expect(viewer.getByText("+0.7 e")).toBeVisible();
  await page.getByRole("button", { name: "CG beads" }).click();
  await expect(viewer.getByText(/Syringyl/)).toBeVisible();
  await page.screenshot({ path: "test-results/structure.png" });

  // force field tab
  await page.getByRole("link", { name: "Force field" }).click();
  await expect(page.getByText("OPLS-AA types in this chain")).toBeVisible();
  await expect(page.getByText("aromatic C–H").first()).toBeVisible();
  await page.getByLabel("Filter atoms").fill("CB");
  await expect(page.locator("tbody tr").last()).toBeVisible();
  await page.screenshot({ path: "test-results/forcefield.png", fullPage: true });

  // charges tab
  await page.getByRole("link", { name: "Charges" }).click();
  await expect(page.getByText(/How charges are set/)).toBeVisible();
  await expect(page.getByText(/Total chain charge/)).toContainText("+0.0000");
  await page.screenshot({ path: "test-results/charges.png", fullPage: true });

  // linkages tab
  await page.getByRole("link", { name: "Linkages" }).click();
  await expect(page.getByRole("img", { name: "Monomer connectivity graph" })).toBeVisible();

  // files tab
  await page.getByRole("link", { name: "Files" }).click();
  await page.getByRole("button", { name: "Generate files" }).click();
  await expect(page.getByText("lig.top")).toBeVisible();
  const [dl] = await Promise.all([page.waitForEvent("download"), page.getByRole("row", { name: /lig\.top/ }).getByRole("link", { name: "Download" }).click()]);
  expect(dl.suggestedFilename()).toBe("lig.top");
});

test("run an atomistic simulation and analyse it", async ({ page }) => {
  await page.goto(`/p/${encodeURIComponent(PROJECT)}/runs/new?chain=c0001`);
  await expect(page.getByRole("heading", { name: "New simulation" })).toBeVisible();
  await page.locator("select").filter({ hasText: "Water (TIP3P)" }).selectOption("none");
  await page.getByLabel("NVT equilibration (ps)").fill("1");
  await page.getByLabel("Production (ns)").fill("0.004");
  await page.getByLabel("Frame interval (ps)").fill("0.5");
  await page.getByLabel("CPU threads").fill("2");
  await page.getByRole("button", { name: "Start simulation" }).click();
  await expect(page).toHaveURL(/runs\/r0001/);
  await waitBadge(page, "done", 180_000);
  await expect(page.getByLabel("Simulation log")).toContainText("GROMACS");
  await expect(page.getByText("prod.xtc")).toBeVisible();

  // analysis: radius of gyration
  await page.getByRole("link", { name: "Analyse" }).click();
  await page.getByLabel("Run", { exact: true }).selectOption("r0001");
  await page.getByRole("button", { name: "Run analysis" }).click();
  await waitBadge(page, "done", 60_000);
  await expect(page.getByRole("img", { name: /Rg .* versus Time/ })).toBeVisible();
  await expect(page.getByText("Summary statistics")).toBeVisible();
  await page.screenshot({ path: "test-results/analysis-rg.png", fullPage: true });

  // analysis: contacts heat map
  await page.goto(`/p/${encodeURIComponent(PROJECT)}/analysis/new?run=r0001`);
  await page.getByLabel("Observable").selectOption("contacts");
  await page.getByRole("button", { name: "Run analysis" }).click();
  await waitBadge(page, "done", 60_000);
  await expect(page.getByRole("img", { name: "Inter-monomer contact map" })).toBeVisible();

  // an invalid analysis reports a clear error instead of crashing
  await page.goto(`/p/${encodeURIComponent(PROJECT)}/analysis/new?run=r0001`);
  await page.getByLabel("Observable").selectOption("rdf");
  await page.getByRole("button", { name: "Run analysis" }).click();
  await waitBadge(page, "failed", 60_000);
  await expect(page.getByRole("alert")).toContainText(/water/i);
});

test("dark theme renders the chain page legibly", async ({ page }) => {
  await page.goto("/");
  if ((await page.locator("html").getAttribute("data-theme")) !== "dark")
    await page.getByRole("button", { name: /Toggle colour theme/ }).click();
  await page.goto(`/p/${encodeURIComponent(PROJECT)}/chains/c0001`);
  await expect(page.getByTestId("viewer3d").locator("canvas")).toBeVisible({ timeout: 20_000 });
  await page.screenshot({ path: "test-results/dark-structure.png" });
  await page.getByRole("button", { name: /Toggle colour theme/ }).click();
});

test("lists reflect the project state", async ({ page }) => {
  const base = `/p/${encodeURIComponent(PROJECT)}`;
  await page.goto(`${base}/chains`);
  await expect(page.getByRole("link", { name: "test-5mer" })).toBeVisible();
  await page.goto(`${base}/runs`);
  await expect(page.getByText("Atomistic MD")).toBeVisible();
  await page.goto(`${base}/analysis`);
  await expect(page.getByText("Radius of gyration")).toBeVisible();
  await page.goto(base);
  await expect(page.getByText("Recent chains")).toBeVisible();
});
