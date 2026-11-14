const { readdirSync } = require("node:fs");
const { spawnSync } = require("node:child_process");
const path = require("node:path");

const root = path.resolve(__dirname, "..");
const suites = readdirSync(root)
  .filter(
    (name) => /^test-.*\.cjs$/.test(name) && name !== "test-gpu-runner.cjs",
  )
  .sort();
let failures = 0;
for (const suite of suites) {
  const result = spawnSync(process.execPath, [suite], {
    cwd: root,
    encoding: "utf8",
    timeout: 120_000,
  });
  if (result.error || result.status !== 0) {
    failures++;
    process.stderr.write(
      `FAIL ${suite}\n${result.error || ""}\n${result.stdout}${result.stderr}`,
    );
  } else {
    console.log(`PASS ${suite}`);
  }
}
console.log(`${suites.length - failures}/${suites.length} suites passed.`);
process.exitCode = failures ? 1 : 0;
