const http = require("node:http");
const fs = require("node:fs");
const path = require("node:path");

const root = path.resolve(__dirname, "..");
const types = {
  ".html": "text/html; charset=utf-8",
  ".js": "text/javascript; charset=utf-8",
  ".css": "text/css; charset=utf-8",
};
const port = Number(process.env.PORT || 8080);
const server = http.createServer((request, response) => {
  let relative;
  try {
    relative =
      decodeURIComponent(
        new URL(request.url, "http://localhost").pathname,
      ).replace(/^\/+/, "") || "index.html";
  } catch {
    response.writeHead(400).end();
    return;
  }
  const filename = path.resolve(root, relative);
  const within = path.relative(root, filename);
  if (
    within.startsWith("..") ||
    path.isAbsolute(within) ||
    within.split(/[\\/]/).some((part) => part.startsWith(".")) ||
    !types[path.extname(filename)]
  ) {
    response.writeHead(404).end();
    return;
  }
  fs.readFile(filename, (error, data) => {
    if (error) response.writeHead(404).end();
    else
      response
        .writeHead(200, {
          "Content-Type": types[path.extname(filename)],
          "Cache-Control": "no-store",
        })
        .end(data);
  });
});
server.listen(port, "127.0.0.1", () =>
  console.log(`Orbit: http://localhost:${server.address().port}`),
);
