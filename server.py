#!/usr/bin/env python3
"""Tiny no-cache static dev server (stdlib only).

Usage: python3 server.py [port]
Serves the repository root. Open http://localhost:8000/ for the app and
http://localhost:8000/tests/test.html for the in-browser tests.
"""
import http.server
import os
import sys


class NoCacheHandler(http.server.SimpleHTTPRequestHandler):
    extensions_map = {
        **http.server.SimpleHTTPRequestHandler.extensions_map,
        ".js": "text/javascript",
        ".mjs": "text/javascript",
        ".css": "text/css",
        ".html": "text/html",
        ".json": "application/json",
    }

    def end_headers(self):
        self.send_header("Cache-Control", "no-store, must-revalidate")
        super().end_headers()


def main():
    port = int(sys.argv[1]) if len(sys.argv) > 1 else 8000
    root = os.path.dirname(os.path.abspath(__file__))
    handler = lambda *a, **kw: NoCacheHandler(*a, directory=root, **kw)
    with http.server.ThreadingHTTPServer(("", port), handler) as httpd:
        print(f"Serving {root} at http://localhost:{port}/")
        try:
            httpd.serve_forever()
        except KeyboardInterrupt:
            pass


if __name__ == "__main__":
    main()
