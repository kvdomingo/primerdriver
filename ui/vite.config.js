import path from "node:path";
import { tanstackRouter } from "@tanstack/router-plugin/vite";
import react from "@vitejs/plugin-react-swc";
import { defineConfig } from "vite";

// https://vite.dev/config/
export default defineConfig({
  plugins: [
    tanstackRouter({ target: "react", generatedRouteTree: "./src/route-tree.gen.ts" }),
    react(),
  ],
  server: {
    host: "0.0.0.0",
    port: 3000,
    proxy: {
      "/api": {
        secure: false,
        changeOrigin: false,
        target: "http://backend:5000",
      },
    },
  },
  resolve: {
    alias: {
      "@": path.resolve(__dirname, "./src"),
    },
  },
});
