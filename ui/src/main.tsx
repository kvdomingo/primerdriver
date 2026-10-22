import { StrictMode } from "react";
import { createRoot } from "react-dom/client";
import "@fortawesome/fontawesome-free/css/all.min.css";
import "bootstrap-css-only/css/bootstrap.min.css";
import "mdbreact/dist/css/mdb.css";
import "@fontsource-variable/figtree/wght.css";
import "@fontsource-variable/figtree/wght-italic.css";
import App from "./App";
import "./index.css";

// biome-ignore lint/style/noNonNullAssertion: #root is in index.html
createRoot(document.getElementById("root")!).render(
  <StrictMode>
    <App />
  </StrictMode>,
);
