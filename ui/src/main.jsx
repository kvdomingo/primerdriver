import { StrictMode } from "react";
import { createRoot } from "react-dom/client";
import "@fortawesome/fontawesome-free/css/all.min.css";
import "bootstrap-css-only/css/bootstrap.min.css";
import "mdbreact/dist/css/mdb.css";
import "@fontsource-variable/figtree/wght.css";
import "@fontsource-variable/figtree/wght-italic.css";
import App from "./App.jsx";
import { PrimerDriverProvider } from "./contexts/PrimerDriverContext";
import "./index.css";

createRoot(document.getElementById("root")).render(
  <StrictMode>
    <PrimerDriverProvider>
      <App />
    </PrimerDriverProvider>
  </StrictMode>,
);
