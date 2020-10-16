import React, { useState, useEffect } from "react";
import "./App.css";
import Landing from "./components/Landing";
import Footer from "./components/Footer";
import Station from "./components/Station";

export default function App() {
  const [program_version, setProgVersion] = useState("");
  const [web_version, setWebVersion] = useState("");

  useEffect(() => {
    fetch("/api/version")
      .then(res => res.json())
      .then(res => {
        setWebVersion(res.web_version);
        setProgVersion(res.program_version);
      });
  }, [program_version, web_version]);

  return (
    <div>
      <Landing />
      <Station id="app" />
      <Footer program_version={program_version} web_version={web_version} />
    </div>
  );
}
