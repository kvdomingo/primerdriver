import { MDBContainer as Container, MDBJumbotron as Jumbotron } from "mdbreact";
import { Suspense, useState } from "react";
import { Route, BrowserRouter as Router, Routes } from "react-router-dom";
import ErrorBoundary from "../utils/ErrorBoundary";
import Characterize from "./menu/Characterize";
import Dna from "./menu/DnaView";
import Menu from "./menu/Menu";
import Protein from "./menu/ProteinView";
import Results from "./menu/Result";
import Loading from "./shared/LoadingScreen";

const styles = {
  appContainer: {
    boxShadow: "none",
    overflow: "auto",
    height: "100vh",
  },
};

function Station() {
  const [key, setKey] = useState(0);

  function handleReset() {
    let newKey;
    do {
      newKey = Math.round(Math.random() * 100);
    } while (newKey === key);
    setKey(newKey);
  }

  return (
    <Container className="px-0">
      <Jumbotron
        className="border border-light mx-0 my-0 py-0 scrollbar scrollbar-primary thin"
        style={styles.appContainer}
      >
        <Router>
          <Suspense fallback={<Loading />}>
            <Routes>
              <Route exact path="/" element={<Menu />} />
              <Route path="/characterize" element={<Characterize key={key} handleReset={handleReset} />} />
              <Route path="/dna" element={<Dna key={key} handleReset={handleReset} />} />
              <Route path="/protein" element={<Protein key={key} handleReset={handleReset} />} />
              <Route
                path="/results"
                element={
                  <ErrorBoundary>
                    <Results />
                  </ErrorBoundary>
                }
              />
            </Routes>
          </Suspense>
        </Router>
      </Jumbotron>
    </Container>
  );
}

export default Station;
