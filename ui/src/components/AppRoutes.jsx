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

function AppRoutes() {
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
        className="scrollbar scrollbar-primary thin mx-0 my-0 border border-light py-0"
        style={styles.appContainer}
      >
        <Router>
          <Suspense fallback={<Loading />}>
            <Routes>
              <Route path="/" element={<Menu />} />
              <Route
                path="/characterize"
                element={<Characterize key={key} handleReset={handleReset} />}
              />
              <Route
                path="/dna"
                element={<Dna key={key} handleReset={handleReset} />}
              />
              <Route
                path="/protein"
                element={<Protein key={key} handleReset={handleReset} />}
              />
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

export default AppRoutes;
